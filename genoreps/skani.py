import os
import tempfile
import logging
import subprocess
import shutil
import pyskani
import pandas as pd
import gzip


def _read_fasta_bytes(filepath):
    """
    A lightweight, pure-Python FASTA parser.
    Yields sequence contigs as ascii-encoded bytes for pyskani.
    """
    current_seq = []
    
    with open(filepath, 'rb') as f_check:
        is_gzipped = (f_check.read(2) == b'\x1f\x8b')
        
    open_func = gzip.open if is_gzipped else open
    
    with open_func(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    yield "".join(current_seq).encode('ascii')
                    current_seq = []
            else:
                current_seq.append(line)
        # Yield the final contig
        if current_seq:
            yield "".join(current_seq).encode('ascii')

def run_skani(assemblies, args):
    """
    Greedily filters a list of assemblies based on an ANI threshold using pyskani.
    Assumes 'assemblies' is pre-sorted in the desired priority order.
    """
    logging.info("Starting in-memory pyskani dereplication for %d assemblies...", len(assemblies))
    
    threshold = args.distance
    if threshold > 1.0:
        threshold /= 100.0
        logging.debug("Converted threshold to percentage: %.2f%%", threshold)

    keep = []
    db = pyskani.Database()

    # record fasta names and full paths
    fasta_name_to_path = {os.path.basename(fasta): fasta for fasta in assemblies}

    for current_fasta in assemblies:
        fasta_name = os.path.basename(current_fasta)
        
        # Load the genome into memory as raw bytes
        contigs = list(_read_fasta_bytes(current_fasta))
        
        if not contigs:
            logging.warning("Skipping empty or unreadable file: %s", fasta_name)
            continue

        # 1. The first assembly is automatically our first representative
        if not keep:
            keep.append(fasta_name)
            # Use the full path as the reference name so we don't lose it
            db.sketch(fasta_name, *contigs)
            logging.debug("Added first representative: %s", fasta_name)
            continue
        
        # 2. Query the current genome against the EXISTING representatives
        # (This is incredibly fast because the database is tiny)
        hits = db.query(fasta_name, *contigs)
        
        is_redundant = False
        
        # 3. Check if it matches any kept representative above the threshold
        for hit in hits:
            if hit.identity >= threshold:
                is_redundant = True
                logging.debug("Skipping %s: Matches existing representative (ANI: %.2f)", fasta_name, hit.identity)
                break # Early exit, no need to check other hits
        
        # 4. If it's novel, add it to our keep list and sketch it into the database
        if not is_redundant:
            keep.append(fasta_name)
            db.sketch(fasta_name, *contigs)
            logging.debug("New representative found: %s (Total reps so far: %d)", fasta_name, len(keep))

    logging.info("Dereplication complete. Kept %d out of %d assemblies.", len(keep), len(assemblies))
    
    # Return the exact paths in the order they were kept
    return [fasta_name_to_path[name] for name in keep]