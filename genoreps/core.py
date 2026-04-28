import os
import sys
import logging
import pathlib
import pandas as pd

from .mash import run_mash
from .skani import run_skani

def load_metadata(metafile, assemblies):
    logging.info("Loading metadata file...")

    # Load the TSV file
    df = pd.read_csv(metafile, sep='\t', low_memory=False)

    # Fix accession
    df["accession"] = df["accession"].str.extract(r'(GC[AF]_\d{9}\.\d)')

    # Set overall score according to GTDB formula. 
    df['rep_score'] = (   
        1e7 * (df['gtdb_representative']=='t').astype(float) +
        1e6 * df['gtdb_type_designation_ncbi_taxa_sources'].str.contains("LPSN", na=False).astype(float) +
        1e5 * df['gtdb_type_designation_ncbi_taxa'].str.contains("type strain of species", na=False).astype(float) +
        1e4 * df['ncbi_refseq_category'].str.contains("representative genome", na=False).astype(float) +
        1e3 * df['gtdb_type_designation_ncbi_taxa'].str.contains("type strain of subspecies", na=False).astype(float) +
        1e2 * df['ncbi_assembly_level'].str.contains("Complete Genome", na=False).astype(float) +
        df['checkm_completeness'] - 5*df['checkm_contamination'] -
        100 * (df['mimag_high_quality']=='f').astype(float) -
        5 * (df['contig_count']/1e2) -
        5 * (df['ambiguous_bases']/1e4) +
        10 * (df['ssu_count']>0).astype(float)
    )

    df = df.sort_values(by='rep_score', ascending=False)
    df['gtdb_species'] = df["gtdb_taxonomy"].str.extract(r's__(.*)')
    df['filename'] = df['accession'].map(assemblies)
    df = df.dropna(subset=['filename'])

    return {species: data for species, data in df.groupby('gtdb_species')}

def dereplication(assemblies, metadata, args):
    if getattr(args, 'gtdb_reps', False):
        # Only output the GTDB representatives
        logging.info("Filtering to GTDB representatives...")
        metadata = pd.concat(metadata.values(), ignore_index=True)
        subset_df = metadata[metadata['gtdb_representative'] == 't']
        logging.info(f"Filtered database size: {subset_df.shape[0]}")
        subset_df.to_csv(os.path.join(args.out_dir, "filtered_gtdb_meta.tsv"), sep='\t', index=False)
    else:    
        # Run dereplication
        reduced_file_set = []
        for species in metadata:
            logging.info(f"Processing species: {species}")
            # reduced_file_set += run_skani(list(metadata[species]['filename']), args)
            if len(list(metadata[species]['filename'])) > 1:
                if args.method == 'skani':
                    reduced_file_set += run_skani(list(metadata[species]['filename']), args)
                else:
                    reduced_file_set += run_mash(list(metadata[species]['filename']), args)
            else:
                reduced_file_set += list(metadata[species]['filename'])
            
        metadata = pd.concat(metadata.values(), ignore_index=True)
        subset_df = metadata.set_index('filename').loc[reduced_file_set].reset_index()

        logging.debug(f"reduced set: {reduced_file_set}")
        subset_df.to_csv(os.path.join(args.out_dir, "filtered_gtdb_meta.tsv"), sep='\t', index=False)

    return os.path.join(args.out_dir, "filtered_gtdb_meta.tsv")

def generate_sylph_db(filtered_meta, outdir):
    df = pd.read_csv(filtered_meta, sep='\t')
    df['filename'].to_csv(os.path.join(outdir, "sylph_DB_file_locations.tsv"), index=False, header=False, sep='\t')
    df[['accession', 'gtdb_taxonomy']].to_csv(os.path.join(outdir, "sylph_DB_taxonomy.tsv"), index=False, header=False, sep='\t')

def find_gtdb_assemblies(in_dir):
    logging.info(f'\nLooking for assembly files in: {in_dir}')

    if os.path.isfile(in_dir):
        with open(in_dir, 'r') as infile:
            all_assemblies = [line.strip() for line in infile if os.path.isfile(line.strip())]
    else:
        all_assemblies = [str(x) for x in sorted(pathlib.Path(in_dir).glob('**/*')) if x.is_file()]

    # valid file formats
    all_assemblies = [x for x in all_assemblies if
                      x.endswith('.fasta') or x.endswith('.fasta.gz') or
                      x.endswith('.fna') or x.endswith('.fna.gz') or
                      x.endswith('.fa') or x.endswith('.fa.gz')]
    
    # valid GTDB assembly names
    all_assemblies = [x for x in all_assemblies if '_genomic' in x or '_ASM' in x]

    # Convert to dictionary
    all_assemblies = {'_'.join(os.path.basename(fasta).split('_', 2)[:2]): fasta for fasta in all_assemblies}

    plural = 'assembly' if len(all_assemblies) == 1 else 'assemblies'
    logging.info(f'Found {len(all_assemblies)} {plural}.')
    if len(all_assemblies) == 0:
        sys.exit(f'Error: no assemblies found in {in_dir}')
    return all_assemblies
