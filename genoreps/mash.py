import os
import tempfile
import logging
import subprocess
import shutil
import pandas as pd

def run_mash(assemblies, args):
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.out_dir), "")

    with open(temp_dir + "assembly_path_lists.txt", 'w') as outfile:
        for a in assemblies:
            outfile.write(a + "\n")
    
    # Build database
    cmd = "mash sketch"
    cmd += " -l " + temp_dir + "assembly_path_lists.txt"
    cmd += " -o " + temp_dir + "mashdb.msh"
    cmd += " -p " + str(args.threads)

    logging.info("Building mash database...")
    logging.debug("Running command: %s", cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Greedily filter genomes at given distance threshold
    keep = set()
    while len(assemblies) > 0:
        current = assemblies.pop(0)
        keep.add(current)

        cmd = "mash dist "
        cmd += " -p " + str(args.threads)
        if len(keep)==1:
            cmd += " -d 1"
        else:
            cmd += " -d " + str(1-args.distance/100.0)
        cmd += " " + temp_dir + "mashdb.msh"
        cmd += " " + current 
        cmd += " > " + temp_dir + "distances.tab"

        logging.info("Remaining sequences to consider: %d", len(assemblies))
        
        logging.debug("Running command: %s", cmd)
        subprocess.run(cmd, shell=True, check=True)

        # Load output and remove matching assemblies from list
        df = pd.read_csv(temp_dir + "distances.tab",  sep='\t', header=None)
        df.columns = ['reference', 'query', 'dist', 'p-value', 'hash']

        if len(keep)==1:
            assembly_order = list(dict.fromkeys([current] + list(df.loc[df['dist'].argsort(), 'reference'])))

        matching_assemblies = set(df.loc[df['dist'] < 1-args.distance/100.0, 'reference'])

        assemblies = [fasta for fasta in assemblies if fasta not in matching_assemblies]

    shutil.rmtree(temp_dir)

    # Reorder "keep" based on the order in "assembly_order"
    ordered_keep = [item for item in assembly_order if item in keep]

    return ordered_keep
