import os
import tempfile
import logging
import subprocess
import shutil
import pandas as pd

def run_scani(assemblies, args):
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.out_dir), "")

    with open(temp_dir + "signature_list.txt", 'w') as siglist:
        with open(temp_dir + "assembly_path_lists.txt", 'w') as outfile:
            for a in assemblies:
                outfile.write(a + "\n")
                siglist.write(temp_dir + "skanidb/" + os.path.basename(a) + ".sketch\n") 
    
    # Build database
    cmd = "skani sketch"
    cmd += " -l " + temp_dir + "assembly_path_lists.txt"
    cmd += " -o " + temp_dir + "skanidb"
    cmd += " -t " + str(args.threads)

    logging.info("Building skani database...")
    logging.debug("Running command: %s", cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Greedily filter genomes at given distance threshold
    keep = []
    while len(assemblies) > 0:
        current = assemblies.pop(0)
        keep.append(current)

        cmd = "skani dist --fast"
        cmd += " -t " + str(args.threads)
        cmd += " --rl " + temp_dir + "signature_list.txt"
        cmd += " -q " + temp_dir + "skanidb/" + os.path.basename(current) + ".sketch"
        cmd += " -s " + str(0.95) 
        cmd += " -o " + temp_dir + "skani_search_out.tsv"

        logging.info("Remaining sequences to consider: %d", len(assemblies))
        
        logging.debug("Running command: %s", cmd)
        subprocess.run(cmd, shell=True, check=True)

        # Load output and remove matching assemblies from list
        df = pd.read_csv(temp_dir + "skani_search_out.tsv",  sep='\t')
        matching_assemblies = set(df.loc[df['ANI'] > args.distance, 'Ref_file'])

        assemblies = [fasta for fasta in assemblies if fasta not in matching_assemblies]

    shutil.rmtree(temp_dir)

    return keep
