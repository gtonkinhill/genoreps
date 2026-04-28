# Rapid GTDB Assembly Dereplicator

Command-line tool for greedy dereplication of GTDB assemblies. The tool follows an approach similar to the dRep algorithm but is specifically designed for large, standardised collections such as GTDB. It ranks genomes using the GTDB assembly scoring system, ensuring that type strains are retained while prioritising higher-quality assemblies.

# genoreps

```
usage: genoreps [-h] -i IN_DIR -o OUT_DIR --meta META [--distance DISTANCE]
                [--method {mash,skani}] [--gtdb_reps] [-t THREADS]
                [--sketch_size SKETCH_SIZE] [--sylphDB]
                [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Rapid GTDB Assembly Dereplicator

options:
  -h, --help            show this help message and exit

Input/output:
  -i, --in_dir IN_DIR   Directory containing all assemblies. Can also provide
                        a text file with paths to assemblies on each line.
  -o, --out_dir OUT_DIR
                        Directory where dereplicated assemblies will be copied
  --meta META           GTDB metadata file

Clustering parameters:
  --distance DISTANCE   Dereplicate until the closest pair has a pairwise
                        identity of at lease this value (default=99)
  --method {mash,skani}
                        Method for pairwise distance estimation (default:
                        mash). Mash is faster but less accurate than skani.
  --gtdb_reps           Simply outputs the current GTDB species representative
                        genomes as specified in the meta data file.

Settings:
  -t, --threads THREADS
                        Number of CPU threads for Mash (default: 1)
  --sketch_size SKETCH_SIZE
                        Mash assembly sketch size (default: 1000)
  --sylphDB             Turns on output of files required to create a Sylph
                        database (default: False)
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging threshold.
```