import argparse
import logging
import os
from .core import find_gtdb_assemblies, load_metadata, dereplication, generate_sylph_db

def get_arguments(args=None):
    parser = argparse.ArgumentParser(description='Rapid GTDB Assembly Dereplicator')

    required_args = parser.add_argument_group('Input/output')

    required_args.add_argument('-i', '--in_dir', type=str, required=True,
                               help='Directory containing all assemblies. '
                               'Can also provide a text file with paths to assemblies on each line.')
    
    required_args.add_argument('-o', '--out_dir', type=str, required=True,
                               help='Directory where dereplicated assemblies will be copied')
    
    required_args.add_argument('--meta', type=str, required=True,
                                 help='GTDB metadata file')

    clustering_args = parser.add_argument_group('Clustering parameters')
    
    clustering_args.add_argument('--distance',
                                 default=99,
                                 type=lambda x: float(x) if 0 <= float(x) <= 100 else argparse.ArgumentTypeError('--distance must be between 0 and 100'),
                                 help='Dereplicate until the closest pair has a pairwise identity of at lease this value (default=99)')
    
    clustering_args.add_argument('--method', choices=['mash', 'skani'], default='mash',
                                 help='Method for pairwise distance estimation (default: mash). '
                                      'Mash is faster but less accurate than skani.')

    clustering_args.add_argument('--gtdb_reps', action='store_true',
                              help='Simply outputs the current GTDB species representative genomes as specified in the meta data file.')

    setting_args = parser.add_argument_group('Settings')


    setting_args.add_argument('-t', '--threads', type=int, default=1,
                              help='Number of CPU threads for Mash (default: 1)')
    
    setting_args.add_argument('--sketch_size', type=int, default=10000,
                              help='Mash assembly sketch size (default: 1000)')
    
    setting_args.add_argument('--sylphDB', action='store_true',
                              help='Turns on output of files required to create a Sylph database (default: False)')

    setting_args.add_argument(
        "--loglevel",
        type=str.upper,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging threshold.",
    )

    return parser.parse_args(args)


def main(args=None):
    parsed_args = get_arguments(args)

    # set logging up
    logging.basicConfig(
        level=parsed_args.loglevel,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    all_assemblies = find_gtdb_assemblies(parsed_args.in_dir)
    os.makedirs(parsed_args.out_dir, exist_ok=True)
    parsed_args.out_dir = os.path.join(parsed_args.out_dir, "")
    
    metadata = load_metadata(parsed_args.meta, all_assemblies)
    derep_meta_file = dereplication(all_assemblies, metadata, parsed_args)

    logging.info("Successfully completed dereplication!")

    if parsed_args.sylphDB:
        generate_sylph_db(derep_meta_file, parsed_args.out_dir)
        logging.info("Generated Sylph database files.")

if __name__ == '__main__':
    main()
