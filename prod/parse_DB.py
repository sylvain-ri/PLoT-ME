#!/usr/bin/env python3
# #############################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-12-11
# Reads Binning Project
#
# #############################################################################
#
#
# Script to divide a Database into bins based on genome's k-mer frequencies
# 


import logging
from util import *


# Import paths and constants for the whole project
PATHS = ProjectPaths()

# #############################################################################
# https://docs.python.org/3/howto/logging-cookbook.html
# create formatter for the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# create file handler which logs even debug messages
fh = logging.FileHandler(PATHS.LOGS)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
ch.setFormatter(formatter)
# create logger with parse_DB.py and add the handlers to the logger
logger = logging.getLogger('parse_DB')
logger.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(ch)




def find_bins_DB(path_database, n_parts=10):
    """ Given a database of genomes in fastq files, split it in n segments """
    NotImplementedError



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        'Take a database of genomes and split it into bins according to the k-mer frequency of each genome. ' \
        'Needs a lot of disk space. Build the part for a specified classifier if specified')
    parser.add_argument('path_database',      help='Input file in fastq format', 
                        type=lambda x: x if osp.isfile(x) else FileNotFoundError(f'The path is not a file : {x}'))
    parser.add_argument('output_folder',      help='Folder for output reports', type=is_valid_directory)
    parser.add_argument('-c', '--classifier', help='choose kraken2 or centrifuge', choices=('kraken2'))
    parser.add_argument('-d', '--database',   default='standard', help='which reference to use', choices=['standard', 'mini', ])
    parser.add_argument('-t', '--threads',    default="10",       help='Number of threads')
    args = parser.parse_args()

    classify_reads(args.input_fastq, args.output_folder, classifier=args.classifier, )
        




