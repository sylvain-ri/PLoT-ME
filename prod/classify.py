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
# Script to classify sequences from fastq file
# 

import argparse
import logging
import multiprocessing
import os
import os.path as osp
import subprocess

# todo: add timing record
from tqdm import tqdm_notebook as tqdm

# Import paths and constants for the whole project
from prod.util import PATHS


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
# create logger with classify.py and add the handlers to the logger
logger = logging.getLogger('classify')
logger.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(ch)


# For display / ipynb only
if False:
    plt.rcParams['figure.figsize'] = 13, 8
    pd.options.display.float_format = '{:,.2f}'.format


# #############################################################################
class FastQClassification:
    """ For a fastq file, bin reads, classify them, and compare results """
    
    def __init__(self, path_original_fastq, db_choice, folder_report, bin_nb=10, classifier="kraken2",
                 cores=multiprocessing.cpu_count(), dry_run=True, verbose=False):
        self.logger = logging.getLogger('classify.FastQClassification')
        assert osp.isfile(path_original_fastq), FileNotFoundError(f"Didn't find original fastq {path_original_fastq}")
        
        self.path_original_fastq    = path_original_fastq
        self.folder, self.file_name = osp.split(osp.splitext(self.path_original_fastq)[0])
        self.path_binned_fastq      = []              # [(<bin i>, <path_file>), ]
        self.folder_report          = folder_report
        
        self.classifier    = classifier
        self.db_choice     = db_choice  # Either full or bins: 2015-standard / 2015-bins10
        self.bin_nb        = bin_nb
        self.folder_out    = f"{self.folder_report}/{self.file_name}"
        if not os.path.isdir(self.folder_out):
            os.makedirs(self.folder_out)
        self.path_out      = f"{self.folder_out}/{self.db_choice}"
        
        self.cores         = cores
        self.dry_run       = dry_run
        self.verbose       = verbose
        self.cmd           = None

    def find_binned_files(self):
        for i in range(self.bin_nb):
            path_bin_i = f"{self.folder}/{self.file_name}.bin-{i}.fastq"
            if osp.isfile(path_bin_i): 
                self.path_binned_fastq.append((i, path_bin_i))
    
    def classify(self):
        if self.classifier == "kraken2":
            if "bins" in self.db_choice:
                self.find_binned_files()
                for i, file in tqdm(self.path_binned_fastq):
                    print(f"bin {i} db loading... {file}")
                    self.kraken2_bins(i, file)
            else:
                print("standard full db loading...")
                self.kraken2_standard()
                
    def kraken2_standard(self):
        self.logger.info('start to classify reads with kraken2 standard DB')
        self.cmd = [
            "kraken2", "--threads", f"{self.cores}",
            "--db", f"{PATHS.kraken2_DB[self.db_choice]}",
            self.path_original_fastq, 
            "--output", f"{self.path_out}.full.kraken2.out",
            "--report", f"{self.path_out}.full.kraken2.report",
        ]
        if self.verbose: print(" ".join(self.cmd))
        if not self.dry_run:            
            results = subprocess.check_output(self.cmd)
            if self.verbose: print(results)
                
    def kraken2_bins(self, i, file):
        self.logger.info('start to classify reads with kraken2 10 bins DB')
        self.cmd = [
            "kraken2", "--threads", f"{self.cores}",
            "--db", f"{PATHS.kraken2_DB[self.db_choice]}/{i}",
            file, 
            "--output", f"{self.path_out}.bin-{i}.kraken2.out",
            "--report", f"{self.path_out}.bin-{i}.kraken2.report",
        ]
        if self.verbose: print(" ".join(self.cmd))
        if not self.dry_run:            
            results = subprocess.check_output(self.cmd)
            if self.verbose: print(results)
            
    def kraken2_report_merging(self):
        self.logger.info('Merging kraken2 reports')
        NotImplementedError
    
    def __repr__(self):
        return f"Fastq file located at <{self.path_original_fastq}>, ready to be classified with " \
               f"{self.classifier} with the DB <{self.db_choice}>"
        

# #############################################################################
# Defaults and main method
path_fastq_comm = "/home/ubuntu/Data/Segmentation/Test-Data/Synthetic_from_Genomes/2019-12-05_100000-WindowReads_20-BacGut/2019-12-05_100000-WindowReads_20-BacGut.fastq"


def classify_reads(path_fastq, path_report, classifier, db):
    """ Should load a file, do all the processing """
    logger.info("let's classify reads!")
    fastq_classifier = FastQClassification(
        path_original_fastq=path_fastq, db_choice=db, folder_report=path_report, 
        bin_nb=10, classifier=classifier)
    
    fastq_classifier.dry_run=False
    fastq_classifier.db_choice = "2015-standard"
    fastq_classifier.classify()
    fastq_classifier.db_choice = "2015-10bins"
    fastq_classifier.classify()
    
    logger.error('Classification not implemented yet')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     'Provided a fastq file with reads, bins the reads into defined bins '
                                     'and launch a classifier loading one bin at the time.')
    parser.add_argument('-i', '--input_fastq',   help='Input file in fastq format', type=is_valid_file, 
                                                 default=path_fastq_comm)
    parser.add_argument('-o', '--output_folder', help='Folder for output reports', type=is_valid_directory,
                        default=folder_today(PATHS.folder_reports))
    parser.add_argument('-c', '--classifier',    help='choose kraken2 or centrifuge', 
                                                 choices=PATHS.classifiers, default=PATHS.classifiers[0])
    parser.add_argument('-d', '--database',      default='standard', help='which reference to use', choices=('standard', 'mini', ))
    parser.add_argument('-t', '--threads',       default="10",       help='Number of threads')
    
    args = parser.parse_args()
    logger.info(f'script called with following arguments: {args.input_fastq}, {args.output_folder}, {args.classifier}')

    classify_reads(args.input_fastq, args.output_folder, classifier=args.classifier, db=args.database)


