#!/usr/bin/env python3
"""
#############################################################################
Script to classify reads/sequences from fastq file, with a binning step to reduce memory consumption.
Bins the reads into defined bins and launch a classifier loading one bin at the time.

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

import argparse
import logging
import multiprocessing
import os
import os.path as osp
import pickle
import subprocess

# todo: add timing record
from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, init_logger


logger = init_logger('classify')


# For display / ipynb only
if False:
    plt.rcParams['figure.figsize'] = 13, 8
    pd.options.display.float_format = '{:,.2f}'.format


# #############################################################################
class MockCommunity:
    """ For a fastq file, bin reads, classify them, and compare results """
    
    def __init__(self, path_original_fastq, db_choice, folder_report, model, bin_nb=10, classifier="kraken2",
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
        raise NotImplementedError
    
    def __repr__(self):
        return f"Fastq file located at <{self.path_original_fastq}>, ready to be classified with " \
               f"{self.classifier} with the DB <{self.db_choice}>"
        

# #############################################################################
# Defaults and main method
path_fastq_comm = ["/home/ubuntu/Data/Segmentation/Test-Data/Synthetic_from_Genomes/2019-12-05_100000-WindowReads_20-BacGut/2019-12-05_100000-WindowReads_20-BacGut.fastq"]


def classify_reads(list_fastq, path_report, classifier, param_folder, db):
    """ Should load a file, do all the processing """
    logger.info("let's classify reads!")

    # Find the model
    path_model = ""
    for file in os.scandir(param_folder):
        if file.name.startswith("model_") and file.name.endswith(".pkl"):
            path_model = file.path
            break
    assert osp.isfile(path_model), FileNotFoundError(f"didn't find the ML model in {param_folder}... {path_model}")

    # Set the folder with hash tables
    kraken2_hash = osp.join(param_folder, "kraken2_hash")

    with open(path_model, 'b') as f:
        model = pickle.load(f)

    # Binning
    counter = 0
    read_to_bin = []
    CustomRead.set_fastq_path(path_fastq_comm)

    for record in tqdm(SeqIO.parse(path_fastq_comm, "fasta")):
        custom_read = CustomRead(record)
        custom_read.count_kmer()
        custom_read.lda_reduce()
        custom_read.find_bin()
        custom_read.to_fastq()

        counter += 1


    # Mock comm object
    for path_fastq in list_fastq:
        fastq_classifier = MockCommunity(
            path_original_fastq=path_fastq, db_choice=db, folder_report=path_report, model=model,
            bin_nb=10, classifier=classifier)
    
        fastq_classifier.dry_run=False
        fastq_classifier.db_choice = "2015-standard"
        fastq_classifier.classify()
        fastq_classifier.db_choice = "2015-10bins"
        fastq_classifier.classify()
    
    logger.error('Classification not implemented yet')


def test_classification():
    """ Should have a toy data set that i can bin, classify, and check the results """
    # todo: toy data set to check if it works
    raise NotImplementedError

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_folder',         help='Folder for output reports', type=is_valid_directory)
    parser.add_argument('-c', '--classifier',    help='choose which metagenomics classifier to use', metavar='',
                                                 choices=PATHS.classifiers, default=PATHS.classifiers[0])
    parser.add_argument('-m', '--model_folder',  help='Folder "clustered_by_<param>" with sub-folders "RefSeq/<bins> '
                                                      'and "model_<name>.pkl" ', metavar='',)
    parser.add_argument('-d', '--database',      help='which reference to use',
                                                 default='standard', metavar='', choices=('standard', 'mini', ))
    parser.add_argument('-t', '--threads',       help='Number of threads', default="10", metavar='', )
    parser.add_argument('-i', '--input_fastq',   help='List of input files in fastq format, space separated.',
                                                 default=path_fastq_comm, type=is_valid_file, nargs="+", metavar='',)

    args = parser.parse_args()
    logger.info(f'script called with following arguments: {args.input_fastq}, {args.output_folder}, {args.classifier}')

    classify_reads(args.input_fastq, args.output_folder, classifier=args.classifier, param_folder=args.model_folder,
                   db=args.database)


