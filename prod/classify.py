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
from multiprocessing import cpu_count
import os
import os.path as osp
import pickle
import subprocess

# todo: add timing record
from Bio import SeqIO
from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, init_logger
from bio import ReadToBin


logger = init_logger('classify')


# For display / ipynb only
if False:
    plt.rcParams['figure.figsize'] = 13, 8
    pd.options.display.float_format = '{:,.2f}'.format


# #############################################################################
class MockCommunity:
    """ For a fastq file, bin reads, classify them, and compare results """
    
    def __init__(self, path_original_fastq, db_path, db_type, folder_report, path_binned_fastq={}, bin_nb=10,
                 classifier_name="kraken2", cores=cpu_count(), dry_run=False, verbose=False):
        self.logger = logging.getLogger('classify.FastQClassification')

        assert osp.isfile(path_original_fastq), FileNotFoundError(f"Didn't find original fastq {path_original_fastq}")
        self.path_original_fastq    = path_original_fastq

        self.folder, self.file_name = osp.split(osp.splitext(self.path_original_fastq)[0])
        self.path_binned_fastq      = path_binned_fastq              # {<bin i>: <path_file>}
        self.folder_report          = folder_report
        
        self.classifier_name = classifier_name
        self.db_path         = db_path    # location of the hash table for the classifier
        self.db_type         = db_type  # Either full or bins: 2015-standard / 2015-bins10
        self.bin_nb          = bin_nb
        self.folder_out      = f"{self.folder_report}/{self.file_name}"
        if not os.path.isdir(self.folder_out):
            os.makedirs(self.folder_out)
        self.path_out        = f"{self.folder_out}/{self.db_type}"
        
        self.cores           = cores
        self.dry_run         = dry_run
        self.verbose         = verbose
        self.cmd             = None

    @property
    def classifier(self):
        if self.classifier_name == "kraken2":
            return self.kraken2
        else:
            NotImplementedError("This classifier hasn't been implemented")

    def classify(self):
        if "bins" in self.db_type:
            for bin_id in tqdm(self.path_binned_fastq.keys()):
                self.classifier(self.path_binned_fastq[bin_id], osp.join(self.db_path, f"{bin_id}"), arg=f"bin-{bin_id}")
        elif "full" in self.db_type:
            self.classifier(self.path_original_fastq, osp.join(self.db_path, "full"), arg="full")
        else:
            NotImplementedError("The database choice is either full or bins")
                
    def kraken2(self, file, path_hash, arg="unknown"):
        self.logger.info('start to classify reads with kraken2')
        self.cmd = [
            "kraken2", "--threads", f"{self.cores}",
            "--db", path_hash,
            file,
            "--output", f"{self.path_out}.{arg}.kraken2.out",
            "--report", f"{self.path_out}.{arg}.kraken2.report",
        ]
        self.logger.info(" ".join(self.cmd))
        if not self.dry_run:
            results = subprocess.check_output(self.cmd)
            if self.verbose: print(results)
            
    def kraken2_report_merging(self):
        self.logger.info('Merging kraken2 reports')
        raise NotImplementedError()
    
    def __repr__(self):
        return f"Fastq file located at <{self.path_original_fastq}>, ready to be classified with " \
               f"{self.classifier_name} with the DB <{self.db_type}> located at {self.db_path}"
        

# #############################################################################
# Defaults and main method
path_fastq_comm = ["/home/ubuntu/Data/Segmentation/Test-Data/Synthetic_from_Genomes/2019-12-05_100000-WindowReads_20-BacGut/2019-12-05_100000-WindowReads_20-BacGut.fastq"]


def classify_reads(list_fastq, path_report, path_database, classifier, db_type):
    """ Should load a file, do all the processing """
    logger.info("let's classify reads!")

    # Find the model
    path_model = ""
    for file in os.scandir(path_database):
        if file.name.startswith("model_") and file.name.endswith(".pkl"):
            path_model = file.path
            break
    assert osp.isfile(path_model), FileNotFoundError(f"didn't find the ML model in {path_database}... {path_model}")

    # Set the folder with hash tables
    path_to_hash = osp.join(path_database, f"{classifier}_hash")

    logger.info("let's classify reads!")
    for file in tqdm(list_fastq):
        # Binning
        if "bins" in db_type:
            ReadToBin.set_fastq_path(file)
            ReadToBin.set_model(path_model)
            fastq_binned = ReadToBin.bin_reads()
        else:
            fastq_binned = {}

        fastq_classifier = MockCommunity(
            path_original_fastq=file, db_path=path_to_hash, db_type=db_type,
            folder_report=path_report, path_binned_fastq=fastq_binned, bin_nb=10, classifier_name=classifier)

        fastq_classifier.classify()


def test_classification():
    """ Should have a toy data set that i can bin, classify, and check the results """
    # todo: toy data set to check if it works
    raise NotImplementedError

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_folder',         help='Folder for output reports', type=is_valid_directory)
    parser.add_argument('database',              help='Folder with the hash table for the classifier, name '
                                                      '"clustered_by_<param>" with sub-folders "RefSeq/<bins> '
                                                      'and "model_<name>.pkl" ', metavar='')
    parser.add_argument('-c', '--classifier',    help='choose which metagenomics classifier to use', metavar='',
                                                 choices=PATHS.classifiers, default=PATHS.classifiers[0])
    parser.add_argument('-t', '--db_type',          help='Choose to use the standard full database or the segmented one',
                                                 default='bins', choices=('full', 'bins',) , metavar='')
    parser.add_argument('-i', '--input_fastq',   help='List of input files in fastq format, space separated.',
                                                 default=path_fastq_comm, type=is_valid_file, nargs="+", metavar='')
    # parser.add_argument('-c', '--cores',         help='Number of cores', default=cpu_count(), metavar='')

    args = parser.parse_args()
    logger.info(f'script called with following arguments: {args.input_fastq}, {args.output_folder}, {args.classifier}')

    classify_reads(args.input_fastq, args.output_folder, args.database,
                   classifier=args.classifier, db_type=args.db_type)


