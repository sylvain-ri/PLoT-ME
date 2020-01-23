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
from multiprocessing.pool import Pool
from os import path as osp
import pickle
import subprocess

# todo: add timing record

import numpy as np
from Bio import SeqRecord, SeqIO
from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, init_logger, scale_df_by_length, is_valid_directory, is_valid_file
from bio import kmers_dic, seq_count_kmer


logger = init_logger('classify')


# #############################################################################
class ReadToBin(SeqRecord.SeqRecord):
    """ General Read. Wrapping SeqIO.Record """
    K = 4
    KMER = kmers_dic(K)
    FASTQ_PATH = None
    BASE_PATH = None
    MODEL = None
    CORES = cpu_count()
    outputs = {}

    def __init__(self, obj):
        # wrap the object
        self._wrapped_obj = obj
        # Additional attributes
        self.cluster = None
        self._kmer_count = None
        self.scaled = None
        self.path_out = None

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self._wrapped_obj, attr)

    @property
    def kmer_count(self, ignore_N=True):
        """ common method """
        if self._kmer_count is None:
            self._kmer_count = seq_count_kmer(self.seq, self.KMER.copy(), self.K, ignore_N=ignore_N)
        return self._kmer_count

    def scale(self):
        logger.debug("scaling the read by it's length and k-mer")
        self.scaled = scale_df_by_length(np.fromiter(self.kmer_count.values(), dtype=int).reshape(-1, 4**self.K),
                                         None, k=self.K, w=len(self.seq), single_row=True)  # Put into 2D one row
        return self.scaled

    def find_bin(self):
        logger.debug('finding bins for each read')
        self.cluster = self.MODEL.predict(self.scaled)[0]
        self.description = f"bin_id={self.cluster}|{self.description}"
        self.path_out = f"{self.BASE_PATH}.bin-{self.cluster}.fastq"
        # Save all output files
        ReadToBin.outputs[self.cluster] = self.path_out
        return self.cluster

    def to_fastq(self):
        assert self.path_out is not None, AttributeError("Path of the fastq file must first be defined")
        with open(self.path_out, "a") as f:
            SeqIO.write(self, f, "fasta")

    @classmethod
    def set_model(cls, path_model):
        # /home/ubuntu/data/Segmentation/4mer_s10000/clustered_by_minikm_4mer_s10000/model_miniKM_4mer_s10000.pkl
        k = path_model.split("/model_")[1].split("mer_")[0].split("_")[1]
        logger.debug(f"got path_model={path_model}, setting k={k}")
        cls.K = int(k)
        cls.KMER = kmers_dic(cls.K)
        with open(path_model, 'rb') as f:
            cls.MODEL = pickle.load(f)

    @classmethod
    def set_fastq_path(cls, path_fastq):
        assert osp.isfile(path_fastq), FileNotFoundError(f"{path_fastq} cannot be found")
        cls.FASTQ_PATH = path_fastq
        cls.BASE_PATH = osp.splitext(path_fastq)[0]
        logger.debug(f"New values: cls.FASTQ_PATH{cls.FASTQ_PATH} and cls.BASE_PATH{cls.BASE_PATH}")

    @classmethod
    def bin_reads(cls):
        """ Bin all reads from provide file """
        logger.info(f"Binning all the read (count kmers, scale, find_bin, copy to file.bin-<cluster>.fastq")
        with Pool(cls.CORES) as pool:
            results = list(tqdm(pool.imap(pll_binning, SeqIO.parse(cls.FASTQ_PATH, "fasta"))))
        # for record in tqdm(SeqIO.parse(cls.FASTQ_PATH, "fasta")):
        #     counter += 1
        logger.info(cls.outputs)
        logger.info(f"{len(results)} reads binned into bins: " + ", ".join(cls.outputs.keys()))
        return cls.outputs


def pll_binning(record):
    """ Parallel processing of read binning """
    custom_read = ReadToBin(record)
    # custom_read.kmer_count
    custom_read.scale()
    custom_read.find_bin()
    custom_read.to_fastq()


# #############################################################################
class MockCommunity:
    """ For a fastq file, bin reads, classify them, and compare results """
    
    def __init__(self, path_original_fastq, db_path, db_type, folder_report, path_binned_fastq={}, bin_nb=10,
                 classifier_name="kraken2", cores=cpu_count(), dry_run=False, verbose=False):
        self.logger = logging.getLogger('classify.MockCommunity')

        assert osp.isfile(path_original_fastq), FileNotFoundError(f"Didn't find original fastq {path_original_fastq}")
        self.path_original_fastq    = path_original_fastq

        self.folder, self.file_name = osp.split(osp.splitext(self.path_original_fastq)[0])
        self.path_binned_fastq      = path_binned_fastq              # {<bin i>: <path_file>}
        self.folder_report          = folder_report
        
        self.classifier_name = classifier_name
        self.db_path         = db_path    # location of the hash table for the classifier
        self.db_type         = db_type    # Either full or bins
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
        self.logger.info(f"Classifying reads with {self.db_type} setting")
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
path_fastq_comm = ["/home/ubuntu/data/Segmentation/Test-Data/Synthetic_from_Genomes/"
                   "2019-12-19_20-WindowReads_EColi_Test/2019-12-19_20-WindowReads_10-EColiTest.fastq"]


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
                                                      'and "model_<name>.pkl" ')
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




