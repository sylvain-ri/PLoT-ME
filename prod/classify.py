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
from time import perf_counter

import numpy as np
from Bio import SeqRecord, SeqIO
from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, init_logger, scale_df_by_length, is_valid_directory, is_valid_file, create_path, \
    ArgumentParserWithDefaults, time_to_h_m_s
from bio import kmers_dic, seq_count_kmer


logger = init_logger('classify')


# Decorator for all these steps
# todo: find to which function to apply it to. How to track multiple files
def timing(func):
    """ Decorator to print steps and check if results have already been computed
        Need the second argument to be the output file/folder: it will check if the file exists / folder isn't empty
    """
    def wrapper(*args, **kwargs):
        # Check arguments for debugging
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ",\t".join(args_repr + kwargs_repr)

        # Time measurement
        start_time = perf_counter()
        try:
            logger.debug(f"START, function \t{func.__name__}({signature})")
            result = func(*args, **kwargs)
            # print time spent
            logger.debug(f"END, {time_to_h_m_s(start_time, perf_counter())}")
        except:
            result = None
        finally:
            timing.records.append({"func": func.__name__, "time": perf_counter()})  # log time for each step
        return result
    return wrapper


timing.records = [{"func": "start", "time": perf_counter()}]


# #############################################################################
class ReadToBin(SeqRecord.SeqRecord):
    """ General Read. Wrapping SeqIO.Record """
    K = 0
    KMER = {}  # kmers_dic(K)
    FASTQ_PATH = None
    FASTQ_BIN_FOLDER = None
    FILEBASE = ""
    MODEL = None
    PARAM = ""
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
        self.cluster = int(self.MODEL.predict(self.scaled)[0])
        self.description = f"bin_id={self.cluster}|{self.description}"
        self.path_out = f"{self.FASTQ_BIN_FOLDER}/{self.FILEBASE}.bin-{self.cluster}.fastq"
        # Save all output files
        ReadToBin.outputs[self.cluster] = self.path_out
        return self.cluster

    def to_fastq(self):
        assert self.path_out is not None, AttributeError("Path of the fastq file must first be defined")
        create_path(self.path_out)
        with open(self.path_out, "a") as f:
            SeqIO.write(self, f, "fasta")

    @classmethod
    def set_fastq_model_and_param(cls, path_fastq, path_model, param):
        assert osp.isfile(path_fastq), FileNotFoundError(f"{path_fastq} cannot be found")
        cls.PARAM = param
        cls.FASTQ_PATH = path_fastq
        folder, file_base = osp.split(osp.splitext(path_fastq)[0])
        cls.FASTQ_BIN_FOLDER = osp.join(folder, param)
        cls.FILEBASE = file_base
        logger.debug(f"New values: cls.FASTQ_PATH{cls.FASTQ_PATH} and cls.BASE_PATH{cls.FASTQ_BIN_FOLDER}")
        # /home/ubuntu/data/Segmentation/4mer_s10000/clustered_by_minikm_4mer_s10000/model_miniKM_4mer_s10000.pkl
        if path_model == "full":
            cls.K = 0
        else:
            k = path_model.split("/model_")[1].split("mer_")[0].split("_")[1]
            logger.debug(f"got path_model={path_model}, setting k={k}")
            cls.K = int(k)
            cls.KMER = kmers_dic(cls.K)
            with open(path_model, 'rb') as f:
                cls.MODEL = pickle.load(f)

    @classmethod
    def bin_reads(cls):
        """ Bin all reads from provide file """
        logger.info(f"Binning the reads (count kmers, scale, find_bin, copy to file.bin-<cluster>.fastq")
        # todo: try to parallelize it, careful of file writing concurrency
        # with Pool(cls.CORES) as pool:
        #     results = list(tqdm(pool.imap(pll_binning, SeqIO.parse(cls.FASTQ_PATH, "fasta"))))
        # counter = len(results)
        counter = 0
        total = 0
        for total, _ in tqdm(enumerate(SeqIO.parse(cls.FASTQ_PATH, "fasta")), desc="Counting number of reads",
                             leave=True):
            pass
        for record in tqdm(SeqIO.parse(cls.FASTQ_PATH, "fasta"), total=total, desc="binning and copying read to bins",
                           leave=True, dynamic_ncols=True):
            counter += 1
            custom_read = ReadToBin(record)
            # custom_read.kmer_count
            custom_read.scale()
            custom_read.find_bin()
            custom_read.to_fastq()
        logger.info(f"{counter} reads binned into bins: [" + ", ".join(map(str, sorted(cls.outputs.keys()))) + "]")
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
                 classifier_name="kraken2", param="", cores=cpu_count(), dry_run=False, verbose=False):
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
        self.folder_out      = osp.join(self.folder_report, self.file_name)
        if not os.path.isdir(self.folder_out):
            os.makedirs(self.folder_out)
        self.path_out        = osp.join(self.folder_out, f"{param}.{self.db_type}")
        
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
            for bin_id in self.path_binned_fastq.keys():
                self.classifier(self.path_binned_fastq[bin_id], osp.join(self.db_path, f"{bin_id}"), arg=f"bin-{bin_id}")
        elif "full" in self.db_type:
            self.classifier(self.path_original_fastq, self.db_path, arg="full")
        else:
            NotImplementedError("The database choice is either full or bins")
                
    def kraken2(self, file, path_hash, arg="unknown"):
        hash_file = osp.join(path_hash, "hash.k2d")
        self.logger.info(f'start to classify reads from file ({osp.getsize(file)/10**6:.2f} MB) {file}')
        self.logger.info(f'with kraken2, {arg}. hash table is ({osp.getsize(hash_file)/10**9:.2f} GB) {path_hash}')
        formatted_out = f"{self.path_out}.{arg}.kraken2" if self.db_type == "bins" else f"{self.path_out}.kraken2"
        self.logger.info(f'output is {formatted_out}.out')
        self.cmd = [
            "kraken2", "--threads", f"{self.cores}",
            "--db", path_hash,
            file,
            "--output", f"{formatted_out}.out",
            "--report", f"{formatted_out}.report",
        ]
        self.logger.debug(" ".join(self.cmd))
        if not self.dry_run:
            results = subprocess.check_output(self.cmd)
            self.logger.debug(results)
            
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


def bin_classify(list_fastq, path_report, path_database, classifier, db_type):
    """ Should load a file, do all the processing """
    logger.info("let's classify reads!")

    # Find the model
    if db_type == "bins":
        path_model = ""
        for file in os.scandir(path_database):
            if file.name.startswith("model_") and file.name.endswith(".pkl"):
                path_model = file.path
                break
        assert osp.isfile(path_model), FileNotFoundError(f"didn't find the ML model in {path_database}... {path_model}")
        path_to_hash = osp.join(path_database, f"{classifier}_hash")
    else:
        path_model = "full"
        path_to_hash = path_database

    # Set the folder with hash tables
    param = osp.basename(path_database)
    if param == "": param = osp.basename(path_database[:-1])
    logger.info(f"Assuming parameters are: {param}")

    for i, file in enumerate(list_fastq):
        logger.info(f"Opening fastq file ({i}/{len(list_fastq)}) {osp.basename(file)}")
        # Binning
        if "bins" in db_type:
            ReadToBin.set_fastq_model_and_param(file, path_model, param)
            fastq_binned = ReadToBin.bin_reads()
        else:
            fastq_binned = {}

        fastq_classifier = MockCommunity(
            path_original_fastq=file, db_path=path_to_hash, db_type=db_type, folder_report=path_report,
            path_binned_fastq=fastq_binned, bin_nb=10, classifier_name=classifier, param=param)

        fastq_classifier.classify()

    times = timing.records
    for i in range(len(times) - 1):
        logger.info(f"timing for function {times[i]['func']} - {time_to_h_m_s(times[i]['time'], times[i + 1]['time'])}")
    logger.info(f"Script ended, total time of {time_to_h_m_s(times[0]['time'], perf_counter())}.")
    print()


bin_classify.classifiers = ('kraken2',)


def test_classification():
    """ Should have a toy data set that i can bin, classify, and check the results """
    # todo: toy data set to check if it works
    raise NotImplementedError

    
if __name__ == '__main__':
    parser = ArgumentParserWithDefaults(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('output_folder',      help='Folder for output reports', type=is_valid_directory)
    parser.add_argument('database',           help='Folder with the hash table for the classifier, name '
                                                   '"clustered_by_<param>" with sub-folders "RefSeq/<bins> '
                                                   'and "model_<name>.pkl" ')
    parser.add_argument('-c', '--classifier', help='choose which metagenomics classifier to use', metavar='',
                                              choices=bin_classify.classifiers, default=bin_classify.classifiers[0])
    parser.add_argument('-t', '--db_type',    help='Choose to use the standard full database or the segmented one',
                                              default='bins', choices=('full', 'bins',) , metavar='')
    parser.add_argument('-i', '--input_fastq',help='List of input files in fastq format, space separated.',
                                              default=path_fastq_comm, type=is_valid_file, nargs="+", metavar='')
    # parser.add_argument('-c', '--cores',         help='Number of cores', default=cpu_count(), metavar='')

    args = parser.parse_args()
    logger.info(f'script called with following arguments: {args.input_fastq}, {args.output_folder}, {args.classifier}')

    bin_classify(args.input_fastq, args.output_folder, args.database,
                 classifier=args.classifier, db_type=args.db_type)




