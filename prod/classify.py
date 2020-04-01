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
import csv
from datetime import datetime as dt
import logging
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
import os
from os import path as osp
import pickle
import subprocess
from time import perf_counter
import re

import numpy as np
from Bio import SeqRecord, SeqIO
from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, init_logger, scale_df_by_length, is_valid_directory, is_valid_file, create_path, \
    ArgumentParserWithDefaults, time_to_hms
from bio import kmers_dic, seq_count_kmer


logger = init_logger('classify')


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
    CORES = 1
    outputs = {}
    NUMBER_BINNED = 0

    def __init__(self, obj):
        # wrap the object
        self._wrapped_obj = obj
        # Additional attributes
        self.cluster = None
        self._kmer_count = None
        self.scaled = None

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

    @property
    def path_out(self, cluster=None):
        return f"{self.FASTQ_BIN_FOLDER}/{self.FILEBASE}.bin-{self.cluster if cluster is None else cluster}.fastq"

    def scale(self):
        logger.log(5, "scaling the read by it's length and k-mer")
        self.scaled = scale_df_by_length(np.fromiter(self.kmer_count.values(), dtype=int).reshape(-1, 4**self.K),
                                         None, k=self.K, w=len(self.seq), single_row=True)  # Put into 2D one row
        return self.scaled

    def find_bin(self):
        logger.log(5, 'finding bins for each read')
        self.cluster = int(self.MODEL.predict(self.scaled)[0])
        self.description = f"bin_id={self.cluster}|{self.description}"
        # self.path_out = f"{self.FASTQ_BIN_FOLDER}/{self.FILEBASE}.bin-{self.cluster}.fastq"
        # Save all output files
        ReadToBin.outputs[self.cluster] = self.path_out
        return self.cluster

    def to_fastq(self):
        assert self.path_out is not None, AttributeError("Path of the fastq file must first be defined")
        with open(self.path_out, "a") as f:
            SeqIO.write(self, f, bin_classify.format)

    @classmethod
    def set_fastq_model_and_param(cls, path_fastq, path_model, param, cores, k):
        assert osp.isfile(path_fastq), FileNotFoundError(f"{path_fastq} cannot be found")
        # todo: load the parameter file from parse_DB.py instead of parsing string.... parameters_RefSeq_binning.txt
        cls.CORES = cores
        cls.PARAM = param
        cls.FASTQ_PATH = path_fastq
        folder, file_base = osp.split(osp.splitext(path_fastq)[0])
        # output folder, will host one file for each bin
        cls.FASTQ_BIN_FOLDER = osp.join(folder, param)
        if osp.isdir(cls.FASTQ_BIN_FOLDER):
            last_modif = dt.fromtimestamp(osp.getmtime(cls.FASTQ_BIN_FOLDER))
            save_folder = f"{cls.FASTQ_BIN_FOLDER}_{last_modif:%Y-%m-%d_%H-%M}"
            logger.warning(f"Folder existing, renaming to avoid losing files: {save_folder}")
            os.rename(cls.FASTQ_BIN_FOLDER, save_folder)
        create_path(cls.FASTQ_BIN_FOLDER)

        cls.FILEBASE = file_base
        logger.debug(f"New values: cls.FASTQ_PATH{cls.FASTQ_PATH} and cls.BASE_PATH{cls.FASTQ_BIN_FOLDER}")
        # /home/ubuntu/data/Segmentation/4mer_s10000/clustered_by_minikm_4mer_s10000/model_miniKM_4mer_s10000.pkl
        if path_model == "full":
            cls.K = 0
        else:
            cls.K = int(k)
            cls.KMER = kmers_dic(cls.K)
            with open(path_model, 'rb') as f:
                cls.MODEL = pickle.load(f)

    @classmethod
    def bin_reads(cls):
        """ Bin all reads from provide file """
        logger.info(f"Binning the reads (count kmers, scale, find_bin, copy to file.bin-<cluster>.fastq")
        # todo: try to parallelize it, careful of file writing concurrency.
        #  Dask ? process to load and count kmers, single one for appending read to fastq ?
        # with Pool(cls.CORES) as pool:
        #     results = list(tqdm(pool.imap(pll_binning, SeqIO.parse(cls.FASTQ_PATH, "fasta"))))
        # counter = len(results)
        counter = 0
        # Count number of lines with bash wc -l and divide by 4 if fastq, otherwise by 2 (fasta)
        total = round(int(subprocess.check_output(["wc", "-l", cls.FASTQ_PATH]).split()[0]) /
                      (4 if bin_classify.format == "fastq" else 2))
        for record in tqdm(SeqIO.parse(cls.FASTQ_PATH, bin_classify.format), total=total,
                           desc="binning and copying reads to bins", leave=True, dynamic_ncols=True):
            counter += 1
            custom_read = ReadToBin(record)
            # custom_read.kmer_count
            custom_read.scale()
            custom_read.find_bin()
            custom_read.to_fastq()
        logger.info(f"{counter} reads binned into bins: [" + ", ".join(map(str, sorted(cls.outputs.keys()))) + "]")
        cls.NUMBER_BINNED = counter
        return cls.outputs

    @classmethod
    def sort_bins_by_sizes(cls):
        """ Sort the bins by their size """
        # todo: drop bins with less than 0.1% of the total size ?
        bin_size = {}
        for f in os.scandir(cls.FASTQ_BIN_FOLDER):
            bin_nb = int(f.name.split('.')[1].split('-')[1])
            size = osp.getsize(f)
            bin_size[size] = bin_nb
        # make a copy first, then empty the dic, and rewrite it in the correct order
        fastq_outputs = ReadToBin.outputs
        ReadToBin.outputs = {}
        for size, bin_nb in sorted(bin_size.items(), reverse=True):
            ReadToBin.outputs[bin_nb] = fastq_outputs[bin_nb]
        return ReadToBin.outputs


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
                 classifier_name="kraken2", param="", cores=1, clf_settings="default", dry_run=False, verbose=False):
        self.logger = logging.getLogger('classify.MockCommunity')

        assert osp.isfile(path_original_fastq), FileNotFoundError(f"Didn't find original fastq {path_original_fastq}")
        logger.debug(f"Path to hash files: db_path = {db_path}")

        self.path_original_fastq    = path_original_fastq

        self.folder, self.file_name = osp.split(osp.splitext(self.path_original_fastq)[0])
        self.path_binned_fastq      = path_binned_fastq              # {<bin i>: <path_file>}
        self.folder_report          = folder_report
        
        self.classifier_name = classifier_name
        self.db_path         = db_path    # location of the hash table for the classifier
        self.db_type         = db_type    # Either full or bins
        self.hash_files      = {}
        self.bin_nb          = bin_nb
        self.folder_out      = osp.join(self.folder_report, self.file_name)
        if not os.path.isdir(self.folder_out):
            os.makedirs(self.folder_out)
        self.path_out        = osp.join(self.folder_out, f"{param}.{classifier_name}.{clf_settings}.{self.db_type}")
        
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
                path_hash_bin = osp.join(self.db_path, f"{bin_id}")
                logger.debug(f"Path of fastq bin : {self.path_binned_fastq[bin_id]}")
                logger.debug(f"Path of hash bin : {path_hash_bin}")
                assert osp.isfile(path_hash_bin), FileNotFoundError(f"Hash table not found ! {path_hash_bin}")
                self.classifier(self.path_binned_fastq[bin_id], path_hash_bin, arg=f"bin-{bin_id}")
        elif "full" in self.db_type:
            self.classifier(self.path_original_fastq, self.db_path, arg="full")
        else:
            NotImplementedError("The database choice is either full or bins")
                
    def kraken2(self, file, path_hash, arg="unknown"):
        hash_file = osp.join(path_hash, "hash.k2d")
        self.hash_files[arg] = hash_file
        self.logger.info(f'start to classify reads from file ({osp.getsize(file)/10**6:.2f} MB) {file}')
        self.logger.info(f'with kraken2, {arg}. hash table is ({osp.getsize(hash_file)/10**9:.2f} GB) {path_hash}')
        formatted_out = f"{self.path_out}.{arg}" if self.db_type == "bins" else f"{self.path_out}"
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


def bin_classify(list_fastq, path_report, path_database, classifier, db_type,
                 cores=cpu_count(), f_record="/home/ubuntu/classify_records.csv", clf_settings=""):
    """ Should load a file, do all the processing """
    print("\n*********************************************************************************************************")
    logger.info("**** Starting script **** \n ")
    logger.info(f"Script {__file__} called with {args}")
    bin_classify.cores = cores

    # preparing csv record file
    if not osp.isfile(f_record):
        with open(f_record, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            headers = ("FILE", "BINS_vs_FULL", "BINNING", "CLASSIFY", "TOTAL", "HASHES_SIZE", "NB_BINS", "HASH_PATH", "HASH_NAME")
            csv_writer.writerow(headers)

    logger.info("let's classify reads!")

    # Find the model
    if db_type == "bins":
        path_model = ""
        for file in os.scandir(path_database):
            if file.name.startswith("model.") and file.name.endswith(".pkl"):
                path_model = file.path
                break
        assert osp.isfile(path_model), FileNotFoundError(f"didn't find the ML model in {path_database}... {path_model}")

        # Parse the model name to find parameters:
        basename = path_model.split("/model.")[1]
        clusterer, bin_nb, k, w, omitted = re.split('_b|_k|_s|_o', basename)

        path_to_hash = osp.join(path_database, classifier, clf_settings)
        logger.info(f"WTH: {path_to_hash}\t{path_database}\t{classifier}\t{clf_settings}")
    else:
        path_model = "full"
        clusterer, bin_nb, k, w, omitted = (None, 1, None, None, "oplant-vertebrate")
        path_to_hash = osp.join(path_database, omitted, classifier, clf_settings)

    # Set the folder with hash tables
    param = osp.basename(path_database)
    if param == "": param = osp.basename(path_database[:-1])
    logger.info(f"Assuming parameters are: {param}")

    t = {}  # recording time at each step
    for i, file in enumerate(list_fastq):
        try:
            assert osp.isfile(file), FileNotFoundError(f"file number {i} not found: {file}")
            if file.lower().endswith(".fastq"):
                bin_classify.format = "fastq"
            elif file.lower().endswith(".fasta"):
                bin_classify.format = "fasta"
            else:
                raise NotImplementedError
            # setting time
            base_name = osp.basename(file)
            key = base_name
            t[key] = {}
            t[key]["start"] = perf_counter()

            logger.info(f"Opening fastq file ({i+1}/{len(list_fastq)}) {base_name}")
            # Binning
            if "bins" in db_type:
                ReadToBin.set_fastq_model_and_param(file, path_model, param, cores, k)
                ReadToBin.bin_reads()
                ReadToBin.sort_bins_by_sizes()
                t[key]["binning"] = perf_counter()
                t[key]["reads_nb"] = ReadToBin.NUMBER_BINNED

            fastq_classifier = MockCommunity(
                path_original_fastq=file, db_path=path_to_hash, db_type=db_type, folder_report=path_report,
                path_binned_fastq=ReadToBin.outputs, bin_nb=bin_nb, classifier_name=classifier, param=param, cores=cores)

            fastq_classifier.classify()
            t[key]["classify"] = perf_counter()
            t[key]["hashes"] = fastq_classifier.hash_files

        except Exception as e:
            logger.exception(e)
            logger.warning(f"script crashed for file: {file}")

    records = []
    for key in t.keys():
        if 'classify' not in t[key].keys():
            break
        if "binning" in t[key]:
            t_binning = time_to_hms(t[key]['start'], t[key]['binning'], short=True)
            t_classify = time_to_hms(t[key]['binning'], t[key]['classify'], short=True)
            t_total = time_to_hms(t[key]['start'], t[key]['classify'], short=True)
            hashes = t[key]["hashes"]
            h_size = sum([osp.getsize(f) for f in hashes.values()])

            logger.info(f"timings for file {key} / binning : {t_binning}, for {t[key]['reads_nb']} reads")
            logger.info(f"timings for file {key} / classify: {t_classify}, "
                        f"{len(hashes)} bins, total size of hashes loaded: {h_size/10**9:.2f} GB")
        else:
            t_binning = time_to_hms(t[key]['start'], t[key]['start'], short=True)
            t_classify = time_to_hms(t[key]['start'], t[key]['classify'], short=True)
            t_total = t_classify
            hashes = t[key]["hashes"]
            h_size = sum([osp.getsize(f) for f in hashes.values()])
            logger.info(f"timings for file {key} / classify: {time_to_hms(t[key]['start'], t[key]['classify'])}")

        # to CSV
        # todo: add precision / sensitivity / abundance
        row = (key, db_type, t_binning, t_classify, t_total, f"{h_size/10**9:.2f}GB", f"{len(hashes)}",
               path_database, osp.basename(path_database))
        records.append(row)

    # Timings and to csv
    with open(f_record, 'a', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerows(records)

    logger.info(f"Script ended, {len(t)} files processed")
    print()


bin_classify.classifiers = ('kraken2',)
bin_classify.cores = 1
bin_classify.format = "fastq"


def test_classification():
    """ Should have a toy data set that i can bin, classify, and check the results """
    # todo: toy data set to check if it works
    raise NotImplementedError

    
if __name__ == '__main__':
    parser = ArgumentParserWithDefaults(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('output_folder',      help='Folder for output reports', type=is_valid_directory)
    parser.add_argument('database',           help='Folder with the hash table for the classifier, name '
                                                   '"minikm_<param>" with sub-folders "RefSeq/<bins> '
                                                   'and "model_<name>.pkl" ')
    parser.add_argument('-c', '--classifier', help='choose which metagenomics classifier to use', metavar='',
                                              choices=bin_classify.classifiers, default=bin_classify.classifiers[0])
    parser.add_argument('-s', '--clf_settings', help="detailed settings, such as 'k35_l31_s7' for kraken2",
                                              metavar='', default='k35_l31_s7')
    parser.add_argument('-b', '--full_DB',  help='Choose to use the standard full database or the segmented one',
                                              default='bins', choices=('full', 'bins',), metavar='')
    parser.add_argument('-t', '--threads',    help='Number of threads', default=cpu_count(), type=int, metavar='')
    parser.add_argument('-i', '--input_fastq',help='List of input files in fastq format, space separated.',
                                              default=path_fastq_comm, type=is_valid_file, nargs="+", metavar='')
    parser.add_argument('-r', '--record',     help='Record the time spent for each run in CSV format',
                                              default="/home/ubuntu/classify_records.csv", type=str, metavar='')
    # parser.add_argument('-c', '--cores',         help='Number of cores', default=cpu_count(), metavar='')

    args = parser.parse_args()

    bin_classify(args.input_fastq, args.output_folder, args.database,
                 classifier=args.classifier, db_type=args.full_DB, cores=args.threads, f_record=args.record,
                 clf_settings=args.clf_settings)




