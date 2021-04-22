#!/usr/bin/env python3
"""
#############################################################################
Project to divide a Database of Genomes (RefSeq) according to their
 k-mer frequency, for a lower RAM requirement of taxonomic classifiers.
This script pre-classifies reads/sequences from a fastq file, and calls a
 classifier (kraken2/centrifuge) with an index built on a subset of the
 entire database, for a lower memory need.
Clusters are defined by plot-me.pre-process / plot-me/parse_DB.py
The higher the number of clusters, the lower the memory requirement.

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain RIONDET <sylvainriondet@gmail.com>
PLoT-ME: Pre-classification of Long-reads for Memory Efficient Taxonomic assignment
https://github.com/sylvain-ri/PLoT-ME
#############################################################################
"""

import argparse
import csv
from datetime import datetime as dt
from glob import glob
import logging
from multiprocessing import cpu_count
# from multiprocessing.pool import Pool
import os
from os import path as osp
import pandas as pd
import pickle
import shutil
import subprocess
from pathlib import Path
from time import perf_counter
import re

import numpy as np
from Bio import SeqRecord, SeqIO
from tqdm import tqdm

# Import paths and constants for the whole project
from plot_me import RECORDS
from plot_me.tools import init_logger, scale_df_by_length, is_valid_directory, is_valid_file, create_path, \
    time_to_hms, f_size, bash_process, import_cython_mod
from plot_me.bio import kmers_dic, seq_count_kmer, combine_counts_forward_w_rc, n_dim_rc_combined, \
    codons_without_rev_comp

logger = logging.getLogger(__name__)
# If the total size of the reads, assigned to one bin, is below this percentage of the total fastq file, those reads are dropped

cython_is_there    = False
cyt_ext            = ImportError

THREADS            = 1
CLASSIFIERS        = (('kraken2', 'k35_l31_s7'),
                      ("centrifuge", ''))
K                  = None
BIN_NB             = None
DROP_BIN_THRESHOLD = -1  # by default, will be set as 1% / BIN_NB


def reads_in_file(file_path):
    """ Find the number of reads in a file.
        Count number of lines with bash wc -l and divide by 4 if fastq, otherwise by 2 (fasta) """
    return round(int(subprocess.check_output(["wc", "-l", file_path]).split()[0]) /
                     (4 if bin_classify.format == "fastq" else 2))


class Binner:
    """ Direct binner of reads, for dev purposes. """

    def __init__(self, _f_model, k=4, w=10**4):
        """ load a clustering model for the read binner """
        self.model     = None
        self.k         = k
        self.w         = w
        self.kmer_cols = []
        self._f_model  = _f_model

        self._load_model()
        self._set_cols()

    def _load_model(self):
        with open(self._f_model, 'rb') as f:
            self.model = pickle.load(f)

    def _set_cols(self):
        self.kmer_cols = codons_without_rev_comp(self.k)

    def scale(self, df):
        divider = df[self.kmer_cols].sum(1) - self.k + 1
        ratio = 4.0 ** self.k / divider
        return df[self.kmer_cols] * ratio

    def classify_df(self, df):
        assert isinstance(df, pd.DataFrame), TypeError("Only Pandas DataFrame data type supported at the moment")
        return self.model.predict(self.scale(df))


# #############################################################################
class ReadToBin(SeqRecord.SeqRecord):
    """ General Read. Wrapping SeqIO.Record """
    logger = logging.getLogger(__name__)
    KMER = {}  # kmers_dic(K)
    FASTQ_PATH = None
    FASTQ_BIN_FOLDER = None
    FILEBASE = ""
    MODEL = None
    PARAM = ""
    outputs = {}
    total_reads = 0
    file_has_been_binned = False
    NUMBER_BINNED = 0
    DIM_COMBINED = None

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
            if cython_is_there:
                self._kmer_count = cyt_ext.kmer_counter(str(self.seq), k=K, dictionary=True, combine=True, length=len(self.seq))
            else:
                self._kmer_count = combine_counts_forward_w_rc(seq_count_kmer(self.seq, self.KMER.copy(), K, ignore_N=ignore_N), k=K)
        return self._kmer_count

    @property
    def path_out(self, cluster=None):
        return f"{self.FASTQ_BIN_FOLDER}/{self.FILEBASE}.bin-{self.cluster if cluster is None else cluster}.fastq"

    def scale(self):
        self.logger.log(5, "scaling the read by it's length and k-mer")
        if cython_is_there:
            counts = cyt_ext.kmer_counter(str(self.seq), k=K, dictionary=False, combine=True, length=len(self.seq))
            cyt_ext.scale_counts(counts, K, len(self.seq))
            self.scaled = counts.reshape(-1, self.DIM_COMBINED)
        else:
            self.scaled = scale_df_by_length(np.fromiter(self.kmer_count.values(), dtype=np.float32)\
                                             .reshape(-1, self.DIM_COMBINED),
                                             None, k=K, w=len(self.seq), single_row=True)  # Put into 2D one row
        return self.scaled

    def find_bin(self):
        self.logger.log(5, 'finding bins for each read')
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
    def set_fastq_model_and_param(cls, path_fastq, path_model, param, force_binning):
        assert osp.isfile(path_fastq), FileNotFoundError(f"{path_fastq} cannot be found")
        # todo: load the parameter file from parse_DB.py instead of parsing string.... parameters_RefSeq_binning.txt
        cls.PARAM = param
        cls.FASTQ_PATH = path_fastq
        folder, file_base = osp.split(osp.splitext(path_fastq)[0])
        # output folder, will host one file for each bin
        cls.FASTQ_BIN_FOLDER = osp.join(folder, param)

        # Compute expected length
        cls.DIM_COMBINED = n_dim_rc_combined(K)

        cls.total_reads = reads_in_file(cls.FASTQ_PATH)

        # skip if reads already binned
        if osp.isdir(cls.FASTQ_BIN_FOLDER):
            total_binned_reads = 0
            if not force_binning:
                # Compute total reads count if it hasn't been forced
                for path in Path(cls.FASTQ_BIN_FOLDER).rglob("*bin-*.fastq"):
                    str_path = path.as_posix()
                    total_binned_reads += reads_in_file(str_path)
                    _, key, _ = re.split('.bin-|.fastq', str_path)
                    cls.outputs[int(key)] = str_path
                cls.logger.debug(f"A folder has been detected, and holds in total {total_binned_reads} reads, "
                                 f"compared to the {cls.total_reads} in the original fastq file.")

            if force_binning or cls.total_reads != total_binned_reads:
                last_modif = dt.fromtimestamp(osp.getmtime(cls.FASTQ_BIN_FOLDER))
                save_folder = f"{cls.FASTQ_BIN_FOLDER}_{last_modif:%Y-%m-%d_%H-%M}"
                cls.logger.warning(f"Folder existing, renaming to avoid losing files: {save_folder}")
                os.rename(cls.FASTQ_BIN_FOLDER, save_folder)
            else:
                # Flag up if read counts are equal, and no forcing to recount
                cls.file_has_been_binned = True
        create_path(cls.FASTQ_BIN_FOLDER)

        cls.FILEBASE = file_base
        if not path_model == "full":
            cls.KMER = kmers_dic(K)
            with open(path_model, 'rb') as f:
                cls.MODEL = pickle.load(f)

    @classmethod
    def bin_reads(cls):
        """ Bin all reads from provided file """
        # Skip binning if already done. Count total number of lines in each binned fastq
        if cls.file_has_been_binned:
            cls.logger.info(f"Fastq has already been binned, skipping reads binning: {cls.FASTQ_PATH}")
            return

        cls.logger.info(f"Binning the reads (count kmers, scale, find_bin, copy to file.bin-<cluster>.fastq")
        # todo: try to parallelize it, careful of file writing concurrency.
        #  Dask ? process to load and count kmers, single one for appending read to fastq ?
        # with Pool(cls.CORES) as pool:
        #     results = list(tqdm(pool.imap(pll_binning, SeqIO.parse(cls.FASTQ_PATH, "fasta"))))
        # counter = len(results)
        # TODO Cyt: use Cython file reader, k-mer counter and binner, fallback on Python otherwise
        counter = 0
        for record in tqdm(SeqIO.parse(cls.FASTQ_PATH, bin_classify.format), total=cls.total_reads,
                           desc="binning and copying reads to bins", leave=True, dynamic_ncols=True):
            counter += 1
            custom_read = ReadToBin(record)
            # custom_read.kmer_count
            custom_read.scale()
            custom_read.find_bin()
            custom_read.to_fastq()
        cls.logger.info(f"{counter} reads binned into bins: [" + ", ".join(map(str, sorted(cls.outputs.keys()))) + "]")
        cls.NUMBER_BINNED = counter
        return cls.outputs

    @classmethod
    def sort_bins_by_sizes_and_drop_smalls(cls):
        """ Sort the fastq bins by their size. drop_bins is the *percentage* below which a bin is ignored """
        bin_size = {}
        for f in os.scandir(cls.FASTQ_BIN_FOLDER):
            bin_nb = int(f.name.split('.')[1].split('-')[1])
            size = osp.getsize(f)
            bin_size[size] = bin_nb
        full_fastq_size = osp.getsize(cls.FASTQ_PATH)
        minimum_size = full_fastq_size * DROP_BIN_THRESHOLD / 100

        # make a copy first, then empty the dic, and rewrite it in the correct order
        fastq_outputs = ReadToBin.outputs.copy()
        ReadToBin.outputs = {}
        dropped_bins = []
        dropped_size = 0
        for size, bin_nb in sorted(bin_size.items(), reverse=True):
            if size > minimum_size:
                ReadToBin.outputs[bin_nb] = fastq_outputs[bin_nb]
            else:
                dropped_bins.append(bin_nb)
                dropped_size += size
                cls.logger.debug(f"Reads in bin {bin_nb} has a size of {f_size(size)}, and will be dropped "
                                 f"(less than {DROP_BIN_THRESHOLD}% of all binned reads {f_size(full_fastq_size)})")
        cls.logger.warning(f"Dropped bins {dropped_bins}, with total file size of {f_size(dropped_size)}. "
                        f"Lower parameter drop_bin_threshold to load all bins despite low number of reads in a bin.")
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
    
    def __init__(self, path_original_fastq, db_path, full_DB, folder_report, path_binned_fastq={},
                 classifier_name="kraken2", param="", clf_settings="default", dry_run=False, verbose=False):
        self.logger = logging.getLogger('classify.MockCommunity')

        assert osp.isfile(path_original_fastq), FileNotFoundError(f"Didn't find original fastq {path_original_fastq}")
        self.logger.debug(f"Path to hash files: db_path = {db_path}")

        self.path_original_fastq    = path_original_fastq

        self.folder, self.file_name = osp.split(osp.splitext(self.path_original_fastq)[0])
        self.path_binned_fastq      = path_binned_fastq              # {<bin i>: <path_file>}
        self.folder_report          = folder_report
        
        self.classifier_name = classifier_name
        self.db_path         = db_path    # location of the hash table for the classifier
        self.db_type         = "full" if full_DB else "bins"    # Either full or bins
        self.hash_size      = {}
        self.folder_out      = osp.join(self.folder_report, self.file_name)
        self.path_out        = osp.join(self.folder_out, f"{param}.{classifier_name}.{clf_settings}.{self.db_type}")

        self.dry_run         = dry_run
        self.verbose         = verbose
        self.cmd             = None

        # Initialization functions
        os.makedirs(self.folder_out, exist_ok=True)
        self.archive_previous_reports()

    @property
    def classifier(self):
        if self.classifier_name == "kraken2":
            return self.kraken2
        elif self.classifier_name == "centrifuge":
            return self.centrifuge
        else:
            NotImplementedError("This classifier hasn't been implemented")

    def archive_previous_reports(self):
        """ move existing reports to _archive """
        archive_folder = osp.join(self.folder_out, "_archive")
        self.logger.info("archiving previous reports into: " + archive_folder)
        os.makedirs(archive_folder, exist_ok=True)
        for file in glob(self.path_out + "*"):
            shutil.move(file, osp.join(archive_folder, osp.basename(file)))

    def classify(self):
        self.logger.info(f"Classifying reads with {self.db_type} setting")
        if "bins" in self.db_type:
            for bin_id in self.path_binned_fastq.keys():
                folder_hash = osp.join(self.db_path, f"{bin_id}")
                self.logger.debug(f"Path of fastq bin : {self.path_binned_fastq[bin_id]}")
                self.logger.debug(f"Path of folder of hash bin : {folder_hash}")
                self.classifier(self.path_binned_fastq[bin_id], folder_hash, arg=f"bin-{bin_id}")
            # todo: combine reports to Kraken2 format
        elif "full" in self.db_type:
            self.classifier(self.path_original_fastq, self.db_path, arg="full")
        else:
            NotImplementedError("The database choice is either full or bins")
                
    def centrifuge(self, fastq_input, folder_hash, arg="unknown"):
        """ Centrifuge calls
            https://ccb.jhu.edu/software/centrifuge/manual.shtml#command-line
        """
        hashes_file = [osp.join(folder_hash, f"cf_index.{i}.cf") for i in range(1, 4)]
        hash_root = osp.join(folder_hash, "cf_index")
        assert osp.isfile(hashes_file[0]), FileNotFoundError(f"Hash table not found ! {hash_root}*")
        self.hash_size[arg] = sum(osp.getsize(f) for f in hashes_file)
        self.logger.info(f'start to classify reads from file ({f_size(fastq_input)}) {fastq_input}')
        self.logger.info(f'with centrifuge, {arg}. hash table is ({f_size(self.hash_size[arg])}) {hash_root}*')
        out_path = f"{self.path_out}.{arg}" if self.db_type == "bins" else f"{self.path_out}"
        out_file = f"{out_path}.out"
        self.logger.info(f'output is {out_file}')
        self.cmd = [
            "centrifuge", "-x", hash_root, "-U", fastq_input,
            "-S", out_file, "--report-file", f"{out_path}.centrifuge-report.tsv",
            "--time", "--threads", f"{THREADS}",
        ]
        if self.dry_run:
            self.logger.debug(" ".join(self.cmd))
        else:
            bash_process(self.cmd, f"launching centrifuge classification on {fastq_input}")
            # Then do the kraken2 report
            cmd2 = ["centrifuge-kreport", "-x", hash_root, out_file, ">", f"{out_path}.report"]
            bash_process(" ".join(cmd2), f"launching centrifuge kreport on {fastq_input}")

    def kraken2(self, fastq_input, folder_hash, arg="unknown"):
        if "hash.k2d" in folder_hash: folder_hash = osp.dirname(folder_hash)
        hash_file = osp.join(folder_hash, "hash.k2d")
        assert osp.isfile(hash_file), FileNotFoundError(f"Hash table not found ! {hash_file}")
        self.hash_size[arg] = osp.getsize(hash_file)
        self.logger.info(f'start to classify reads from file ({f_size(fastq_input)}) {fastq_input}')
        self.logger.info(f'with kraken2, {arg}. hash table is ({f_size(hash_file)}) {hash_file}')
        formatted_out = f"{self.path_out}.{arg}" if self.db_type == "bins" else f"{self.path_out}"
        self.logger.info(f'output is {formatted_out}.out')
        self.cmd = [
            "kraken2", "--threads", f"{THREADS}",
            "--db", folder_hash,
            fastq_input,
            "--output", f"{formatted_out}.out",
            "--report", f"{formatted_out}.report",
        ]
        if self.dry_run:
            self.logger.debug(" ".join(self.cmd))
        else:
            bash_process(self.cmd, f"launching kraken2 classification on {fastq_input}")
            
    def kraken2_report_merging(self):
        self.logger.info('Merging kraken2 reports')
        # todo: merging reports to species level
        raise NotImplementedError()

    def report_to_csv(self):
        raise NotImplementedError()
        # todo: create .csv per species

    def __repr__(self):
        return f"Fastq file located at <{self.path_original_fastq}>, ready to be classified with " \
               f"{self.classifier_name} with the DB <{self.db_type}> located at {self.db_path}"
        

# #############################################################################
# Defaults and main method

def bin_classify(list_fastq, path_report, path_database, classifier, full_DB=False, threads=cpu_count(),
                 f_record="~/logs/classify_records.csv", clf_settings="", drop_bin_threshold=DROP_BIN_THRESHOLD,
                 skip_clas=False, force_binning=False, no_cython=False):
    """ Should load a file, do all the processing """
    _ = init_logger(__package__)  # initialize the global logger
    logger.info("\n*********************************************************************************************************")
    logger.info("**** Starting script **** \n ")
    global THREADS
    THREADS = threads
    if not no_cython:
        global cyt_ext, cython_is_there
        cyt_ext, cython_is_there = import_cython_mod()

    # preparing csv record file
    if not osp.isfile(f_record):
        os.makedirs(osp.dirname(f_record), exist_ok=True)
        with open(f_record, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            headers = ("FILE", "BINS_vs_FULL", "BINNING", "CLASSIFY", "TOTAL", "HASHES_SIZE", "NB_BINS", "HASH_PATH", "HASH_NAME")
            csv_writer.writerow(headers)

    logger.info("let's classify reads!")

    # Find the model
    global K, BIN_NB, DROP_BIN_THRESHOLD
    if full_DB:
        path_model = "full"
        K          = 0
        BIN_NB     = 1
        # clusterer, bin_nb, k, w, omitted = (None, 1, None, None, None)
        path_to_hash = path_database
        if "hash.k2d" in path_to_hash:
            path_to_hash = osp.dirname(path_to_hash)
        if "hash.k2d" not in os.listdir(path_to_hash):
            FileNotFoundError(f"hash.k2d not found in folder: {path_to_hash}")
    else:
        path_model = ""
        for file in os.scandir(path_database):
            if file.name.startswith("model.") and file.name.endswith(".pkl"):
                path_model = file.path
                break
        assert osp.isfile(path_model), FileNotFoundError(f"didn't find the ML model in {path_database}... {path_model}")

        # Parse the model name to find parameters:
        basename = path_model.split("/model.")[1]
        clusterer, bin_nb, k, w, omitted, _ = re.split('_b|_k|_s|_o|.pkl', basename)
        K      = int(k)
        BIN_NB = int(bin_nb)
        DROP_BIN_THRESHOLD = drop_bin_threshold if drop_bin_threshold != -1 else 1. / BIN_NB
        path_to_hash = osp.join(path_database, classifier, clf_settings)
        logger.debug(f"path_to_hash: {path_to_hash}")
        logger.debug(f"Found parameters: clusterer={clusterer}, bin number={BIN_NB}, k={K}, w={w}, omitted={omitted}")
        if cython_is_there:
            cyt_ext.set_verbosity(logging.INFO)
            cyt_ext.init_variables(K)

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
                raise NotImplementedError("The file is neither ending with .fasta nor with .fastq")
            # setting time
            base_name = osp.basename(file)
            key = base_name
            t[key] = {}
            t[key]["start"] = perf_counter()

            logger.info(f"Opening fastq file ({i+1}/{len(list_fastq)}) {f_size(file)}, {base_name}")
            # Binning
            if not full_DB:
                ReadToBin.set_fastq_model_and_param(file, path_model, param, force_binning)
                ReadToBin.bin_reads()
                ReadToBin.sort_bins_by_sizes_and_drop_smalls()
                t[key]["binning"] = perf_counter()
                t[key]["reads_nb"] = ReadToBin.NUMBER_BINNED

            if not skip_clas:
                fastq_classifier = MockCommunity(
                    path_original_fastq=file, db_path=path_to_hash, full_DB=full_DB, folder_report=path_report,
                    path_binned_fastq=ReadToBin.outputs, classifier_name=classifier, param=param)

                fastq_classifier.classify()
                t[key]["classify"] = perf_counter()
                t[key]["hashes"] = fastq_classifier.hash_size
            # todo: process reports to have one clean one

        except Exception as e:
            logger.exception(e)
            logger.error(f"script crashed for file: {file}")

    records = []
    for key in t.keys():
        if 'classify' not in t[key].keys():
            break
        if "binning" in t[key]:
            t_binning = time_to_hms(t[key]['start'], t[key]['binning'], short=True)
            t_classify = time_to_hms(t[key]['binning'], t[key]['classify'], short=True)
            t_total = time_to_hms(t[key]['start'], t[key]['classify'], short=True)
            hashes = t[key]["hashes"]
            h_size = sum(hashes.values())

            logger.info(f"timings for file {key} / binning : {t_binning}, for {t[key]['reads_nb']} reads")
            logger.info(f"timings for file {key} / classify: {t_classify}, "
                        f"{len(hashes)} bins, total size of hashes loaded: {f_size(h_size)}")
        else:
            t_binning = time_to_hms(t[key]['start'], t[key]['start'], short=True)
            t_classify = time_to_hms(t[key]['start'], t[key]['classify'], short=True)
            t_total = t_classify
            hashes = t[key]["hashes"]
            h_size = sum(hashes.values())
            logger.info(f"timings for file {key} / classify: {time_to_hms(t[key]['start'], t[key]['classify'])}")

        # to CSV
        # todo: add precision / sensitivity / abundance
        row = (key, "full" if full_DB else "bins", t_binning, t_classify, t_total, f"{h_size / 10 ** 9:.2f}GB", f"{len(hashes)}",
               path_database, osp.basename(path_database))
        records.append(row)

    # Timings and to csv
    with open(f_record, 'a', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerows(records)

    logger.info(f"Script ended, {len(t)} files processed \n")


bin_classify.format = "fastq"


def test_classification():
    """ Should have a toy data set that i can bin, classify, and check the results """
    # todo: toy data set to check if it works
    raise NotImplementedError

    
def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_plot_me',         help='Sub-folder generated by PLoT-ME, containing the pre-classifier '
                                                     'model (model*.pkl), as well as the classifier\'s hash tables. '
                                                     'In the form: `../PLoT-ME-data/k4_s10000/minikm_b10_k4_s10000_oAllRefSeq`. '
                                                     'If using the full_index, provide the path up to '
                                                     '`.../PLoT-ME-data/no-binning/o<omitted>`.')
    parser.add_argument('path_reports',         help='Folder for output reports', type=is_valid_directory)

    parser.add_argument('-i', '--input_fastq',  help='List of input files in fastq format, space separated.',
                                                default=[], type=is_valid_file, nargs="+", metavar='')
    parser.add_argument('-c', '--classifier',   help="classifier's name and its parameters, space separated. "
                                                     "Ex: '--classifier kraken k35_l31_s7', or '-c centrifuge'. "
                                                     "For unsupported classifiers, you can stop after "
                                                     "step 3, and build their index based on 'RefSeq_binned'",
                                                default=CLASSIFIERS[0], type=str, nargs="+", metavar='')
    parser.add_argument('-f', '--full_index',   help='Use the full index', action='store_true')
    parser.add_argument('-t', '--threads',      help='Number of threads (default=%(default)d)',
                                                default=cpu_count(), type=int, metavar='')
    parser.add_argument('-d', '--drop_bin_threshold', help='Drop fastq bins smaller than x percent of the initial '
                                                           'fastq. Helps to avoid loading hash tables for very few '
                                                           'reads (default = 1%% / <number of bins>)',
                                                default=DROP_BIN_THRESHOLD, type=float, metavar='')
    parser.add_argument('-r', '--record',       help='Record the time spent for each run in CSV format (default=%(default)s)',
                                                default=RECORDS, type=str, metavar='')
    parser.add_argument('--skip_classification',help='Skip the classification itself '
                                                     '(for benchmarking or to use other classifiers)',
                                                action='store_true')
    parser.add_argument('--force_binning',      help='If reads have already been binned, binning is skipped, unless '
                                                     'this flag is activated',
                                                action='store_true')
    parser.add_argument('--no_cython',          help='Disable Cython', action='store_true')

    args = parser.parse_args()
    logger.debug(f"Script {__file__} called with {args}")
    if len(args.classifier) == 1:
        args.classifier.append('')

    bin_classify(args.input_fastq, args.path_reports, args.path_plot_me,
                 classifier=args.classifier[0], full_DB=args.full_index, threads=args.threads, f_record=args.record,
                 drop_bin_threshold=args.drop_bin_threshold, skip_clas=args.skip_classification,
                 clf_settings=args.classifier[1], force_binning=args.force_binning, no_cython=args.no_cython)


if __name__ == '__main__':
    arg_parser()

