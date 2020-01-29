#!/usr/bin/env python3
"""
#############################################################################
common resources for multiple scripts

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""
import argparse
import shutil
from datetime import datetime
import os
import os.path as osp
from multiprocessing import cpu_count
from multiprocessing.pool import Pool

import numpy as np
import pandas as pd
from pathlib import Path
import logging

# import sys
# print('I am being imported by', sys._getframe(1).f_globals.get('__name__'))
# print(sys.argv[0])
from tqdm import tqdm


# #############################################################################
# Paths
class ProjectPaths:
    def __init__(self):
        self.home = str(Path.home())
        self.data = f"{self.home}/Data"
        n = datetime.now()
        self.LOGS = f"{self.home}/logs/{n:%Y}-{n:%m}-{n:%d}_{n:%H}-{n:%M}.log"

        # self.RefSeq_DB = f"{self.data}/NCBI/20190704/refseq"
        # self.RefSeq_kmer_freq = f"{self.data}/kmer_freq"
        # self.RefSeq_4mer_freq = f"{self.RefSeq_kmer_freq}/4mer"

        # self.classifiers = ('kraken2', )
        # self.classifier_DB = f"{self.home}/database/kraken2"
        # self.kraken2_DB = {
        #     "2015-10bins"  : f"{self.classifier_DB}/2015-10bins",
        #     "2015-standard": f"{self.classifier_DB}/2015-standard",
        # }
        self.folder_reports = f"{self.data}/Reports"
        
        # self.models = f"{self.data}/kmer_freq/4mer/V4"
        # self.lda_model = f"{self.models}/LDA/lda_model_20_int.pd"
        # self.kmeans_model = f"{self.models}/clustering/10_kmeans_2019-05-09_04-08.pkl"


PATHS = ProjectPaths()


# #############################################################################
# https://docs.python.org/3/howto/logging-cookbook.html
def init_logger(logger_name='reads_binning', verbose=True):
    # create formatter for the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # create file handler which logs even debug messages
    # todo: find better log name: name of the script attached to the date ?
    fh = logging.FileHandler(PATHS.LOGS)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.DEBUG)
    ch.setFormatter(formatter)
    # create logger with parse_DB.py and add the handlers to the logger
    new_logger = logging.getLogger(logger_name)
    new_logger.setLevel(logging.DEBUG)
    new_logger.addHandler(fh)
    new_logger.addHandler(ch)
    return new_logger


logger = init_logger('tools')


# #############################################################################
# File directory checking
def is_valid_directory(x):
    if osp.isdir(x):
        return x
    else:
        reply = input(f'Folder not found, would like to create it ? y/[n] \n{x}')
        if 'y' in reply.lower():
            os.makedirs(x)
        else:
            logger.error('directory does not exist and has not been created ' + x)
            raise NotADirectoryError(f'The path is not a folder : {x}')
        return x


def is_valid_file(x):
    if osp.isfile(x):
        return x
    else:
        logger.error('file does not exist ' + x)
        raise FileNotFoundError(f'The path is not a file : {x}')


def create_path(path, with_filename=True):
    """ Create the intermediate folders if not existing. """
    folder = osp.split(path)[0] if with_filename else path
    if not osp.isdir(folder):
        logger.log(5, f"created folder {folder}")
        os.makedirs(folder)


def delete_folder_if_exists(path_dir):
    if osp.isdir(path_dir):
        logger.warning(f"Folder exists, DELETE IT ? (need to delete to redo a clean install): {path_dir}")
        user_in = input("y/[n]").lower()
        logger.debug(f"user entered: {user_in}")
        if 'y' in user_in:
            shutil.rmtree(path_dir)


def folder_today(path):
    s_today = f"{date.today()}"
    final_path = osp.join(path, s_today)
    if not osp.isdir(final_path):
        os.makedirs(final_path)
    return final_path


def div_z(n, d):
    return n / d if d else 0


def time_to_h_m_s(start, end, fstring=True):
    assert start < end, ArithmeticError(f"The start time is later than the end time: {start} > {end}")
    delay = int(end - start)
    m, s = divmod(delay, 60)
    h, m = divmod(m, 60)
    if fstring:
        return f"{h:d} hours, {m:02d} minutes, {s:02d} seconds"
    else:
        return h, m, s


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    """ Customized Argparser to get both formatted docstring and defaults arguments
        https://stackoverflow.com/a/52025430/4767645 """
    def add_argument(self, *args, help=None, default=None, **kwargs):
        if help is not None:
            kwargs['help'] = help
        if default not in (None, '') and args[0] != '-h':
            kwargs['default'] = default
            if help is not None:
                if default in (None, ''):
                    pass  # No default value to add
                if isinstance(default, list) or isinstance(default, tuple):
                    formatted = " ".join(default)
                    kwargs['help'] += f' ({type(default).__name__} - default: "{formatted})"'
                else:
                    kwargs['help'] += f' ({type(default).__name__} - default: {default} )'
        super().add_argument(*args, **kwargs)


def pll_scaling(serie):
    serie = pd.to_numeric(serie, downcast='float')
    serie *= pll_scaling.ratio
    return serie


pll_scaling.ratio = 0


def scale_df_by_length(data, kmer_cols, k, w, single_row=False, cores=cpu_count()):
    """ Divide the kmer counts by the length of the segments, and multiply by the number kmer choices"""
    ratio = 4**k / (w - k + 1)
    ratio = np.float32(ratio)
    if single_row:
        return data * ratio
    else:
        logger.info(f"Scaling the dataframe {data.shape}, converting to float32")
        logger.debug(f"{data}")

        pll_scaling.ratio = ratio
        with Pool(cores) as pool:  # file copy don't need many cores (main.cores)
            results = list(tqdm(pool.imap(pll_scaling, (data.loc[:, col] for col in kmer_cols)),
                                total=len(kmer_cols)))
        logger.debug(f"{data}")
        logger.debug(f"results len{len(results)}, {results[0]}")
        logger.debug(f"dataframe has been scaled {data.shape}")
        logger.debug(f"{results}")
        logger.debug(f"{results[0]}")
        logger.debug(f"{type(results[0])}")
        logger.debug(f"{len(results[0])}")
        logger.debug(f"{results[0].shape}")

        # for col in tqdm(kmer_cols):
        #     data.loc[:, col] = pd.to_numeric(data.loc[:, col], downcast='float')
        #     data.loc[:, col] *= ratio
        # logger.debug(f"{data}")
        # logger.debug(f"results, {results[0]}")
        # logger.debug(f"dataframe has been scaled {data.shape}")


class ScanFolder:
    """ Set class attributes, root & target folder, extensions to find and create
        tqdm scan the folder and create abs, rel, target path
    """
    obj_id        = 0
    folder_root   = ""
    folder_target = ""
    count_files   = None
    ext_find      = ()
    ext_check     = ""
    ext_create    = ""
    skip_folders  = ()

    def __init__(self, path):
        ScanFolder.obj_id += 1
        self.logger = logging.getLogger('tools.ScanFolder')

        self.path_abs      = os.path.abspath(path)
        self.path_rel      = osp.relpath(self.path_abs, self.folder_root)
        self.base          = osp.splitext(osp.split(self.path_abs)[1])[0]

    @property
    def path_check(self):
        """ Check if a file in the same folder, but different extension, is also in the same folder """
        assert self.ext_check != "", logger.error(f"No file extension provided to check files "
                                                  f"(define with ScanFolder.ext_check")
        return osp.splitext(self.path_abs)[0] + self.ext_check

    @property
    def path_target(self):
        if ScanFolder.folder_root == "":
            self.logger.warning("no root folder, set it with ScanFolder.folder_root = <path>")
            return ""
        elif ScanFolder.ext_create == "":
            self.logger.warning("no extension specified for the target file name")
            return ""
        else:
            path_to_target = osp.join(ScanFolder.folder_target, self.path_rel)
            res = osp.splitext(path_to_target)[0] + ScanFolder.ext_create
            create_path(res)
            return res

    def file_matches_ext(self):
        """ does the folder contains the file we are looking for (=with these extensions) """
        return self.path_rel.lower().endswith(self.ext_find)

    def file_complies(self, log=True):
        """ Find files with the extension to find, check if related file (check) """
        if not self.file_matches_ext():
            return False
        if self.ext_check != "" and not osp.isfile(self.path_check):
            self.logger.warning(f"Related file with extension {self.ext_check} not found in root directory for {self}")
            return False
        if log:  self.logger.log(5, f"file complies {self}")
        return True

    @classmethod
    def set_folder_scan_options(cls, scanning="", target="", ext_find=(), ext_check="", ext_create="", skip_folders=()):
        """ Set the options to scan a folder, filter files to find, files to check, and create the target path """
        assert osp.isdir(scanning), logger.error(f"the provided path to scan is not a directory {scanning}")
        assert target == "" or osp.isdir(target), logger.error(f"the provided path as target is not a directory {target}")
        cls.folder_root   = scanning
        cls.folder_target = target
        cls.ext_find      = ext_find
        cls.ext_check     = ext_check
        cls.ext_create    = ext_create
        cls.skip_folders  = skip_folders

    @classmethod
    def tqdm_scan(cls, folder="", with_tqdm=True):
        """ replicated os.walk, with total file count, for a folder (default root folder)
            yields a ScanFolder object
        """
        if folder != "":
            cls.folder_root = folder
        assert osp.isdir(cls.folder_root), logger.error(f"the provided path to scan is not a directory {cls.folder_root}")

        n = 0
        if with_tqdm:
            if cls.count_files is None:
                cls.count_root_files()
            logger.info(f"Yielding the {cls.count_files} files found in folder {cls.folder_root}")
            for obj in tqdm(cls.walk_dir(log=False), total=cls.count_files):
                n += 1
                yield obj
        else:
            for obj in cls.walk_dir(log=False):
                n += 1
                yield obj
        logger.debug(f"{n} have been processed")

    @classmethod
    def walk_dir(cls, log=True):
        """ Walk through every files in a directory (default root folder) and yield FileInDir """
        for dir_path, dirs, files in os.walk(cls.folder_root):
            # Skip folders
            rel_path = osp.relpath(dir_path, cls.folder_root)
            if any((name_to_skip in rel_path for name_to_skip in cls.skip_folders)):
                logger.debug(f"omitting folder {rel_path}")
                continue

            for filename in files:
                file = ScanFolder(os.path.join(dir_path, filename))
                if file.file_complies(log):
                    yield file

    @classmethod
    def count_root_files(cls):
        logger.debug(f"counting matching files in {cls.folder_root}")
        file_count = 0
        for _ in tqdm(cls.walk_dir()):
            file_count += 1
        cls.count_files = file_count
        return file_count

    def __repr__(self):
        return self.path_abs


# #############################################################################
# Save for programming
# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical message')



