#!/usr/bin/env python3
"""
#############################################################################
common resources for multiple scripts

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain RIONDET <sylvainriondet@gmail.com>
PLoT-ME: Pre-classification of Long-reads for Memory Efficient Taxonomic assignment
https://github.com/sylvain-ri/PLoT-ME
#############################################################################
"""
import argparse
from datetime import datetime
import logging
from multiprocessing import cpu_count
# from multiprocessing.pool import Pool
import numpy as np
import os
import os.path as osp
import pandas as pd
from pathlib import Path
import shutil
import subprocess
from tqdm import tqdm

from plot_me import LOGS

logger = logging.getLogger(__name__)


# #############################################################################
# https://docs.python.org/3/howto/logging-cookbook.html
# loggers = {}

def init_logger(logger_name='reads_binning', verbose=False):
    # if loggers.get(logger_name):
    #     return loggers.get(logger_name)
    # else:
    # create formatter for the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s --- %(message)s')
    # create file handler which logs even debug messages
    fh = logging.FileHandler(LOGS)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(formatter)
    # create logger with parse_DB.py and add the handlers to the logger
    new_logger = logging.getLogger(logger_name)
    new_logger.setLevel(logging.DEBUG)
    new_logger.addHandler(fh)
    new_logger.addHandler(ch)
    return new_logger


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


def create_path(path):
    """ Create the intermediate folders if not existing. """
    # consider that it's a file if the string after the "." is shorter than 4 character
    folder = osp.dirname(path) if "." in osp.basename(path) and len(osp.splitext(osp.basename(path))[1]) <= 4 else path
    if not osp.isdir(folder):
        logger.log(5, f"created folder {folder}")
        os.makedirs(folder, exist_ok=True)


def delete_folder_if_exists(path_dir):
    if osp.isdir(path_dir):
        logger.warning(f"Folder exists, DELETE IT ? (need to delete to redo a clean install): {path_dir}")
        user_in = input("y/[n]").lower()
        logger.debug(f"user entered: {user_in}")
        if 'y' in user_in:
            shutil.rmtree(path_dir)


def folder_today(path):
    s_today = f"{datetime.today()}"
    final_path = osp.join(path, s_today)
    if not osp.isdir(final_path):
        os.makedirs(final_path)
    return final_path


def f_size(path_or_size):
    """ If supplied a string, try to get the file size (otherwise size can be directly feed),
        then format the file size with MB/GB/TB and return it as a string """
    if isinstance(path_or_size, str):
        assert osp.isfile(path_or_size), FileNotFoundError(f"checking for file size, but file not found: {path_or_size}")
        size = osp.getsize(path_or_size)
    elif isinstance(path_or_size, (int, float)):
        assert path_or_size >= 0, ValueError(f"this function doesn't work with non positive value: {path_or_size}. supposed to be a file size")
        size = path_or_size
    else:
        raise NotImplementedError(f"Received neither a path (string) nor a number: {path_or_size}, can't return a file size")

    for threshold in f_size.splits.keys():
        if size > threshold:
            return f"{size/threshold:.2f} {f_size.splits[threshold]}"
        elif size == 0:
            return "0 B"
    raise


f_size.splits = {
    10**12: "TB",
    10**9 : "GB",
    10**6 : "MB",
    10**3 : "kB",
    1     : "B",
}


def bash_process(cmd, msg=""):
    """ execute a bash command (list of string), redirect stream into logger
        encoding=utf-8 to have text stream (somehow text=True not accepted by PyCharm),
        redirecting all stream to the Pipe, shell on for commands with bash syntax like wild cards
    """
    # https://docs.python.org/3/library/subprocess.html#subprocess.Popen
    if isinstance(cmd, str):
        shell = True
    else:
        shell = False
        assert isinstance(cmd, (list, tuple)), \
            TypeError(f"the input should be a list or tuple, but got type:{type(cmd)}, {cmd}")
    logger.info((msg if msg != "" else "launching bash command")
                + ": " + (cmd.split()[0] if shell else cmd[0]))
    logger.debug(cmd if shell else " ".join(cmd))

    # Combine stdout and stderr into the same stream, both as text (non binary)
    proc = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    for line in iter(proc.stdout.readline, ''):
        logger.debug(line.replace("\n", ""))
    # Check that the process ended successfully
    proc.wait(60*60*24)  # wait 24 hours max
    if proc.returncode == 123:
        logger.warning(f"Process {proc.pid} exited with exit status {proc.returncode}")
    elif proc.returncode != 0:
        logger.warning(f"Process {proc.pid} exited with exit status {proc.returncode}")
        raise ChildProcessError(f"see log file, bash command raised errors: " +
                                cmd if isinstance(cmd, str) else " ".join(cmd))


def div_z(n, d):
    return n / d if d else 0


def time_to_hms(start, end, fstring=True, short=False):
    assert start <= end, ArithmeticError(f"The start time is later than the end time: {start} > {end}")
    delay = int(end - start)
    m, s = divmod(delay, 60)
    h, m = divmod(m, 60)
    if short:
        return f"{h:d}:{m:02d}:{s:02d}"
    elif fstring:
        return f"{h:d} hours, {m:02d} minutes, {s:02d} seconds"
    else:
        return h, m, s


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    """ Customized Argparser to get both formatted docstring and defaults arguments
        https://stackoverflow.com/a/52025430/4767645 """
    def add_argument(self, *args, help=None, default=None, choices=None, **kwargs):
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
        if choices not in (None, [], ()) and args[0] != '-h':
            kwargs['default'] = default
            if help is not None:
                kwargs['help'] += " (choices: " + ", ".join(choices) + ")"
        super().add_argument(*args, **kwargs)


def scale_df_by_length(data, kmer_cols, k, w, single_row=False, ):
    """ Divide the kmer counts by the length of the segments, and multiply by the number kmer choices"""
    # todo: should scale by the actual number of columns (palindromes and reverse complemented k-mers)
    divider = w - k + 1
    ratio = 4**k / divider if divider > 1 else 4**k  # avoid divide by 0
    ratio = np.float32(ratio)
    if single_row:
        return data * ratio
    else:
        logger.info(f"Scaling the dataframe {data.shape}")
        logger.debug(f"{data.iloc[:5,:20]}")

        for col in tqdm(kmer_cols):
            data[col] *= ratio

        logger.debug(f"{data.iloc[:5,:20]}")


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


def import_cython_mod():
    """ Dirty way of importing cyt_ext """
    cyt_ext = ImportError
    cython_is_there = False
    try:
        try:
            from .cython_module import cyt_ext
        except:
            try:
                from plot_me.cython_module import cyt_ext
            except:
                from cython_module import cyt_ext
        cython_is_there = True
        logger.info("Cython has been imported")
    except ModuleNotFoundError:
        logger.warning("Module not found: 'from plot_me.cython_module import cython_module'")
    except ImportError:
        logger.warning("Import error 'from plot_me.cython_module import cython_module'")
    except Exception as e:
        logger.warning(e)
        logger.warning("\n ************************************************************ \n"
                       "Failed to import Cython extension, falling back to pure Python code. \n"
                       "Check the following: https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility \n"
                       "If this didn't solve your issue, Please consider raising an issue on github.")
    return cyt_ext, cython_is_there

# #############################################################################
# Save for programming
# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical message')



