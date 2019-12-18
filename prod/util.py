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
# common resources for multiple scripts
# 

from datetime import date
import os
import os.path as osp
import logging

# import sys
# print('I am being imported by', sys._getframe(1).f_globals.get('__name__'))
# print(sys.argv[0])
from tqdm import tqdm

util_logger = logging.getLogger('classify.util')
util_logger.debug('write messages')


# #############################################################################
# File directory checking
def is_valid_directory(x):
    if osp.isdir(x):
        return x
    else:
        reply = input('Folder not found, would like to create it ? y/[n]')
        if 'y' in reply.lower():
            os.makedirs(x)
        else:
            util_logger.error('directory does not exist and has not been created ' + x)
            raise NotADirectoryError(f'The path is not a folder : {x}')
        return x


def is_valid_file(x):
    if osp.isfile(x):
        return x
    else:
        util_logger.error('file does not exist ' + x)
        raise FileNotFoundError(f'The path is not a file : {x}')


def folder_today(path):
    s_today = f"{date.today()}"
    final_path = osp.join(path, s_today)
    if not osp.isdir(final_path):
        os.makedirs(final_path)
    return final_path


def div_z(n, d):
    return n / d if d else 0


# #############################################################################
# Paths
class ProjectPaths:
    def __init__(self):
        self.data = "/home/ubuntu/Data"
        self.LOGS = f"/home/ubuntu/logs/{date.today()}.log"

        self.RefSeq_DB = f"{self.data}/NCBI/20190704/refseq"
        self.RefSeq_kmer_freq = f"{self.data}/kmer_freq"
        self.RefSeq_4mer_freq = f"{self.RefSeq_kmer_freq}/4mer"

        self.classifiers = ('kraken2', )
        self.classifier_DB = "/home/ubuntu/database/kraken2"
        self.kraken2_DB = {
            "2015-10bins"  : f"{self.classifier_DB}/2015-10bins",
            "2015-standard": f"{self.classifier_DB}/2015-standard",
        }
        self.folder_reports = f"{self.data}/Reports"
        
        self.models = f"{self.data}/kmer_freq/4mer/V4"
        self.lda_model = f"{self.models}/LDA/lda_model_20_int.pd"
        self.kmeans_model = f"{self.models}/clustering/10_kmeans_2019-05-09_04-08.pkl"


PATHS = ProjectPaths()


class FilesInDir:
    """ Provide the root path, then scan the entire folder for all matching files
        comes with file checking (extension, matching string, special methods like matching extension)
        applied for genome files
    """
    obj_id = 0
    file_count_root = None
    root_folder = ""
    target_folder = ""
    matching_ext = (".fastq", ".fq", ".fna")
    preset_expected_ext = {"root"  : [],
                           "target": []}
    _preset_genomes = {"root"  : [".taxon", ],
                       "target": [".kmer_count.pd", ]}

    def __init__(self, path):
        FilesInDir.obj_id += 1
        self.logger = logging.getLogger('parse_DB.Paths')

        self.abs_path = os.path.abspath(path)
        self.rel_path = osp.relpath(self.abs_path, self.root_folder)

        self.root_file   = {}
        self.target_file = {}
        self.check_root = True
        self.check_target = False
        self.presets()

    @property
    def kmer_count_path(self):
        kmer_path = osp.join(FilesInDir.target_folder, self.rel_path)
        kmer_dir = osp.split(kmer_path)[0]
        if not osp.isdir(kmer_dir):
            os.makedirs(kmer_dir)
        return kmer_path

    @staticmethod
    def create_target_path(path):
        folder = osp.split(path)[0]
        if not osp.isdir(folder):
            os.makedirs(folder)

    def add_tied_root_file(self, new_ext):
        """ file that is related, in the same root directory, with different extension """
        self.root_file[new_ext] = osp.splitext(self.abs_path)[0] + new_ext

    def add_tied_target_file(self, new_ext):
        """ file that is related, in the the target directory, with different extension. ensure the directory exists """
        self.target_file[new_ext] = osp.join(FilesInDir.target_folder,
                                             osp.splitext(self.rel_path)[0] + new_ext)
        self.create_target_path(self.target_file)

    def presets(self):
        """ associated files expected, in root and target directories """
        for v in self.preset_expected_ext["root"]:
            self.add_tied_root_file(v)
        for v in self.preset_expected_ext["target"]:
            self.add_tied_target_file(v)

    def file_matches_ext(self):
        """ does the folder contains the file we are looking for (=with these extensions) """
        return self.rel_path.lower().endswith(FilesInDir.matching_ext)

    def file_complies(self):
        """ Flags to not proceed """
        if not self.file_matches_ext():
            return False
        if self.check_root:
            for key in self.root_file.keys():
                if not osp.isfile(self.root_file[key]):
                    self.logger.warning(f"Related file {key} not found in root directory for {self}")
                    return False
        if self.check_target:
            for key in self.root_file.keys():
                if not osp.isfile(self.root_file[key]):
                    self.logger.warning(f"Related file {key} not found in target directory "
                                        f"({self.target_folder}) for {self}")
                    return False
        self.logger.info(f"Processing file {self}")
        return True

    @classmethod
    def set_defaults_parse_RefSeq(cls, root, target, preset=None):
        """ set where the root and target directories are (usually files are read from somewhere
            and out files somewhere else
            Creates expected files extensions. These extensions will be set in the object
        """
        cls.root_folder = root
        cls.target_folder = target
        cls.preset_expected_ext = cls._preset_genomes if preset is None else preset

    @classmethod
    def tqdm_scan(cls, folder=""):
        """ replicated os.walk, with total file count, for a folder (default root folder)
            yields a FileInDir object
        """
        if cls.file_count_root is None:
            cls.count_root_files(folder)
        for obj in tqdm(cls.walk_dir(folder), total=cls.file_count_root):
            yield obj


    @classmethod
    def walk_dir(cls, folder=""):
        """ Walk through every files in a directory (default root folder) and yield FileInDir """
        for dir_path, dirs, files in os.walk(cls.root_folder if folder == "" else folder):
            for filename in files:
                file = cls(os.path.join(dir_path, filename))
                if file.file_complies():
                    yield file

    @classmethod
    def count_root_files(cls, folder):
        file_count = 0
        for _ in tqdm(cls.walk_dir(folder)):
            file_count += 1
        cls.file_count_root = file_count

    def __repr__(self):
        return self.abs_path

    # @property
    # def taxon_path(self):
    #     return osp.splitext(self.abs_path)[0] + ".taxon"
    # @property
    # def abs_path(self):
    #     return osp.join(self.folder, self.file)
    # @property
    # def rel_path(self):
    #     return osp.relpath(self.abs_path, self.NCBI_path)
    # @property
    # def kmer_count_path(self):
    #     return osp.join(self.abs_path, self.NCBI_path)


# #############################################################################
# Save for programming
# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical message')



