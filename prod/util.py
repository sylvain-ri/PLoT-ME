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

util_logger = logging.getLogger('classify.util')
util_logger.debug('write messages')


# #############################################################################
# Paths
class ProjectPaths:
    def __init__(self):
        self.LOGS = f"/home/ubuntu/logs/{date.today()}.log"
        self.classifiers = ('kraken2', )
        self.DB = "/home/ubuntu/database/kraken2"
        self.KRAKEN2_DB = {
            "2015-10bins"  : f"{self.DB}/2015-10bins",
            "2015-standard": f"{self.DB}/2015-standard",
        }
        self.FOLDER_REPORTS = "/home/ubuntu/Data/Reports"
        
        self.models = "/home/ubuntu/Data/kmer_freq/4mer/V4"
        self.lda_model = f"{self.models}/LDA/lda_model_20_int.pd"
        self.kmeans_model = f"{self.models}/clustering/10_kmeans_2019-05-09_04-08.pkl"

# PATHS = ProjectPaths()


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
        

# #############################################################################
# Methods for nucleotides manipulations
nucleotides = "ACGT"

def kmers_dic(n, choice=nucleotides):
    return {a:0 for a in combinaisons(choice, n)}

def combinaisons(combi, n, instances=nucleotides):
    if n == 1:
        return combi
    else:
        return [f"{a}{n}" for a in combinaisons(combi, n-1) for n in instances]

def window(seq, window_size=4):
    """ Return a sliding window from a string """
    for i in range(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def seq_count_kmer(seq, kmer_count, k=4, ignore_N=True):
    """ Count all kmers, ignore kmers with N or other undecided nucleotides 
        seq: string nucleotide input
        kmer_count: dict with all combinations of nucleotides, initialized with zeros (kmer_template["AAAA"] = 0, kmer_count["AAAC"] = 0...)
        the new string hashing behaviour, BiopythonWarning: Using str(seq) to use the new behaviour
    """
    util_logger.info('counting kmers')
    wrong_base = "N"*k
    kmer_count[wrong_base] = 0
    
    try:
        for kmer in window(str(seq), k):
            try:
                kmer_count[kmer] += 1
            except:
                kmer_count[wrong_base] += 1
                
        if ignore_N:
            kmer_count.pop(wrong_base)
        return kmer_count
    except Exception as e:
        print("type error: " + str(e))
        print(traceback.format_exc())
        return kmer_count

# #############################################################################



# #############################################################################
# Save for programming
# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical message')



