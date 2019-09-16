#!/usr/bin/env python3
# ######################################################################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-06-18
#
# ######################################################################################################################
# Synthesise reads from genomes .fna files.
#


import argparse
import gzip
import os
from os.path import getsize, join as path_join
import pandas as pd
import subprocess

from tqdm import tqdm

path_genomes   = "/home/sjriondet/Data/ncbi-2019-complete/refseq"
path_out_reads = "/home/sjriondet/Data/synthetic_reads"


def preprocess_genomes(folder):
    """ Read through all genomes files to have a summary of the available data """
    cache_file_name = "_genomes.pd"
    # meta_file_name = ""
    cache_file_path = path_join(folder, cache_file_name)

    if os.path.isfile(cache_file_path):
        return pd.read_pickle(cache_file_path)
    else:
        # todo create the summary file
        number_genomes = sum([len(files) for _, _, files in os.walk(path_genomes)])
        for dir_path, dir_names, files in tqdm(os.walk(path_genomes), total=number_genomes):
            for file in files:
                if file.endswith((".fna", ".fa", ".fastq", )):
                    file_path = path_join(dir_path, file)
                    file_size = getsize(file_path)
                    with gzip.open(file_path, 'rb') as f:
                        content = f.read().decode("utf-8")
                        content.count("\n")
                        # what is inside, example:
                        # b'>NC_002162.1 Ureaplasma parvum serovar 3 str. ATCC 700970, complete genome\nATGGCTAATAATTATCAAACTTTATATGATTCAGCAATAAAAAGGATTCCATACGATCTTATTTCTGATCAAGCTTATGC\nAATTCTACAAAATGCTAAAACTCATAAAGTTTGCGATGGTGTTT'

        df = pd.DataFrame()
        df.to_pickle(cache_file_path)
        return df


def reads_from_genomes(number_reads):
    """ Create synthetic reads from intact genomee """
    # todo: add a counter (like want x reads)

    number_genomes = sum([len(files) for _, _, files in os.walk(path_genomes)])
    reads


    # todo: add a counter (like want x reads)
    for dir_path, dir_names, files in tqdm(os.walk(path_genomes), total=number_genomes):
        print(path_out_reads)




read_synthesizers = {"genome": reads_from_genomes,
                     }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create reads for benchmarking')
    parser.add_argument('read_type', help='Select which kind of read you want', choices=read_synthesizers.keys())
    # parser.add_argument('-s', '--selected_tasks', type=lambda s: {int(item) for item in s.split(',')} )
    # parser.add_argument('-z', '--gzip', help='zgip input files', action='store_true')
    # parser.add_argument('-m', '--max', help='Maximum files to run', type=int, default=0)
    # parser.add_argument('-n', '--nanopore_custom', help='Enable Nanopore customed settings (Centrifuge ONLY). ',
    #                             type=lambda s: [item for item in s.split(',')][:2], )
    args = parser.parse_args()


