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
import subprocess
import os

from tqdm import tqdm

path_genomes   = "/home/sjriondet/Data/ncbi-2019-complete/refseq"
path_out_reads = "/home/sjriondet/Data/synthetic_reads"


def reads_from_genomes():

    genomes_count = sum([len(files) for r, d, files in os.walk(path_genomes)])

    for dir_path, dir_names, files in tqdm(os.walk(path_genomes), total=genomes_count):
        print(path_out_reads)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run multiple classifiers')
    parser.add_argument('folder', help='folder to check')
    parser.add_argument('-s', '--selected_tasks', type=lambda s: {int(item) for item in s.split(',')} )
    parser.add_argument('-z', '--gzip', help='zgip input files', action='store_true')
    parser.add_argument('-m', '--max', help='Maximum files to run', type=int, default=0)
    parser.add_argument('-n', '--nanopore_custom', help='Enable Nanopore customed settings (Centrifuge ONLY). ',
                                type=lambda s: [item for item in s.split(',')][:2], )
    args = parser.parse_args()
