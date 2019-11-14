#!/usr/bin/env python3
# ######################################################################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-06-18
#
# ######################################################################################################################
# Scrap 
# Synthesise reads from genomes .fna files.
#

# what is inside, example:
# b'>NC_002162.1 Ureaplasma parvum serovar 3 str. ATCC 700970, complete genome\nATGGCTAATAATTATCAAACTTTATATGATTCAGCAATAAAAAGGATTCCATACGATCTTATTTCTGATCAAGCTTATGC\nAATTCTACAAAATGCTAAAACTCATAAAGTTTGCGATGGTGTTT'


import argparse
import gzip
import os
from os.path import getsize, join as path_join
import pandas as pd
import subprocess
from tqdm import tqdm
from CONSTANTS import paths


# todo: check what in the NCBI database, count length and the line number of each genome, each chromosome, save it.
#  save the path, taxo id, lineage at various ranks
#  to be able to filter, select species at some rank, with some length of genome, ...

path_genomes   = "/home/ubuntu/Data/NCBI/20190704/refseq"
path_cache     = "/home/ubuntu/Data/NCBI"
path_out_reads = "/home/ubuntu/Data/Segmentation/Test-Data/Synthetic_from_Genomes"

meta_data = {"taxon", "file_path", "", "", "", }  # todo: add info for each field


def preprocess_genomes(folder, is_zip=True, update_summary=False):
    """ Read through all genomes files to have a summary of the available data """
    cache_file_name = "_genomes.pd"
    # meta_file_name = ""
    cache_file_path = path_join(path_cache, cache_file_name)
    fna_extensions = (".fna", ".fa", ".fastq", )
    if is_zip: fna_extensions = (ext + gz for ext in fna_extensions for gz in (".gz", ".gzip"))

    if os.path.isfile(cache_file_path) and update_summary is False:
        return pd.read_pickle(cache_file_path)
    
    else:
        # todo create the summary file
        summary_genomes = []
        
        number_genomes = sum([len(files) for _, _, files in os.walk(path_genomes)])
        for dir_path, dir_names, files in tqdm(os.walk(path_genomes), total=number_genomes):
            for file in files:
                    
                if file.endswith(fna_extensions):
                    file_path = path_join(dir_path, file)
                    file_size = getsize(file_path)
                    
                    if is_zip:
                        with gzip.open(file_path, 'rb') as f:
                            content = f.read().decode("utf-8")
                    else:
                        with open(file_path, 'rb') as f:
                            content = f.read()
                    summary_genomes.append(scrap_genomes(content))
                        
                        
        df = pd.DataFrame(summary_genomes, columns=["taxon", "file_path", "", "", "", ])
        df.to_pickle(cache_file_path)
        return df

def scrap_genomes(content):
    
    return ["something"]
    

def reads_from_genomes(number_reads):
    """ Create synthetic reads from intact genome """
    # todo: add a counter (like want x reads)

    number_genomes = sum([len(files) for _, _, files in os.walk(path_genomes)])



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


    
    
    
    
    
    
