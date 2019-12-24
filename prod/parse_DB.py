#!/usr/bin/env python3
"""
#############################################################################
Script to divide a Database of genomes (RefSeq) and split it into bins
according to the k-mer frequency of each genome. Needs a lot of disk space,
and RAM according to the largest genome to process.

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

import argparse
from copy import deepcopy
import logging
import os
import os.path as osp
import pandas as pd
import re
from time import time

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Import paths and constants for the whole project
from prod.util import PATHS, FilesInDir, is_valid_directory
from prod.bio import kmers_dic, ncbi, Read, seq_count_kmer


# #############################################################################
# https://docs.python.org/3/howto/logging-cookbook.html
# create formatter for the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# create file handler which logs even debug messages
fh = logging.FileHandler(PATHS.LOGS)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
ch.setFormatter(formatter)
# create logger with parse_DB.py and add the handlers to the logger
logger = logging.getLogger('parse_DB')
logger.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(ch)


class Genome:
    """ Genome from RefSeq. Methods to split it into plasmid/genome and into segments """
    categories = ["plasmid", "chloroplast", "scaffold", "contig",
                  "chromosome", "complete genome", "whole genome shotgun sequence", ]
    kmer_count_zeros = kmers_dic(4)

    def __init__(self, fna_file, kmer_file, taxon, window, k=4):
        self.path_fna = fna_file
        self.path_kmers  = kmer_file
        self.taxon       = taxon
        self.window      = window
        self.k           = k
        # records is a dict of SeqRecord
        self.records = {cat: [] for cat in self.categories}
        # self.splits  = {cat: [] for cat in self.categories}

    def load_genome(self):
        """ Loop through all chromosomes/plasmids/genomes/.. and parse them into the object """
        for record in SeqIO.parse(self.path_fna, "fasta"):
            for cat in self.categories:
                if cat in record.description:
                    self.records[cat].append(record)
                    break

    def yield_genome_split(self):
        """ Split a genome/plasmid into multiple windows, to count the k-mer
            or to create the .fna files for kraken2-build
        """
        for cat in self.categories:
            for record in self.records[cat]:

                full_seq = record.seq
                len_genome = len(full_seq)
                for start in range(0, len_genome - self.window, self.window):
                    end = min(start + self.window, len_genome - 1)
                    # Include the taxonomy id, start and end of the segment into the description
                    description = f"|kraken:taxid|{self.taxon}|s:{start}-e:{end-1}|{record.description}"
                    segment = SeqRecord(full_seq[start:end],
                                        record.id, record.name, description, record.dbxrefs,
                                        record.features, record.annotations, record.letter_annotations)
                    yield (segment, self.taxon, cat, start, end)

    def count_kmers_to_df(self):
        """ Take all splits, count the kmer distribution and save to the kmer folder as pandas DataFrame """
        for_csv = []
        for segment, taxon, cat, start, end in self.yield_genome_split():
            kmer_count = seq_count_kmer(str(segment.seq), deepcopy(self.kmer_count_zeros), k=self.k)
            for_csv.append((taxon, cat, start, end, segment.name, segment.description, self.path_fna,
                            *kmer_count.values() ))
        kmer_keys = list(self.kmer_count_zeros.keys())
        df = pd.DataFrame(for_csv, columns=["taxon", "category", "start", "end", "name", "description", "fna_path"]
                                           + kmer_keys)
        df.taxon       = df.taxon.astype('category')
        df.category    = df.category.astype('category')
        df.description = df.description.astype('category')
        df.to_pickle(self.path_kmers)


def idkwhat():
        if "write_fna":
            raise NotImplementedError


def scan_RefSeq_to_kmer_counts(scanning, folder_kmers, k=4, window=10000, stop=3, ):
    """ Scan through RefSeq, split genomes into windows, count their k-mer, save in similar structure
        Compatible with 2019 RefSeq format hopefully
    """
    start = time()

    # scanning folder Class set up:
    preset = {"root"  : [".taxon", ],
              "target": [".kmer_count.pd", ]}
    FilesInDir.set_defaults_parse_RefSeq(scanning, folder_kmers, preset=preset)
    FilesInDir.folder_kmer = folder_kmers

    logger.info("scanning through all genomes in refseq " + scanning)
    for i, fastq in enumerate(FilesInDir.tqdm_scan()):
        if osp.isfile(fastq.target_file[".kmer_count.pd"]):
            continue
        with open(fastq.taxon) as f:
            taxon = int(f.read())
        genome = Genome(fastq.abs_path, fastq.kmer_count_pd, taxon, window=window, k=k)
        genome.load_genome()
        genome.count_kmers_to_df()
        if i > stop:
            logger.warning("Early stop of the scanning")
            break

    # Ending
    elapsed_time = time() - start
    print(f"\n{FilesInDir.file_count_root} folders have been scanned\n"
          f"Took {elapsed_time:,.1f}s / {elapsed_time / 60:.1f}min  to complete. "
          f"{FilesInDir.file_count_root / elapsed_time:,.0f} genome/s")


def combine_genome_kmer_counts(folder_kmers, path_df):
    """ Combine single dataframes into one. Might need high memory """
    logger.info("loading all kmer frequencies into a single file from " + folder_kmers)
    dfs = []
    for file in os.scandir(folder_kmers):
        dfs.append(pd.read_pickle(file.path))
    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Saving dataframe to {path_df}")
    df.to_pickle(path_df)


def find_bins_DB(path_kmer_counts, n_parts=10):
    """ Given a database of genomes in fastq files, split it in n segments """
    raise NotImplementedError


def write_split_to_bins(path_df_bins, path_db_bins):
    pass


def kraken_build(path_db_bins, path_bins_hash):
    pass


def main(folder_database, folder_intermediate_files, n_clusters, cores):
    """ Pre-processing of RefSeq database to split genomes into windows, then count their k-mers
        Second part, load all the k-mer counts into one single Pandas dataframe
        Third train a clustering algorithm on the k-mer frequencies of these genomes' windows
        folder_database          : RefSeq root folder
        folder_intermediate_files: empty root folder to store kmer counts
    """
    scan_RefSeq_to_kmer_counts(folder_database, folder_intermediate_files, stop=3)
    path_kmer_counts = osp.join(folder_intermediate_files, "_all_counts.kmer.pd")
    combine_genome_kmer_counts(folder_intermediate_files, path_kmer_counts)
    # todo: find bins and write split genomes into bins
    path_df_bins = "todo"
    find_bins_DB(path_kmer_counts, path_df_bins)
    path_db_bins = "todo"
    write_split_to_bins(path_df_bins, path_db_bins)
    path_bins_hash = "todo"
    kraken_build(path_db_bins, path_bins_hash)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path_database', type=is_valid_directory,
                        help='Database root folder. Support format: RefSeq 2019.')
    parser.add_argument('path_intermediate_files', type=is_valid_directory,
                        help='Folder for the k-mer counts per genome and for clustering models')
    parser.add_argument('-c', '--clusters', default=10, type=int, help='Number of clusters to split the DB into')
    parser.add_argument('-t', '--threads',  default="10", help='Number of threads')
    args = parser.parse_args()

    main(folder_database=args.path_database, folder_intermediate_files=args.path_intermediate_files,
         n_clusters=args.clusters, cores=args.threads)
    print("Not implemented yet")
        





#
# Deprecated method
def kmer_pkl_path(kmer_folder, fna_path, taxo_ext="gff"):
    """ Legacy method, might not use it anymore in the future
        Return a file name based on the taxonomy id instead of the file name.
        We retrieve the taxo id from the .gff file.
        To avoid re-reading file, taxo id is stored into <bac>.taxon
    """
    assert taxo_ext in ("gbk", "gff"), "Only extensions .gbk and .gff are implemented"

    #     bacteria_name = os.path.split(os.path.split(fna_path)[0])[1]
    fna_name = os.path.split(os.path.splitext(fna_path)[0])[1]

    taxo = ""
    path_taxon = fna_path.replace(".fna", ".taxon")
    if os.path.isfile(path_taxon):
        with open(path_taxon) as f:
            taxo = f.read()

    if not str.isdigit(taxo):
        path_gbk = fna_path.replace(".fna", f".{taxo_ext}")
        assert os.path.isfile(path_gbk), f"{fna_path} DOESN'T have a .{taxo_ext} file ??"

        with open(path_gbk) as gbk:
            description = [next(gbk) for i in range(9)][-1]

        if taxo_ext == "gbk":
            identificator = 'db_xref="taxon:'
        elif taxo_ext == "gff":
            identificator = 'Taxonomy/Browser/wwwtax.cgi?id='
        taxo_start = description.find(identificator)
        taxo = description[taxo_start + len(identificator):
                           taxo_start + description[taxo_start:].find('\n')]

        assert 1 <= len(taxo) <= 8, f"The taxo id search failed, found an id of length {len(taxo)}, \n" \
            f"for the file: {path_gbk} \n" \
            f"found string : {taxo[:min(50, len(taxo))]} ..."

        with open(path_taxon, "w") as f:
            f.write(taxo)

    # todo: check here to retrieve species' name
    bacteria_name = ncbi.translate_to_names([taxo])[0]
    # path_taxo_names = "/home/ubuntu/Data/Segmentation/Kraken_10_clusters_V1/Kraken2_building/taxonomy/names.dmp"
    # taxo_table = pd.read_csv(path_taxo_names, sep="\t|\t")
    # query = taxo_table[(taxo_table.taxo == int(taxo)) & (taxo_table.class_name == "scientific name")]
    # assert query.shape[0] == 1, f"Found {query.shape[0]} matches for the scientific name of taxo {taxo}. "
    #                             f"Display the taxo table: \n {taxo_table[taxo_table.taxo == int(taxo)]}"
    # bacteria_name = query.name.iat[0]

    formatted_bacteria = re.sub('[^A-Za-z0-9]+', '_', bacteria_name)
    out_path = osp.join(PATHS.RefSeq_4mer_freq, kmer_folder, f"{taxo}__{fna_name}__{formatted_bacteria}.pd")
    return taxo, bacteria_name, fna_name, out_path


# Deprecated method
def count_all(folder_kmers, scanning=PATHS.RefSeq_DB, k=4, window=1000, stop=3, skip={}):
    start = time()
    n = 0
    nucleotides_counts = []
    dic_template = {"bacteria": "", "fna": "", "start": None, }
    dic_template.update(kmers_dic(k))

    # Looping through each family folder
    for genera in tqdm(os.scandir(scanning), desc="Genera", total=len(os.listdir(scanning))):
        if stop > 0 and n > stop:  # 5400
            break
        # Looping through each bacterial folder
        #         results = Parallel(n_jobs=n_cores)(delayed(extract_folder)(folder, dic_template, ) \
        #             for folder in tqdm(os.scandir(genera), desc="Species", total=len(os.listdir(genera)), leave=False))

        for folder in tqdm(os.scandir(genera), desc=genera.name, total=len(os.listdir(genera)), leave=False):
            if stop > 0 and n > stop:  # 5400
                break
            if genera.name in skip: continue

            if not os.path.isdir(folder): continue
            files = [f for f in os.scandir(folder) if f.name.endswith(".fna")
                     #                      and (f.name.startswith("NC_") or f.name.startswith("AC_"))
                     and "multiisoloate" not in f.path and "multispecies" not in f.path]
            if len(files) == 0: continue

            # Looping through each file for a single bacteria (multiple chromosomes or alternative genomes ?)
            bac_kmers = []
            for file_i in files:
                try:
                    # Check if already done
                    taxo, bacteria_name, fna_name, kmer_freq_path = \
                        kmer_pkl_path(folder_kmers, file_i.path, taxo_ext="gff")
                    if os.path.isfile(kmer_freq_path):
                        continue  # Already done for this folder

                    # Count
                    rec = read_fna(file_i)  # go through all files
                    dic_template["bacteria"] = bacteria_name
                    dic_template["fna"] = fna_name
                    dic_template["len_genome"] = len(rec)
                    success_n, kmer_counts = \
                        count_kmers(rec, dic_template, k, bacteria_name, fna_name, w=window)
                    succ_fail = "Success" if len(rec) - 3 == success_n else "Fail   "
                    #                     print(f"{succ_fail} -> Bacteria: {bacteria_name},\t file: {fna_name},\t len: {len(rec)}")
                    nucleotides_counts.append(success_n)

                    bac_kmers.extend(kmer_counts)
                except Exception as e:
                    print("type error: " + str(e))
                    #                     print(traceback.format_exc())
                    print(file_i.path)

            if len(bac_kmers) > 0:
                # Pandas
                df = to_pandas(bac_kmers, bacteria_name)
                # Save to a file
                df.to_pickle(kmer_freq_path)
                n += 1

    elapsed_time = time() - start
    total = sum(nucleotides_counts)
    print(f"\n{n} folders have been scanned\n"
          f"Took {elapsed_time:,.1f}s / {elapsed_time / 60:.1f}min  to complete. {total / elapsed_time:,.0f} bp/s")
    return nucleotides_counts
