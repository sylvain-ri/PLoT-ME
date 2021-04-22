#!/usr/bin/env python3
"""
#############################################################################
Project to divide a Database of Genomes (RefSeq) according to their
 k-mer frequency, for a lower RAM requirement of taxonomic classifiers
Needs a lot of disk space, and RAM according to the largest genome to process.

*** STEPS ***
0 -> Scan the given RefSeq, count kmer frequencies per segment for each genome
1 -> Combine these counts into a single file
2 -> Scale the values by the length of segments and combination of kmers,
     and apply a clustering algorithm (KMean, mini batch KMeans)
     to find the cluster association of each segment (RAM intensive)
3 -> Copy these segments of genomes into bins (DISK intensive)
4 -> kraken2-build --add-to-library
5 -> kraken2-build --build

*** FULL Index ***
-f/--full_index -> Building the hash for the full refseq, for comparison
  builds the index for the same classifier, same database, but without PLoT-ME.
  Good for comparison and benchmarking. --skip_existing and --early
  then refers to library preparation and index building
  (kraken2 --add_library and --build respectively).

Using large k (5+) and small s (10000-) yield very large kmer counts, costing
 high amounts of RAM (esp. when combining all kmer counts together,
 RAM needs to reach ~30GB or more). Currently working on this (July 2020)


#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain RIONDET <sylvainriondet@gmail.com>
PLoT-ME: Pre-classification of Long-reads for Memory Efficient Taxonomic assignment
https://github.com/sylvain-ri/PLoT-ME
#############################################################################
"""

import argparse
import logging
from glob import glob
import shutil
from copy import deepcopy
from itertools import islice
from multiprocessing import cpu_count, Pool
from pathlib import Path

from numpy import float32
import os
import os.path as osp
import pandas as pd
import pickle
import re
from time import perf_counter
import traceback

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from sklearn.cluster import KMeans, MiniBatchKMeans
# from sklearn.decomposition import PCA

from tqdm import tqdm

# Import paths and constants for the whole project
from plot_me import LOGS
from plot_me.tools import ScanFolder, is_valid_directory, init_logger, create_path, scale_df_by_length, \
    time_to_hms, delete_folder_if_exists, bash_process, f_size, import_cython_mod
from plot_me.bio import kmers_dic, ncbi, seq_count_kmer, combinations, combine_counts_forward_w_rc, \
    n_dim_rc_combined, codons_without_rev_comp


logger = logging.getLogger(__name__)
CLASSIFIERS     = (('kraken2', 'k', '35', 'l', '31', 's', '7'),
                   ("centrifuge", ))

# UGLY but necessary to allow import as raw python code as well as when compiled as a wheel/package
cython_is_there = False
cyt_ext = ImportError


class Genome:
    """ Genome from RefSeq. Methods to split it into plasmid/genome and into segments
        SET K BEFORE ANY INSTANCE IS CREATED, with set_k_kmers()
    """
    categories = ["plasmid", "chloroplast", "scaffold", "contig",
                  "chromosome", "complete genome", "whole genome shotgun sequence", ]
    K = 0
    col_kmers = []
    col_types = {}
    kmer_count_zeros = {}

    def __init__(self, fna_file, taxon, window_size, k=-1):
        logger.log(0, "Created genome object")
        self.path_fna    = fna_file
        self.taxon       = taxon
        self.window_size = window_size
        self.k           = Genome.K if k < 0 else k
        # records is a dict of SeqRecord
        self.records = {cat: {} for cat in self.categories}
        # self.splits  = {cat: [] for cat in self.categories}

    def __repr__(self):
        return f"Genome object from {osp.split(self.path_fna)[1]}"

    def load_genome(self):
        """ Loop through all chromosomes/plasmids/genomes/.. and parse them into the object
            Split them into various categories (plasmids, genomes, ...)
        """
        logger.debug(f"loading genome {f_size(self.path_fna)} {self.path_fna}")
        for record in SeqIO.parse(self.path_fna, "fasta"):
            for cat in self.categories:
                if cat in record.description:
                    self.records[cat][record.name] = record
                    break

    def yield_genome_split(self):
        """ Split a genome/plasmid into multiple segments, to count the k-mer
            or to create the .fna files for kraken2-build
        """
        for cat in self.categories:
            for record in self.records[cat].values():

                full_seq = record.seq
                len_genome = len(full_seq)
                for start in range(0, len_genome - self.window_size, self.window_size):
                    end = min(start + self.window_size, len_genome - 1)
                    # Include the taxonomy id, start and end of the segment into the description
                    description = f"|kraken:taxid|{self.taxon}|s:{start}-e:{end-1}|{record.description}"
                    segment = SeqRecord(full_seq[start:end],
                                        record.id, record.name, description, record.dbxrefs,
                                        record.features, record.annotations, record.letter_annotations)
                    yield (segment, self.taxon, cat, start, end)

    def count_kmers_to_df(self, path_kmers):
        """ Take all splits, count the kmer distribution and save to the kmer folder as pandas DataFrame """
        for_csv = []
        for segment, taxon, cat, start, end in self.yield_genome_split():
            if cython_is_there:
                kmer_count = cyt_ext.kmer_counter(str(segment.seq), k=self.k, dictionary=True, combine=True)
            else:
                kmer_count = combine_counts_forward_w_rc(seq_count_kmer(str(segment.seq), deepcopy(self.kmer_count_zeros), k=self.k),
                                                         k=self.k)
            for_csv.append((taxon, cat, start, end, segment.name, segment.description, self.path_fna,
                            *kmer_count.values() ))
        # kmer_keys = list(self.kmer_count_zeros.keys())
        df = pd.DataFrame(for_csv, columns=main.cols_types)
        df.taxon       = df.taxon.astype('category')
        df.category    = df.category.astype('category')
        df.name        = df.name.astype('category')
        df.fna_path    = df.fna_path.astype('category')
        for col in self.col_kmers:
            df[col] = df[col].astype("uint16")   # todo: not using float32 ??
        df.to_pickle(path_kmers)
        logger.debug(f"saved kmer count to {path_kmers}")

    @classmethod
    def initialize_set_k_mers(cls, k, combine_rc=True):
        cls.K = k
        cls.col_kmers = codons_without_rev_comp(k) if combine_rc else combinations(k)
        cls.kmer_count_zeros = kmers_dic(k)


def create_n_folders(path, n, delete_existing=False):
    """ Create the sub-folders of bins from 0/ to n/ """
    logger.debug(f"creates {n} folder under {path}")
    for i in range(n):
        new_path = osp.join(path, str(i))
        if delete_existing and osp.isdir(new_path):
            shutil.rmtree(new_path)
        create_path(new_path)


def add_file_with_parameters(folder, add_description=""):
    """ Add a file in the folder to remind which parameters were used, which folder were omitted """
    path = osp.join(folder, "parameters_RefSeq_binning.txt")
    # todo: update for identifier: value kind of. Need kraken k=35, reverse complement, comments to describe settings
    #  signature inside.
    with open(path, 'w') as f:
        f.write(f"script = {__file__} \n"
                f"From RefSeq located at: {main.folder_database} \n"
                f"k={main.k}, w={main.w} (segments size), \n"
                f"folders *containing* these strings have been omitted: " + ", ".join(main.omit_folders) + ". \n"
                f"{add_description}")


# Decorator for all these steps
def check_step(func):
    """ Decorator to print steps and check if results have already been computed
        Need the second argument to be the output file/folder: it will check if the file exists / folder isn't empty
    """
    def wrapper(*args, **kwargs):
        # Check arguments for debugging
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ",\t".join(args_repr + kwargs_repr)
        to_check = args[1]  # Output path for file or folder, will check if the output already exists

        # First check if skip has been allowed,
        if check_step.can_skip[check_step.step_nb] == "1" and \
                (osp.isfile(to_check)                                 # and there's already a file
                 or (osp.isdir(to_check) and os.listdir(to_check))):  # or there's a folder, not empty
            logger.info(f"Step {check_step.step_nb} SKIPPING, function \t{func.__name__}({signature}, "
                        f"Output has already been generated.")
            result = None

        # If limiting steps to run, stop it
        elif check_step.step_nb > check_step.early_stop:
            logger.info(f"Step {check_step.step_nb} EARLY STOP, no run for \t{func.__name__}({signature})")
            result = None

        else:
            # Time measurement
            start_time = perf_counter()
            logger.info(f"Step {check_step.step_nb} START, function \t{func.__name__}({signature})")
            result = func(*args, **kwargs)
            # print time spent
            logger.info(f"Step {check_step.step_nb} END, {time_to_hms(start_time, perf_counter())}, "
                        f"function {func.__name__}")

        # Step counter
        check_step.step_nb += 1
        check_step.timings.append(perf_counter())  # log time for each step
        return result
    return wrapper


check_step.timings    = []
check_step.step_nb    = 0         # For decorator to know which steps has been done
check_step.early_stop = -1        # Last step to run, later steps are not ran. Only display arguments
check_step.can_skip   = "111110"  # By default skip step that has been started, except fot kraken2 build (hard to check)


def parallel_kmer_counting(fastq, ):
    if osp.isfile(fastq.path_target):
        logger.debug(f"File already existing, skipping ({fastq.path_target})")
        return
    else:
        logger.debug(f"Processing: {fastq.path_abs}")
    with open(fastq.path_check) as f:
        taxon = int(f.read())
    genome = Genome(fastq.path_abs, taxon, window_size=main.w, k=main.k)
    genome.load_genome()
    genome.count_kmers_to_df(fastq.path_target)


@check_step
def scan_RefSeq_kmer_counts(scanning, folder_kmers, stop=-1):
    """ Scan through RefSeq, split genomes into segments, count their k-mer, save in similar structure
        Compatible with 2019 RefSeq format hopefully
    """
    create_path(folder_kmers)
    # scanning folder Class set up:
    # todo: change the kmer_count into the k_s_ notation
    ScanFolder.set_folder_scan_options(scanning=scanning, target=folder_kmers,
                                       ext_find=(".fastq", ".fq", ".fna"), ext_check=".taxon",
                                       ext_create=f".{main.k}mer_count.pd", skip_folders=main.omit_folders)

    logger.info("scanning through all genomes in refseq to count kmer distributions " + scanning)

    # Count in parallel. islice() to take a part of an iterable
    Genome.initialize_set_k_mers(main.k)
    with Pool(main.cores) as pool:
        results = list(tqdm(pool.imap(parallel_kmer_counting, islice(ScanFolder.tqdm_scan(with_tqdm=False),
                                                                     stop if stop>0 else None)),
                            total=ScanFolder.count_root_files(), dynamic_ncols=True))

    logger.info(f"{len(results)} genomes have been scanned and kmer counted.")


@check_step
def combine_genome_kmer_counts(folder_kmers, path_df):
    """ DEPRECATED Combine single dataframes into a single Dataframe. Might need high memory """
    logger.info("loading all kmer frequencies into a single file from " + folder_kmers)
    dfs = []
    added = 0
    ScanFolder.set_folder_scan_options(scanning=folder_kmers, target="", ext_find=(f".{main.k}mer_count.pd", ),
                                       ext_check="", ext_create="", skip_folders=main.omit_folders)
    for file in ScanFolder.tqdm_scan():
        dfs.append(pd.read_pickle(file.path_abs))
        added += 1
    logger.info(f"{added} {main.k}-mer distributions have been added. now concatenating")
    df = pd.concat(dfs, ignore_index=True)
    # Need to set again as categories
    df.taxon       = df.taxon.astype('category')
    df.category    = df.category.astype('category')
    df.name        = df.name.astype('category')
    df.fna_path    = df.fna_path.astype('category')
    # Save output file
    df.to_pickle(path_df)
    logger.info(f"Combined file of all kmer counts ({osp.getsize(path_df)/10**9:.2f} GB) save at: {path_df}")


@check_step
def append_genome_kmer_counts(folder_kmers, path_df):
    """ Combine single dataframes into one. Might need high memory """
    logger.info(f"Appending all kmer frequencies from {folder_kmers} into a single file {path_df}")
    added = 0
    ScanFolder.set_folder_scan_options(scanning=folder_kmers, target="", ext_find=(f".{main.k}mer_count.pd", ),
                                       ext_check="", ext_create="", skip_folders=main.omit_folders)
    # Append all the df. Don't write the index. Write the header only for the first frame
    for file in ScanFolder.tqdm_scan():
        if added == 0:
            pd.read_pickle(file.path_abs).to_csv(path_df, mode='w', index=False, header=True)
        else:
            pd.read_pickle(file.path_abs).to_csv(path_df, mode='a', index=False, header=False)
        added += 1
    logger.info(f"Combined file of {added} {main.k}-mer counts ({osp.getsize(path_df)/10**9:.2f} GB) save at {path_df}")


@check_step
def clustering_segments(path_kmer_counts, output_pred, path_model, n_clusters, model_name="minikm", combined_rc=True):
    """ Given a database of segments of genomes in fastq files, split it in n clusters/bins """
    assert model_name in clustering_segments.models, f"model {model_name} is not implemented"
    # Paths
    create_path(output_pred)
    create_path(path_model)
    k = main.k
    w = main.w
    nb_col = n_dim_rc_combined(k) if combined_rc else 4**k

    # https://www.codementor.io/@guidotournois/4-strategies-to-deal-with-large-datasets-using-pandas-qdw3an95k
    # filename = "data.csv"
    # n = sum(1 for line in open(filename)) - 1  # Calculate number of rows in file
    # s = n // 10  # sample size of 10%
    # skip = sorted(random.sample(range(1, n + 1), n - s))  # n+1 to compensate for header
    # df = pandas.read_csv(filename, skiprows=skip)

    path_pkl_kmer_counts = path_kmer_counts.replace(".csv", ".pd")
    if osp.isfile(path_pkl_kmer_counts):
        logger.info(f"Clustering the genomes' segments into {n_clusters} bins. Loading combined kmer counts "
                    f"(file size: {osp.getsize(path_pkl_kmer_counts)/10**9:.2f} GB) ...")
        df = pd.read_pickle(path_pkl_kmer_counts)
    else:
        logger.info(f"Clustering the genomes' segments into {n_clusters} bins. Loading combined kmer counts "
                    f"(file size: {osp.getsize(path_kmer_counts)/10**9:.2f} GB) ...")
        df = pd.read_csv(path_kmer_counts, dtype=main.cols_types)
        logger.info(f"save pickle copy for faster loading {path_pkl_kmer_counts}")
        # Soon DEPRECATED, set my the loading type
        # (Need to set again as categories)
        df.category    = df.category.astype('category')
        df.name        = df.name.astype('category')
        df.fna_path    = df.fna_path.astype('category')
        df.description = df.description.astype('category')
        df.to_pickle(path_pkl_kmer_counts)

    cols_kmers = df.columns[-nb_col:]
    cols_spe = df.columns[:-nb_col]
    logger.debug(f"cols_kmers={cols_kmers[:5]} {cols_kmers[-5:]}")

    # ## 1 ## Scaling by length and kmers
    df_mem = df.memory_usage(deep=False).sum()
    logger.info(f"Kmer counts loaded, scaling the values to the length of the segments. "
                f"DataFrame size: {df_mem/10**9:.2f} GB - shape: {df.shape}")

    # todo: save intermediate data
    scale_df_by_length(df, cols_kmers, k, w)

    # ## 2 ## Could add PCA

    # Model learning
    logger.info(f"Data takes {df_mem/10**9:.2f} GB. Training {model_name}...")
    if model_name == "kmeans":
        ml_model = KMeans(n_clusters=n_clusters, n_jobs=main.cores, random_state=3)
    elif model_name == "minikm":
        ml_model = MiniBatchKMeans(n_clusters=n_clusters, random_state=3, batch_size=1000, max_iter=100)
    else:
        logger.error(f"No model defined for {model_name}.")
        raise NotImplementedError

    ml_model.fit(df[cols_kmers])

    # Model saving
    with open(path_model, 'wb') as f:
        pickle.dump(ml_model, f)
    logger.info(f"{model_name} model saved for k={k} s={w} at {path_model}, now predicting bins for each segment...")

    # ## 3 ##
    predicted = ml_model.predict(df[cols_kmers])
    df["cluster"] = predicted

    df[list(cols_spe) + ["cluster"]].to_pickle(output_pred)
    logger.info(f"Defined {n_clusters} clusters, assignments here: {output_pred} with ML model {model_name}.")
    return


clustering_segments.models = ("minikm", "kmeans")


def pll_copy_segments_to_bin(args):
    """ Function for parallel copying of segments of genomes to a bin, file path and bin number in a dataframe
        Input is only ONE .fna file, which has to be split into segments, but these might be recombined
        if their bin association are consecutive.
    """
    genome_path, df = args
    taxon = df.taxon.iloc[0]
    logger.debug(f"Got the segments clustering: {df.shape} (nb of segments, nb of bins) "
                 f"for the genome {osp.split(genome_path)[1]}")

    # Load the entire genome
    genome = Genome(genome_path, taxon, window_size=main.w, k=main.k)
    genome.load_genome()
    # todo: optimize these dataframes. Looks veeery slow.
    # First get the real segmentation depending on cluster continuity of the segments
    # Aggregate segments with same cluster (consecutive values of cluster), get start, end and description updated
    count = 0
    for i, df_split in df.groupby([(df.cluster != df.cluster.shift()).cumsum()]):
        cluster_id = df_split.cluster.iloc[0]
        description= df_split.description.iloc[0]
        category   = df_split.category.iloc[0]
        name       = df_split.name.iloc[0]
        start      = df_split.start.iloc[0]
        end        = df_split.end.iloc[-1]

        path_bin_segment = osp.join(pll_copy_segments_to_bin.path_db_bins, str(cluster_id), f"{taxon}.fna")

        descr = description.replace(" ", "_")  # To avoid issues with bash. Space and non-breaking space
        descr_splits = descr.split("|", maxsplit=4)
        description_new = "|".join(descr_splits[:3]) + f"s:{start}-e:{end-1}|{descr_splits[4]}|"

        # Need to find the genome/plasmid/ and the right chromosome
        seq = genome.records[category][name]
        logger.log(5, f"Adding combined segment {i}, start={start}, end={end-1}, id={seq.id}, "
                      f"from {(end-start)/main.w} seqs, to bin {cluster_id}, file: {path_bin_segment}")

        segment = SeqRecord(seq.seq[start:end], seq.id, seq.name, description_new, seq.dbxrefs,
                            seq.features, seq.annotations, seq.letter_annotations)
        # Append the combined segment to avoid multiple files for the same taxon
        with open(path_bin_segment, "a") as f:
            SeqIO.write(segment, f, "fasta")
        count += 1
    return


pll_copy_segments_to_bin.path_db_bins = ""


@check_step
def split_genomes_to_bins(path_bins_assignments, path_db_bins, clusters):
    """ Write .fna files from the clustering into n bins """
    logger.info(f"deleting existing sub-folders to avoid duplicates by append to existing files at: {path_db_bins}")
    create_n_folders(path_db_bins, clusters, delete_existing=True)

    # Load bin assignment of each segment
    logger.info(f"loading cluster/bin assignment for each genomes' window "
                f"({osp.getsize(path_bins_assignments)/10**9:.2f} GB): {path_bins_assignments}")
    df = pd.read_pickle(path_bins_assignments)

    # Split it per file to allow parallel processing
    logger.debug(f"Split the DF of segments assignments per fna file ({path_bins_assignments}")
    df_per_fna = df.groupby(["fna_path"])

    # Copy in parallel
    pll_copy_segments_to_bin.path_db_bins = path_db_bins
    add_file_with_parameters(path_db_bins, add_description=f"cluster number = {clusters}")

    logger.info(f"Copy genomes segments to their respective bin into {path_db_bins}")
    Genome.initialize_set_k_mers(main.k)
    try:
        with Pool(main.cores) as pool:  # file copy don't need many cores (main.cores)
            results = list(tqdm(pool.imap(pll_copy_segments_to_bin, df_per_fna), total=len(df_per_fna), dynamic_ncols=True))
    except KeyboardInterrupt:
        logger.warning(f"Stopped by user")
    except:
        logger.warning(f"Multiprocessing failed, launching single core version")
        results = []
        for part in tqdm(df_per_fna, total=len(df_per_fna), dynamic_ncols=True):
            results.append(pll_copy_segments_to_bin(part))

    logger.info(f"{len(results)} genomes have been split into {path_db_bins}")


def classifier_param_checker(l_param):
    """ check kraken2-build --help. Default values to feed in, default is ["kraken2", "35", "31", "7"] """
    assert isinstance(l_param, (list, tuple)), TypeError
    assert len(l_param) > 0, "Empty list"
    assert l_param[0] in [i[0] for i in CLASSIFIERS], f"{l_param[0]} does not correspond to existing supported classifiers"

    params = {"name": l_param[0], }
    buffer = []
    for i in range(2, len(l_param), 2):
        params[l_param[i - 1]] = l_param[i]
        buffer.append(l_param[i - 1] + l_param[i])
    s = "_".join(buffer)
    return params, s


@check_step
def add_library(path_refseq_binned, path_bins_hash, n_clusters, classifier):
    """ launch kraken2-build add-to-library. DELETE EXISTING FOLDER !!
        https://htmlpreview.github.io/?https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html#custom-databases
    """
    create_n_folders(path_bins_hash, n_clusters)
    add_file_with_parameters(path_bins_hash, add_description=f"cluster number = {n_clusters}")

    logger.info(f"{classifier} add_to_library, {n_clusters} clusters, under {path_bins_hash} ")
    for cluster in tqdm(range(n_clusters), dynamic_ncols=True):
        bin_id = f"{cluster}/"

        if "kraken2" in classifier:
            # if library exist in another folder (other classifier parameters, but same binning param), make a link to it !
            existing_lib = glob(f"{osp.dirname(path_bins_hash)}/*/{bin_id}/library")
            path_new_lib = osp.join(path_bins_hash, bin_id, "library")

            # If library has already been done, skip it
            if osp.isdir(path_new_lib):
                logger.debug(f"Library {bin_id} already existing. Delete folder if reinstall needed: {path_new_lib}")
            # If done with other parameters, k25, can reuse it
            elif len(existing_lib) > 0:
                os.symlink(existing_lib[0], path_new_lib)
            else:
                cmd = ["find", osp.join(path_refseq_binned, bin_id), "-name", "'*.fna'", "-print0", "|",
                       "xargs", "-P", f"{main.cores}", "-0", "-I{}", "-n1",
                       "kraken2-build", "--add-to-library", "{}", "--db", osp.join(path_bins_hash, bin_id)]
                bash_process(" ".join(cmd), "Adding genomes to kraken2 library")

        elif "centrifuge" in classifier:
            # Concat all .fna files in a bin into one file.
            path_fnas = osp.join(path_bins_hash, bin_id, "library.fna")
            if osp.isfile(path_fnas):
                logger.info(f"Library file for centrifuge, bin {cluster} exists, skipping step")
                continue
            with open(path_fnas, 'w') as concatenated_fna:
                logger.debug(f"for centrifuge library, concatenated fna files into {path_fnas}")
                for path in tqdm(Path(path_refseq_binned, bin_id).rglob("*.fna"), leave=False):
                    with open(path) as fna:
                        concatenated_fna.write(fna.read())
        else:
            raise NotImplementedError(f"classifier unsupported {classifier}")


@check_step
def build_indexes(path_taxonomy, path_classifier, n_clusters, p):
    """ launch kraken build on each bin
        https://htmlpreview.github.io/?https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html#custom-databases
        Skip skipping by checking if folder exists: **check_step NO FOLDER CHECK** (DON'T REMOVE)
    """
    assert osp.isdir(path_taxonomy), logger.error(f"Path to taxonomy doesn't seem to be a directory: {path_taxonomy}")
    add_file_with_parameters(path_classifier, add_description=f"cluster = {n_clusters} \ntaxonomy = {path_taxonomy}")

    logger.info(f"{p['name']} build its {n_clusters} indexes, will take lots of time. Under: {path_classifier}")
    for cluster in tqdm(range(n_clusters), dynamic_ncols=True):
        bin_id = f"{cluster}/"

        if "kraken2" in p['name']:
            # check if hash has already been done
            path_kraken2 = osp.join(path_classifier, bin_id)
            path_kraken2_hash = osp.join(path_kraken2, "hash.k2d")
            if osp.isfile(path_kraken2_hash) and not any([fname.endswith('.tmp') for fname in os.listdir(path_kraken2)]):
                logger.debug(f"Hash table already created, skipping bin {cluster}")
                continue

            # add link to taxonomy
            taxon_in_cluster = osp.join(path_classifier, bin_id, "taxonomy")
            if osp.islink(taxon_in_cluster):
                logger.debug(f"removing existing link at {taxon_in_cluster}")
                os.unlink(taxon_in_cluster)
            os.symlink(path_taxonomy, taxon_in_cluster)

            # Build
            cmd = ["kraken2-build", "--build", "--threads", f"{main.cores}", "--db", path_kraken2,
                   "--kmer-len", p['k'], "--minimizer-len", p['l'], "--minimizer-spaces", p['s'], ]
            bash_process(cmd, "launching kraken2-build")

        elif "centrifuge" in p['name']:
            path_bin = osp.join(path_classifier, bin_id)
            p_seqtxid = Path(path_classifier).parent.joinpath("kraken2/k35_l31_s7", bin_id, "seqid2taxid.map").as_posix()
            path_lib = osp.join(path_bin, "library.fna")
            path_cf = osp.join(path_bin, "cf_index")

            # if one cf_index.1.cf exists, and there's no more *.sa files, and all *.cf files are not empty...
            if osp.isfile(f"{path_cf}.1.cf") and not list(Path(path_bin).rglob("*.sa")) \
                    and all([f.stat().st_size > 0 for f in Path(path_bin).rglob("*.cf")]):
                logger.info(f"index has already been generated, skipping bin {cluster}")
                continue

            cmd = ["centrifuge-build", "-p", f"{main.cores}",
                   "--conversion-table", p_seqtxid,
                   "--taxonomy-tree", osp.join(path_taxonomy, "nodes.dmp"),
                   "--name-table", osp.join(path_taxonomy, "names.dmp"),
                   path_lib, path_cf, ]
            bash_process(cmd, "launching centrifuge-build. Expect very long run time (in hours)")

    logger.info(f"{p['name']} finished building hash tables. " +
                ("You can clean the intermediate files with: kraken2-build --clean {path_bins_hash}/<bin number>"
                 if "kraken2" in p['name'] else "All files, except the index *.[123].cf, can be removed"))


@check_step
def kraken2_full_add_lib(path_refseq, path_output):
    """ Build the hash table with the same genomes, but without binning, for comparison """
    # todo: adapt for centrifuge as well
    delete_folder_if_exists(path_output)
    create_path(path_output)
    add_file_with_parameters(path_output, add_description=f"no binning database for comparison")

    logger.warning(f"DO NOT INTERRUPT this process, you will have restart from scratches.")
    # Add genomes to
    for folder in os.scandir(path_refseq):
        if not osp.isdir(folder.path):
            continue
        if any([to_omit in folder.name for to_omit in main.omit_folders]):
            logger.info(f"skipping {folder.name}")
            continue
        else:
            cmd = ["find", folder.path, "-name", "'*.fna'", "-print0", "|",
                   "xargs", "-P", f"{main.cores}", "-0", "-I{}", "-n1",
                   "kraken2-build", "--add-to-library", "{}", "--db", path_output]
            bash_process(" ".join(cmd), "adding genomes for kraken2 libraries")


@check_step
def kraken2_full_build_hash(taxonomy, path_output, p):

    # Build hash table
    taxon_link = osp.join(path_output, "taxonomy")
    if osp.islink(taxon_link):
        logger.debug(f"removing existing link at {taxon_link}")
        os.unlink(taxon_link)
    os.symlink(taxonomy, taxon_link)
    cmd = ["kraken2-build", "--build", "--threads", f"{main.cores}", "--db", path_output,
           "--kmer-len", p['k'], "--minimizer-len", p['l'], "--minimizer-spaces", p['s'], ]
    bash_process(cmd, f"Launching CMD to build KRAKEN2 Hash, will take lots of time and memory: ")


def kraken2_clean(path_bins_hash, n_clusters):
    """ Use of kraken2-build --clean option to remove temporary files.
        No cleaning by default because the library is the same for various values of k, l and s
    """
    if n_clusters <= 1:
        logger.info(f"kraken2-build --clean, for all the hashes under {path_bins_hash}")
        cmd = ["kraken2-build", "--clean", "--threads", f"{main.cores}", "--db", path_bins_hash]
        bash_process(cmd, "Launching cleaning with kraken2-build --clean")

    else:
        logger.info(f"kraken2-build --clean, for all the hashes under {path_bins_hash}")
        for cluster in tqdm(range(n_clusters), dynamic_ncols=True):
            bin_id = f"{cluster}/"
            cmd = ["kraken2-build", "--clean", "--threads", f"{main.cores}", "--db", osp.join(path_bins_hash, bin_id)]
            bash_process(cmd, "Launching cleaning with kraken2-build --clean")
    logger.info(f"Cleaning done")


#   **************************************************    MAIN   **************************************************   #
def main(folder_database, folder_output, n_clusters, k, window, cores=cpu_count(), skip_existing="111110",
         early_stop=len(check_step.can_skip)-1, omit_folders=("plant", "vertebrate"),
         path_taxonomy="", full_DB=False, k2_clean=False,
         ml_model=clustering_segments.models[0], classifier_param=CLASSIFIERS[0], verbose_lvl=30, no_cython=False,
         combine_rc=True):
    """ Pre-processing of RefSeq database to split genomes into windows, then count their k-mers
        Second part, load all the k-mer counts into one single Pandas dataframe
        Third train a clustering algorithm on the k-mer frequencies of these genomes' windows
        folder_database : RefSeq root folder
        folder_output   : empty root folder to store kmer counts
    """
    _ = init_logger(__package__)  # initialize the global logger
    logger.info("\n*********************************************************************************************************")
    logger.info("**** Starting script **** \n ")
    try:
        if not no_cython:
            global cyt_ext, cython_is_there
            cyt_ext, cython_is_there = import_cython_mod()
            if cython_is_there:
                logger.info(f"Cython is available, initializing variables")
                cyt_ext.set_verbosity(verbose_lvl)
                cyt_ext.init_variables(k)
        else:
            logger.info(f"Falling back on pure python")
        # Common folder name keeping parameters
        param_k_s = f"k{k}_s{window}"
        o_omitted = "" if len(omit_folders) == 0 else "o" + "-".join(omit_folders)
        folder_intermediate_files = osp.join(folder_output, param_k_s, "kmer_counts")
        # Parameters
        main.folder_database= folder_database
        main.omit_folders   = omit_folders
        main.k              = k
        main.w              = window
        main.cores          = cores
        # Set all columns type
        cols_types = {
            "taxon": int, "category": 'category',
            "start": int, "end": int,
            "name": 'category', "description": 'category', "fna_path": 'category',
        }
        codons = codons_without_rev_comp(main.k) if combine_rc else combinations(main.k)
        for key in codons:
            cols_types[key] = float32
        main.cols_types = cols_types

        check_step.timings    = [perf_counter(), ]  # log time spent
        check_step.step_nb    = 0         # For decorator to know which steps has been
        check_step.early_stop = early_stop
        check_step.can_skip   = skip_existing        # Set the skip variable for the decorator of each step
        # Check classifier/kraken2's parameters
        param, s_param = classifier_param_checker(classifier_param)
        # Check that taxonomy wasn't forgotten
        if '0' in check_step.can_skip[5:] and check_step.early_stop >= 5:
            assert osp.isdir(path_taxonomy), NotADirectoryError

        if full_DB:
            # Run kraken2 on the full RefSeq, without binning, for reference
            path_full_hash = osp.join(folder_output, "no-binning", o_omitted, param['name'], s_param)
            kraken2_full_add_lib(folder_database, path_full_hash)
            kraken2_full_build_hash(path_taxonomy, path_full_hash, param)
            if k2_clean: kraken2_clean(path_full_hash, 1)

        else:
            #    KMER COUNTING
            # get kmer distribution for each window of each genome, parallel folder with same structure
            path_individual_kmer_counts = osp.join(folder_intermediate_files, f"counts.k{k}_s{window}")
            scan_RefSeq_kmer_counts(folder_database, path_individual_kmer_counts)

            # combine all kmer distributions into one single file
            path_stacked_kmer_counts = osp.join(folder_intermediate_files, f"all-counts.k{k}_s{window}_{o_omitted}.csv")
            append_genome_kmer_counts(path_individual_kmer_counts, path_stacked_kmer_counts)

            #    CLUSTERING
            # From kmer distributions, use clustering to set the bins per segment
            string_param = f"{ml_model}_b{n_clusters}_k{main.k}_s{main.w}_{o_omitted}"
            folder_by_model = osp.join(folder_output, param_k_s, string_param)
            path_model = osp.join(folder_by_model, f"model.{string_param}.pkl")
            path_segments_clustering = osp.join(folder_by_model, f"segments-clustered.{string_param}.pd")
            clustering_segments(path_stacked_kmer_counts, path_segments_clustering, path_model, n_clusters, ml_model)

            #    CREATING THE DATABASES
            # create the DB for each bin (copy parts of each .fna genomes into a folder with taxonomy id)
            path_refseq_binned = osp.join(folder_by_model, f"RefSeq_binned")
            split_genomes_to_bins(path_segments_clustering, path_refseq_binned, n_clusters)

            # Run kraken2-build add libray
            path_bins_hash = osp.join(folder_by_model, param['name'], s_param)
            add_library(path_refseq_binned, path_bins_hash, n_clusters, param['name'])

            # Run kraken2-build make hash tables
            build_indexes(path_taxonomy, path_bins_hash, n_clusters, param)

            # Cleaning
            if k2_clean and "kraken2" in param['name']: kraken2_clean(path_bins_hash, n_clusters)

    except KeyboardInterrupt:
        check_step.timings.append(perf_counter())  # log time for the last step that has been interrupted
        logger.error("User interrupted")
        logger.error(traceback.format_exc())
    except Exception as e:
        check_step.timings.append(perf_counter())  # log time for the last step that has been interrupted
        logger.error(f"Trace the log file there: {LOGS}")
        logger.exception(e)

    finally:
        # End
        times = check_step.timings
        for i in range(len(times)-1):
            logger.info(f"timing for STEP {i} - {time_to_hms(times[i], times[i+1])}")
        logger.info(f"Script ended, total time of {time_to_hms(times[0], perf_counter())}. \n")


main.folder_database = ""
main.omit_folders    = ""
main.k               = 0
main.w               = 0
main.cores           = 0
main.cols_types      = {}


def arg_parser():
    # Option to display default values, metavar='' to remove ugly capitalized option's names
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('path_database',    help='Database root folder. Supported format: RefSeq 2019',
                                            type=is_valid_directory)
    parser.add_argument('taxonomy',         help='path to taxonomy (absolute path)',    type=is_valid_directory)
    parser.add_argument('path_plot_me',     help="Data folder for PLoT-ME. For each setting (k, s, b), will store the "
                                                 "k-mer counts, bins with genomes'segments, ML models and hash tables "
                                                 "of the classifiers. Allocate around twice as much space as NCBI RefSeq. ",
                                            type=is_valid_directory)

    parser.add_argument('-k', '--kmer',     help='Size of the kmers (default=%(default)d)',
                                            default=4,          type=int, metavar='')
    parser.add_argument('-w', '--window',   help='Segments/windows size to split genomes into (default=%(default)d)',
                                            default=10000,      type=int, metavar='')
    parser.add_argument('-b', '--bins',     help='Number of bins/clusters to split the DB into (default=%(default)d)',
                                            default=10,         type=int, metavar='')

    parser.add_argument('-t', '--threads',  help='Number of threads (default=%(default)d). Use One (1) if your disk is slow',
                                            default=cpu_count(), type=int,  metavar='')
    parser.add_argument('-e', '--early',    help="Early stop. Index of last step to run. "
                                                 "Use -1 to display all steps and paths (DRY RUN)",
                                            default=len(check_step.can_skip)-1, type=int, metavar='',)
    parser.add_argument('-o', '--omit',     help='Omit some folder/families containing these names. '
                                                 'Write names with spaces (defaults=plant vertebrate)',
                                            default=("plant", "vertebrate"), nargs="+", type=str, metavar='')
    parser.add_argument('-s', '--skip_existing', help="By default, skip already existing files/folders (1:skip, 0:no skipping). "
                                                      "111110 means that steps 0 to 4 will be skipped if file/folders "
                                                      "exist. If the script has been stopped in the middle, use 0 to "
                                                      "redo that step. (default=%(default)s)'",
                                            default=check_step.can_skip, type=str, metavar='')

    parser.add_argument('--clean',          help='Make use of kraken2-build --clean to remove temporary files '
                                                 '(library/added/ and others)',
                                            action='store_true',)
    parser.add_argument('--verbosity',      help='Set the verbosity level (default=%(default)d)',
                                            default=30, type=int,  metavar='')
    parser.add_argument('--no_cython',      help='Disable Cython', action='store_true')
    parser.add_argument('-f', '--full_index', help='Build the index with the full RefSeq database for the chosen '
                                                   'classifier, without binning, omitting the directories set by '
                                                   '--omit. Skips all the other steps/processes '
                                                   '(unused: -e, -n, -s). Used for comparison/benchmarking',
                                            action='store_true')
    parser.add_argument('-c', '--classifier', help="classifier's name and its parameters, space separated. "
                                                   "Ex: '--classifier kraken k 35 l 31 s 7' (default), or "
                                                   "'-c centrifuge'. For unsupported classifiers, you can stop after "
                                                   "step 3 (--early), and build their index based on the .fna files "
                                                   "in PLoT-ME/<param>/'RefSeq_binned'",
                                            default=CLASSIFIERS[0], type=str, nargs="+", metavar='')
    # parser.add_argument('-m', '--ml_model', help='name of the model to use for clustering',
    #                                         choices=clustering_segments.models, type=str, metavar='',
    #                                         default=clustering_segments.models[0])
    args = parser.parse_args()

    logger.info(f"Script {__file__} called with {args}")
    main(folder_database=args.path_database, folder_output=args.path_plot_me, n_clusters=args.bins,
         k=args.kmer, window=args.window, cores=args.threads, skip_existing=args.skip_existing,
         early_stop=args.early, omit_folders=tuple(args.omit), path_taxonomy=args.taxonomy,
         full_DB=args.full_index, classifier_param=args.classifier, k2_clean=args.clean, verbose_lvl=args.verbosity,
         no_cython=args.no_cython)


# python ~/Scripts/Reads_Binning/plot_me/classify.py -t 4 -d bins /hdd1000/Reports/ /ssd1500/Segmentation/3mer_s5000/clustered_by_minikm_3mer_s5000_omitted_plant_vertebrate/ -i /ssd1500/Segmentation/Test-Data/Synthetic_from_Genomes/2019-12-05_100000-WindowReads_20-BacGut/2019-12-05_100000-WindowReads_20-BacGut.fastq /ssd1500/Segmentation/Test-Data/Synthetic_from_Genomes/2019-11-26_100000-SyntReads_20-BacGut/2019-11-26_100000-SyntReads_20-BacGut.fastq /ssd1500/Segmentation/Test-Data/ONT_Silico_Communities/Mock_10000-uniform-bacteria-l1000-q8.fastq /ssd1500/Segmentation/Test-Data/ONT_Silico_Communities/Mock_100000-bacteria-l1000-q10.fastq


if __name__ == '__main__':
    arg_parser()


#
# Deprecated method
def kmer_pkl_path(kmer_folder, fna_path, taxo_ext="gff"):
    """ Legacy method, might not use it anymore in the future
        Return the taxonomy id from a description file of a genome.
        We retrieve the taxo id from the .gff file.
        To avoid re-reading file, taxo id is stored into <genome>.taxon
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
    out_path = osp.join(".", kmer_folder, f"{taxo}__{fna_name}__{formatted_bacteria}.pd")
    return taxo, bacteria_name, fna_name, out_path

