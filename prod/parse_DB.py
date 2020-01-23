#!/usr/bin/env python3
"""
#############################################################################
Script to divide a Database of genomes (RefSeq), split them into segments and
cluster them into bins according to their k-mer frequency.
Needs a lot of disk space, and RAM according to the largest genome to process.
*** STEPS ***
0 -> Scan the given RefSeq, count kmer frequencies per segment for each genome
1 -> Combine these counts into a single file (RAM intensive)
2 -> Scale the values by the length of segments and combination of kmers,
     and apply a clustering algorithm (KMean, mini batch KMeans)
     to find the cluster association of each segment
3 -> Copy these segments of genomes into bins (DISK intensive)
4 -> kraken2-build --add-to-library
5 -> kraken2-build --build

For 17GB file of combined kmer counts, combining counts took up to 55GB,
loading the file up to 35GB, and KMeans crashed when reaching the 60GB RAM.
Using AWS R4.2XLarge instance with 60GB RAM

** kmer counts DataFrames are under this format:
taxon	category	start	end	name	description	fna_path	AAAA .... TTTT
** cluster/bin assignments trade the nucleotides columns to a "cluster" column

Once the bins created, tmp files (kmer counts) can be removed (read_binning_tmp/),
as well as classifier's tmp files (for kraken2: kraken2-build --clean)

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

import argparse
import shutil
import subprocess
from copy import deepcopy
from itertools import islice
from multiprocessing import cpu_count, Pool
import os
import os.path as osp
import pandas as pd
import pickle
import re
from time import time, perf_counter
import traceback

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA

from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, ScanFolder, is_valid_directory, init_logger, create_path, scale_df_by_length
from bio import kmers_dic, ncbi, seq_count_kmer


logger = init_logger('parse_DB')


class Genome:
    """ Genome from RefSeq. Methods to split it into plasmid/genome and into segments """
    categories = ["plasmid", "chloroplast", "scaffold", "contig",
                  "chromosome", "complete genome", "whole genome shotgun sequence", ]
    kmer_count_zeros = kmers_dic(4)

    def __init__(self, fna_file, taxon, window_size, k=4):
        logger.log(0, "Created genome object")
        self.path_fna    = fna_file
        self.taxon       = taxon
        self.window_size = window_size
        self.k           = k
        # records is a dict of SeqRecord
        self.records = {cat: [] for cat in self.categories}
        # self.splits  = {cat: [] for cat in self.categories}

    def __repr__(self):
        return f"Genome object from {osp.split(self.path_fna)[1]}"

    def load_genome(self):
        """ Loop through all chromosomes/plasmids/genomes/.. and parse them into the object
            Split them into various categories (plasmids, genomes, ...)
        """
        logger.debug(f"loading genome {self.path_fna}")
        for record in SeqIO.parse(self.path_fna, "fasta"):
            for cat in self.categories:
                if cat in record.description:
                    self.records[cat].append(record)
                    break

    def yield_genome_split(self):
        """ Split a genome/plasmid into multiple segments, to count the k-mer
            or to create the .fna files for kraken2-build
        """
        for cat in self.categories:
            for record in self.records[cat]:

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
        # todo: consider combining forward and backward kmers as well as complements.
        #  Single counter for AAAT, TAAA, TTTA and ATTT

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
        df.name        = df.name.astype('category')
        df.fna_path    = df.fna_path.astype('category')
        df.to_pickle(path_kmers)
        logger.debug(f"saved kmer count to {path_kmers}")


def create_n_folders(path, n, delete_existing=False):
    """ Create the sub-folders of bins from 0/ to n/ """
    logger.info(f"creates {n} folder under {path}")
    for i in range(n):
        new_path = osp.join(path, str(i))
        if delete_existing and osp.isdir(new_path):
            shutil.rmtree(new_path)
        create_path(new_path, with_filename=False)


def time_to_h_m_s(start, end, fstring=True):
    assert start < end, ArithmeticError(f"The start time is later than the end time: {start} > {end}")
    delay = int(end - start)
    m, s = divmod(delay, 60)
    h, m = divmod(m, 60)
    if fstring:
        return f"{h:d} hours, {m:02d} minutes, {s:02d} seconds"
    else:
        return h, m, s


def add_file_with_parameters(folder, add_description=""):
    """ Add a file in the folder to remind which parameters were used, which folder were omitted """
    path = osp.join(folder, "parameters_RefSeq_binning.txt")
    with open(path, 'w') as f:
        f.write(f"Files created by {__file__} \n"
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
        if check_step.step_nb > check_step.early_stop:
            logger.info(f"Step {check_step.step_nb} EARLY STOP, no run for \t{func.__name__}({signature})")
            result = None

        elif check_step.can_skip[check_step.step_nb] == "1" and \
                (osp.isfile(to_check)                                 # and there's already a file
                 or (osp.isdir(to_check) and os.listdir(to_check))):  # or there's a folder, not empty
            logger.info(f"Step {check_step.step_nb} SKIPPING, function \t{func.__name__}({signature}, "
                        f"Output has already been generated.")
            result = None

        else:
            # Time measurement
            start_time = perf_counter()
            logger.info(f"Step {check_step.step_nb} START, function \t{func.__name__}({signature})")
            create_path(to_check, with_filename=True if "." in osp.split(to_check)[1] else False)
            result = func(*args, **kwargs)
            # print time spent
            logger.info(f"Step {check_step.step_nb} END, {time_to_h_m_s(start_time, perf_counter())}, "
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
    if osp.isfile(fastq.path_target) and parallel_kmer_counting.force_recount is False:
        logger.debug(f"File already existing, skipping ({fastq.path_target})")
        return
    with open(fastq.path_check) as f:
        taxon = int(f.read())
    genome = Genome(fastq.path_abs, taxon,
                    window_size=main.w, k=main.k)
    genome.load_genome()
    genome.count_kmers_to_df(fastq.path_target)


parallel_kmer_counting.force_recount = True


@check_step
def scan_RefSeq_kmer_counts(scanning, folder_kmers, stop=-1, force_recount=False):
    """ Scan through RefSeq, split genomes into segments, count their k-mer, save in similar structure
        Compatible with 2019 RefSeq format hopefully
    """
    # scanning folder Class set up:
    ScanFolder.set_folder_scan_options(scanning=scanning, target=folder_kmers,
                                       ext_find=(".fastq", ".fq", ".fna"), ext_check=".taxon",
                                       ext_create=f".{main.k}mer_count.pd", skip_folders=main.omit_folders)

    logger.info("scanning through all genomes in refseq to count kmer distributions " + scanning)

    # Set constants to avoid arguments passing
    parallel_kmer_counting.force_recount = force_recount
    # Count in parallel. islice() to take a part of an iterable
    with Pool(main.cores) as pool:
        results = list(tqdm(pool.imap(parallel_kmer_counting, islice(ScanFolder.tqdm_scan(with_tqdm=False),
                                                                     stop if stop>0 else None)),
                            total=ScanFolder.count_root_files()))

    # Traditional sequential run
    # for i, fastq in enumerate(ScanFolder.tqdm_scan()):
    #     if osp.isfile(fastq.path_target) and force_recount is False:
    #         logger.debug(f"File already existing, skipping ({fastq.path_target})")
    #         return
    #     with open(fastq.path_check) as f:
    #         taxon = int(f.read())
    #     genome = Genome(fastq.path_abs, fastq.path_target, taxon, segments=segments, k=k)
    #     genome.load_genome()
    #     genome.count_kmers_to_df()
    #     if i > stop >= 0:
    #         logger.warning("Early stop of the scanning")
    #         break

    logger.info(f"{len(results)} genomes have been scanned and kmer counted.")


@check_step
def combine_genome_kmer_counts(folder_kmers, path_df):
    """ Combine single dataframes into one. Might need high memory """
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
    logger.info(f"Combined file of all kmer counts save at: {path_df}")


@check_step
def clustering_segments(path_kmer_counts, output_pred, path_model, n_clusters, model_name="minikm"):
    """ Given a database of segments of genomes in fastq files, split it in n clusters/bins """
    assert model_name in clustering_segments.models, f"model {model_name} is not implemented"
    logger.info(f"Clustering the genomes' segments into {n_clusters} bins. Loading combined kmer counts...")
    k = main.k
    w = main.w

    df = pd.read_pickle(path_kmer_counts)
    cols_kmers = df.columns[-256:]
    cols_spe = df.columns[:-256]

    # ## 1 ## Scaling by length and kmers
    df_mem = df.memory_usage(deep=False).sum()
    logger.info(f"Model loaded, scaling the values to the length of the segments. "
                f"DataFrame size: {df_mem/10**9:.2f} GB.")

    # todo: save intermediate data
    scale_df_by_length(df, cols_kmers, k, w)

    # ## 2 ## Could add PCA

    # Paths
    create_path(path_model)

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
    logger.info(f"{model_name} model saved for k={k} s={w} at {path_model}")

    # ## 3 ##
    predicted = ml_model.predict(df[cols_kmers])
    df["cluster"] = predicted

    df[list(cols_spe) + ["cluster"]].to_pickle(output_pred)
    logger.info(f"Defined {n_clusters} clusters, assignments here: {output_pred} with ML model {model_name}.")
    return


clustering_segments.models = ("minikm", "kmeans")


def pll_copy_segments_to_bin(df):
    """ Function for parallel copying of segments of genomes to a bin, file path and bin number in a dataframe
        Input is only ONE .fna file, which has to be split into segments, but these might be recombined
        if their bin association are consecutive.
    """
    # todo: call the Genome methods, split into segments, recombine consecutive segments,
    #  write the file with taxo to the appropriate bin
    taxon = df.taxon.iloc[0]
    genome_path = df.fna_path.iloc[0]
    logger.debug(f"Got the segments clustering: {df.shape} (nb of segments, nb of bins) "
                 f"for the genome {osp.split(genome_path)[1]}")

    # Load the entire genome
    genome = Genome(genome_path, taxon, window_size=main.w, k=main.k)
    genome.load_genome()
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

        descr = description.replace(" ", "_")  # To avoid issues with bash
        descr_splits = descr.split("|")
        description_new = "|".join(descr_splits[:3] + [f"s:{start}-e:{end-1}"] + descr_splits[4:])

        # Need to find the genome/plasmid/ and the right chromosome
        for seq in genome.records[category]:
            if seq.name == name:
                logger.debug(f"Adding combined segment {i}, start={start}, end={end-1}, id={seq.id}, "
                             f"from {(end-start)/main.w} seqs, to bin {cluster_id}, file: {path_bin_segment}")

                segment = SeqRecord(seq.seq[start:end], seq.id, seq.name, description_new, seq.dbxrefs,
                                    seq.features, seq.annotations, seq.letter_annotations)
                # Append the combined segment to avoid multiple files for the same taxon
                with open(path_bin_segment, "a") as f:
                    SeqIO.write(segment, f, "fasta")
                count += 1
                break
    return


pll_copy_segments_to_bin.path_db_bins = ""


@check_step
def split_genomes_to_bins(path_bins_assignments, path_db_bins, clusters, stop=-1):
    """ Write .fna files from the clustering into n bins """
    logger.info(f"deleting existing sub-folders to avoid duplicates by append to existing files at: {path_db_bins}")
    create_n_folders(path_db_bins, clusters, delete_existing=True)

    # Load bin assignment of each segment
    df = pd.read_pickle(path_bins_assignments)

    # Split it per file to allow parallel processing
    logger.info(f"Split the DF of segments assignments per fna file ({path_bins_assignments}")
    df_per_fna = []
    for file in tqdm(df.fna_path.unique()):
        df_per_fna.append(df[df.fna_path == file])

    # Copy in parallel
    pll_copy_segments_to_bin.path_db_bins = path_db_bins
    add_file_with_parameters(path_db_bins, add_description=f"cluster number = {clusters}")

    logger.info(f"Copy genomes segments to their respective bin into {path_db_bins}")
    with Pool(min(2, main.cores)) as pool:  # file copy don't need many cores (main.cores)
        results = list(tqdm(pool.imap(pll_copy_segments_to_bin, islice(df_per_fna, stop if stop > 0 else None)),
                            total=len(df_per_fna)))

    logger.info(f"got {len(results)} results...")


@check_step
def kraken2_add_lib(path_refseq_binned, path_bins_hash, n_clusters):
    """ launch kraken2-build add-to-library
        https://htmlpreview.github.io/?https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html#custom-databases
    """
    add_file_with_parameters(path_bins_hash, add_description=f"cluster number = {n_clusters}")
    create_n_folders(path_bins_hash, n_clusters)

    logger.info(f"kraken2 add_to_library, {n_clusters} clusters.... ")
    for cluster in tqdm(range(n_clusters)):
        bin_id = f"{cluster}/"
        cmd = ["find", osp.join(path_refseq_binned, bin_id), "-name", "'*.fna'", "-print0", "|",
               "xargs", "-P", f"{main.cores}", "-0", "-I{}", "-n1",
               "kraken2-build", "--add-to-library", "{}", "--db", osp.join(path_bins_hash, bin_id)]
        res = subprocess.call(" ".join(cmd), shell=True, stderr=subprocess.DEVNULL)
        logger.debug(res)
        # with Pool(min(4, main.cores)) as pool:  # file copy don't need many cores (main.cores)
        #     list_files = os.listdir(osp.join(path_refseq_binned, bin_id)
        #     results = list(tqdm(pool.imap(kraken2_build_lib, list_files), total=len(list_files)))


@check_step
def kraken2_build_hash(path_taxonomy, path_bins_hash, n_clusters):
    """ launch kraken build on each bin
        https://htmlpreview.github.io/?https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html#custom-databases
        Skip skipping by checking if folder exists: **check_step NO FOLDER CHECK** (DON'T REMOVE)
    """
    assert osp.isdir(path_taxonomy), logger.error(f"Path to taxonomy doesn't seem to be a directory: {path_taxonomy}")
    add_file_with_parameters(path_bins_hash, add_description=f"cluster = {n_clusters} \ntaxonomy = {path_taxonomy}")

    logger.info(f"kraken2 build its hash tables, {n_clusters} clusters, will take lots of time.... ")
    for cluster in tqdm(range(n_clusters)):
        bin_id = f"{cluster}/"
        taxon_in_cluster = osp.join(path_bins_hash, bin_id, "taxonomy")
        if osp.islink(taxon_in_cluster):
            logger.debug(f"removing existing link at {taxon_in_cluster}")
            os.unlink(taxon_in_cluster)
        os.symlink(path_taxonomy, taxon_in_cluster)

        cmd = ["kraken2-build", "--build", "--threads", f"{main.cores}", "--db", osp.join(path_bins_hash, bin_id)]
        logger.debug(f"Launching CMD to build KRAKEN2 Hash: " + " ".join(cmd))
        res = subprocess.call(cmd)
        logger.debug(res)

    logger.info(f"Kraken2 finished building hash tables. You can clean the intermediate files with: "
                f"kraken2-build --clean {path_bins_hash}/<bin number>")


#   **************************************************    MAIN   **************************************************   #
def main(folder_database, folder_output, n_clusters, k, window, cores=cpu_count(), skip_existing="11111",
         force_recount=False, early_stop=len(check_step.can_skip)-1, omit_folders=("plant", "vertebrate"),
         path_taxonomy="", ml_model=clustering_segments.models[0]):
    """ Pre-processing of RefSeq database to split genomes into windows, then count their k-mers
        Second part, load all the k-mer counts into one single Pandas dataframe
        Third train a clustering algorithm on the k-mer frequencies of these genomes' windows
        folder_database : RefSeq root folder
        folder_output   : empty root folder to store kmer counts
    """
    print("\n*********************************************************************************************************")
    logger.info("**** Starting script **** \n ")
    logger.info(f"Script {__file__} called with {args}")
    try:
        # Common folder name keeping parameters
        parameters = f"{k}mer_s{window}"
        folder_intermediate_files = osp.join(folder_output, parameters, "kmer_counts")
        # Parameters
        main.folder_database= folder_database
        main.omit_folders   = omit_folders
        main.k              = k
        main.w              = window
        main.cores          = cores
        # If force recount of the kmer, disable the skip of the step
        if force_recount:
            skip_existing = "0" + skip_existing[1:]
        check_step.can_skip = skip_existing        # Set the skip variable for the decorator of each step
        check_step.early_stop = early_stop
        check_step.timings.append(perf_counter())  # log time spent

        #    KMER COUNTING
        # get kmer distribution for each window of each genome, parallel folder with same structure
        path_individual_kmer_counts = osp.join(folder_intermediate_files, f"counts_{k}mer_s{window}")
        scan_RefSeq_kmer_counts(folder_database, path_individual_kmer_counts, force_recount=force_recount)

        # combine all kmer distributions into one single file
        omitted = "" if len(omit_folders) == 0 else "_omitted_" + "_".join(omit_folders)
        path_stacked_kmer_counts = osp.join(folder_intermediate_files, f"_all_counts{omitted}.{k}mer_s{window}.pd")
        combine_genome_kmer_counts(path_individual_kmer_counts, path_stacked_kmer_counts)

        #    CLUSTERING
        # From kmer distributions, use clustering to set the bins per segment
        folder_by_model = osp.join(folder_output, parameters, f"clustered_by_{ml_model}_{main.k}mer_s{main.w}{omitted}")
        path_model = osp.join(folder_by_model, f"model_{ml_model}_{main.k}mer_s{main.w}.pkl")
        path_segments_clustering = osp.join(folder_by_model, f"segments_clustered.{main.k}mer_s{main.w}.pd")
        clustering_segments(path_stacked_kmer_counts, path_segments_clustering, path_model, n_clusters, ml_model)

        #    CREATING THE DATABASES
        # create the DB for each bin (copy parts of each .fna genomes into a folder with taxonomy id)
        path_refseq_binned = osp.join(folder_by_model, f"RefSeq_binned")
        split_genomes_to_bins(path_segments_clustering, path_refseq_binned, n_clusters)

        # Run kraken2-build add libray
        path_bins_hash = osp.join(folder_by_model, "kraken2_hash")  # Separate hash tables by classifier
        kraken2_add_lib(path_refseq_binned, path_bins_hash, n_clusters)

        # Run kraken2-build make hash tables
        kraken2_build_hash(path_taxonomy, path_bins_hash, n_clusters)

        # Run kraken2 on the full RefSeq, without binning, for reference
        # todo: Build the full database from RefSeq

    except KeyboardInterrupt:
        logger.error("User interrupted")
        logger.error(traceback.format_exc())
    except Exception as e:
        logger.exception(e)

    finally:
        # End
        times = check_step.timings
        for i in range(len(times)-1):
            logger.info(f"timing for STEP {i} - {time_to_h_m_s(times[i], times[i+1])}")
        logger.info(f"Script ended, total time of {time_to_h_m_s(times[0], perf_counter())}.")
        print()



main.folder_database = ""
main.omit_folders    = ""
main.k               = 0
main.w               = 0
main.cores           = 0


if __name__ == '__main__':
    # Option to display default values, metavar='' to remove ugly capitalized option's names
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('path_database', type=is_valid_directory,
                        help='Database root folder. Support format: RefSeq 2019.')
    parser.add_argument('path_output_files', type=is_valid_directory,
                        help="Folder for the k-mer counts, bins with genomes'segments, ML models and final hash tables")

    parser.add_argument('-t', '--taxonomy', default="", type=str, help='path to the taxonomy', metavar='')
    parser.add_argument('-m', '--ml_model', choices=clustering_segments.models, type=str, metavar='',
                        help='name of the model to use for clustering', default=clustering_segments.models[0])
    parser.add_argument('-k', '--kmer',   default=4, type=int, help='Size of the kmers', metavar='')
    parser.add_argument('-w', '--window', default=10000, type=int, help='Size of each segments/windows of the genomes', metavar='')
    parser.add_argument('-n', '--number_bins', default=10, type=int, help='Number of bins to split the DB into', metavar='')
    parser.add_argument('-c', '--cores', default=cpu_count(), type=int, help='Number of threads', metavar='')

    parser.add_argument('-e', '--early', default=len(check_step.can_skip)-1, type=int, metavar='',
                        help='Early stop. Index of last step to run. Use -1 to display all steps and paths (DRY RUN)')
    parser.add_argument('-o', '--omit', nargs="+", type=str, help='Omit some folder/families. Write names with spaces',
                        default=("plant", "vertebrate"), metavar='')
    parser.add_argument('-f', '--force', help='Force recount kmers (set skip to 0)', action='store_true')
    parser.add_argument('-s', '--skip_existing', type=str, default=check_step.can_skip,
                        help="By default, skip files/folders that already exist. Write 110000 to skip steps 0 and 1. "
                             "To recount all kmers, and stop after combining the kmer dataframes, write -s 011111, "
                             "add option -f, and option -e 0.", metavar='')
    args = parser.parse_args()

    main(folder_database=args.path_database, folder_output=args.path_output_files, n_clusters=args.number_bins,
         k=args.kmer, window=args.window, cores=args.cores, skip_existing=args.skip_existing,
         force_recount=args.force, early_stop=args.early, omit_folders=tuple(args.omit),
         path_taxonomy=args.taxonomy, ml_model=args.ml_model)





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

