#!/usr/bin/env python3
"""
#############################################################################
Script to divide a Database of genomes (RefSeq), split them into segments and
cluster them into bins according to their k-mer frequency.
Needs a lot of disk space, and RAM according to the largest genome to process.

For 17GB file of combined kmer counts, combining dataframes took up to 55GB,
loading the file up to 35GB.
KMeans crashed because it reached the 60GB RAM.
Using AWS R4.2XLarge instance

** kmer counts DataFrames are under this format:
taxon	category	start	end	name	description	fna_path	AAAA .... TTTT

** cluster/bin assignments trade the nucleotides columns to a "cluster" column

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

import argparse
import subprocess
from copy import deepcopy
from itertools import islice
from multiprocessing import cpu_count, Pool
import os
import os.path as osp
import pandas as pd
import pickle
import re
from time import time, process_time
import traceback

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA

from tqdm import tqdm

# Import paths and constants for the whole project
from tools import PATHS, ScanFolder, is_valid_directory, init_logger, create_path
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


def create_n_folders(path, n):
    """ Create the sub-folders of bins from 0/ to n/ """
    logger.info(f"creates {n} folder under {path}")
    for i in range(n):
        create_path(osp.join(path, str(i)), with_filename=False)


def scale_df_by_length(df, kmer_cols, k, w):
    """ Divide the kmer counts by the length of the segments, and multiply by the number kmer choices"""
    logger.info(f"Scaling the dataframe, converting to float32")
    ratio = 4**k / (w - k + 1)
    for col in tqdm(kmer_cols):
        df[col] = pd.to_numeric(df[col], downcast='float')
        df[col] *= ratio


# Decorator for all these steps
def check_step(func):
    """ Decorator to print steps and check if results have already been computed
        Need the second argument to be the output file/folder: it will check if the file exists / folder isn't empty
    """
    def wrapper(*args, **kwargs):
        # Check arguments for debugging
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)

        # If step already done, skip it
        to_check = args[1]
        if check_step.can_skip[check_step.step_nb] == "1" and \
                (osp.isfile(to_check) or                           # there's already a file or
                 (osp.isdir(to_check) and os.listdir(to_check))):  # there's a folder, and not empty
            logger.info(f"Step {check_step.step_nb} SKIPPING, function {func.__name__}({signature}, "
                        f"Output has already been generated : {to_check}")
            result = None

        else:
            # Time measurement
            start_time = process_time()
            logger.info(f"Step {check_step.step_nb} START, function {func.__name__}({signature})")
            create_path(to_check, with_filename=True if "." in osp.split(to_check)[1] else False)
            result = func(*args, **kwargs)
            # print time spent
            logger.info(f"Step {check_step.step_nb} END, {process_time() - start_time:.3f}s, function {func.__name__}")

        # Step counter
        check_step.step_nb += 1
        return result
    return wrapper


check_step.step_nb = 0
check_step.can_skip = "11111"


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


# parallel_kmer_counting.k = 4
# parallel_kmer_counting.segments = 10000
parallel_kmer_counting.force_recount = True


@check_step
def scan_RefSeq_kmer_counts(scanning, folder_kmers, stop=30, force_recount=False):
    """ Scan through RefSeq, split genomes into segments, count their k-mer, save in similar structure
        Compatible with 2019 RefSeq format hopefully
    """
    # scanning folder Class set up:
    ScanFolder.set_folder_scan_options(scanning=scanning, target=folder_kmers,
                                       ext_find=(".fastq", ".fq", ".fna"), ext_check=".taxon",
                                       ext_create=f".{main.k}mer_count.pd", skip_folders=main.omit_folders)

    logger.info("scanning through all genomes in refseq to count kmer distributions " + scanning)

    # Set constants to avoid arguments passing
    # parallel_kmer_counting.k = k
    # parallel_kmer_counting.segments = segments
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
def define_cluster_bins(path_kmer_counts, output, path_models, n_clusters):
    """ Given a database of segments of genomes in fastq files, split it in n clusters/bins """
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

    # Model learning
    if df_mem < 5*10**9:
        logger.warning(f"df takes {df_mem/10**9:.2f} GB, choosing KMeans")
        name = "KMeans"
        ml_model = KMeans(n_clusters=n_clusters, n_jobs=main.cores, random_state=3)
    else:
        logger.warning(f"df takes {df_mem/10**9:.2f} GB, choosing Mini Batch KMeans")
        name = "miniKM"
        ml_model = MiniBatchKMeans(n_clusters=n_clusters, random_state=3, batch_size=1000, max_iter=100)

    ml_model.fit(df[cols_kmers])

    # Model saving
    path_model = osp.join(path_models, f"{name}_{k}mer_s{w}.pkl")
    create_path(path_model)
    with open(path_model, 'wb') as f:
        pickle.dump(ml_model, f)
    logger.info(f"{name} model saved for k={k} s={w} at {path_model}")

    # ## 3 ##
    predicted = ml_model.predict(df[cols_kmers])
    df["cluster"] = predicted

    df[list(cols_spe) + ["cluster"]].to_pickle(output)
    logger.info(f"Defined {n_clusters} clusters, assignments here: {output} with ML model {name}.")


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
    """ Write .fna files from the binning for kraken build """
    create_n_folders(path_db_bins, clusters)

    # Load bin assignment of each segment
    df = pd.read_pickle(path_bins_assignments)

    # Split it per file to allow parallel processing
    logger.info(f"Split the DF of segments assignments per fna file ({path_bins_assignments}")
    df_per_fna = []
    for file in tqdm(df.fna_path.unique()):
        df_per_fna.append(df[df.fna_path == file])

    # Copy in parallel
    pll_copy_segments_to_bin.path_db_bins = path_db_bins

    logger.info(f"Copy genomes segments to their respective bin into {path_db_bins}")
    with Pool(main.cores) as pool:
        results = list(tqdm(pool.imap(pll_copy_segments_to_bin, islice(df_per_fna, stop if stop > 0 else None)),
                            total=len(df_per_fna)))

    logger.info(f"got {len(results)} results...")


def pll_kraken2_add_lib(cluster_n):

    cmd = ["do", "find", f"{pll_kraken2_add_lib.path_bins_segments}/{cluster_n}/*/",
           "-name", "'*.fna'", "", "",
           "done", "done", ]
    # do
    # find ~/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/$cluster/*/
    #  -name '*.fa' -print0 | xargs -0 -I{} -n1 -P$cores
    #  kraken2-build --add-to-library {}
    #  --db /home/ubuntu/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/Kraken2_building/$cluster/
    # done
    return subprocess.check_output(cmd)


pll_kraken2_add_lib.path_bins_segments = ""
pll_kraken2_add_lib.path_lib = ""


@check_step
def kraken_build(path_db_bins, path_bins_hash, n_clusters):
    """ launch kraken build on each bin """
    tmp_library = osp.join(path_bins_hash, "kraken2_library")
    create_n_folders(tmp_library, n_clusters)
    # Give values for parallel processing
    pll_kraken2_add_lib.path_bins_segments = path_db_bins
    pll_kraken2_add_lib.path_lib = tmp_library

    # todo: run kraken2-build on each subfolder (add to library and build)
    # Add to library
    with Pool(main.cores) as pool:
        results = list(tqdm(pool.imap(pll_kraken2_add_lib, range(n_clusters)),
                            total=len(n_clusters)))

    create_n_folders(path_bins_hash, n_clusters)

    # Build hashes
    with Pool(main.cores) as pool:
        results = list(tqdm(pool.imap(build_hashes, range(n_clusters)),
                            total=len(n_clusters)))

    return


#   **************************************************    MAIN   **************************************************   #
def main(folder_database, folder_output, n_clusters, k, window, cores,
         skip_existing="11111", force_recount=False, early_stop=-1, omit_folders=("plant", "vertebrate")):
    """ Pre-processing of RefSeq database to split genomes into windows, then count their k-mers
        Second part, load all the k-mer counts into one single Pandas dataframe
        Third train a clustering algorithm on the k-mer frequencies of these genomes' windows
        folder_database : RefSeq root folder
        folder_output   : empty root folder to store kmer counts
    """
    # Common folder name keeping parameters
    parameters = f"{k}mer_s{window}"
    folder_intermediate_files = osp.join(folder_output, "read_binning_tmp")
    # Parameters
    check_step.can_skip = skip_existing
    main.omit_folders   = omit_folders
    main.k              = k
    main.w              = window
    main.cores          = cores

    # get kmer distribution for each window of each genome, parallel folder with same structure
    path_individual_kmer_counts = osp.join(folder_intermediate_files, parameters, f"counts_{k}mer_s{window}")
    scan_RefSeq_kmer_counts(folder_database, path_individual_kmer_counts, stop=early_stop, force_recount=force_recount)

    # combine all kmer distributions into one single file
    path_stacked_kmer_counts = osp.join(folder_intermediate_files, parameters, f"_all_counts.{k}mer_s{window}.pd")
    combine_genome_kmer_counts(path_individual_kmer_counts, path_stacked_kmer_counts)

    # From kmer distributions, use clustering to set the bins per segment
    path_segments_to_bins = osp.join(folder_intermediate_files, parameters, f"_genomes_bins_{k}mer_s{window}.pd")
    path_models           = osp.join(folder_intermediate_files, parameters, "models")
    define_cluster_bins(path_stacked_kmer_counts, path_segments_to_bins, path_models, n_clusters)

    # create the DB for each bin (copy parts of each .fna genomes into a folder with taxonomy id)
    path_DB_bins = osp.join(folder_output, parameters, f"_bins_DB")
    split_genomes_to_bins(path_segments_to_bins, path_DB_bins, n_clusters, stop=early_stop)

    # Run kraken2-build into database folder
    path_bins_hash = osp.join(folder_output, parameters, "bins_kraken2_DB")  # Separate hash tables by classifier
    kraken_build(path_DB_bins, path_bins_hash, n_clusters)


main.omit_folders = ""
main.k            = 0
main.w            = 0
main.cores        = 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path_database', type=is_valid_directory,
                        help='Database root folder. Support format: RefSeq 2019.')
    parser.add_argument('path_output_files', type=is_valid_directory,
                        help="Folder for the k-mer counts, bins with genomes'segments, ML models and final hash tables")
    parser.add_argument('-k', '--kmer', default=4, type=int, help='Size of the kmers. Usual value : 4')
    parser.add_argument('-w', '--window', default=10000, type=int, help='Size of each segments/windows of the genomes')
    parser.add_argument('-n', '--number_bins', default=10, type=int, help='Number of bins to split the DB into')
    parser.add_argument('-c', '--cores', default=cpu_count(), type=int, help='Number of threads')
    parser.add_argument('-x', '--debug', default=-1, type=int, help='For debug purpose')
    parser.add_argument('-f', '--force', help='Force recount kmers', action='store_true')
    parser.add_argument('-o', '--omit', nargs="+", type=str, help='Omit some folder/families',
                        default=("plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other"))
    parser.add_argument('-s', '--skip_existing', type=str, default=check_step.can_skip,
                        help="By default, don't redo files that already exist. "
                             "Write 11011 to force redo the 2rd step, 0-indexed. "
                             "To continue counting kmers, write 01111. To recount all kmers, also add option -f")
    args = parser.parse_args()

    # Set the skip variable for the decorator of each step
    check_step.can_skip = args.skip_existing

    # If force recount of the kmer, disable the skip of the step
    if args.force:
        args.skip_existing = "0" + args.skip_existing[1:]

    logger.warning("**** Starting script ****")
    try:
        main(folder_database=args.path_database, folder_output=args.path_output_files, n_clusters=args.number_bins,
             k=args.kmer, window=args.window, cores=args.cores, skip_existing=args.skip_existing,
             force_recount=args.force, early_stop=args.debug,
             omit_folders=tuple(args.omit))
    except KeyboardInterrupt:
        logger.error("User interrupted")
        logger.error(traceback.format_exc())
    except Exception as e:
        logger.exception(e)





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
