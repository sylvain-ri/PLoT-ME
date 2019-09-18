#!/usr/bin/env python3
# ######################################################################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-09-18
#
# ######################################################################################################################
# Methods used for various scripts
#
import os
import re
import traceback
import pandas as pd


def read_fna(file_path):
    """ Read a file with DNA only, returns a string with all 80 char long sequences concatenated"""
    with open(file_path) as f:
        rec = f.readlines()
        return "".join(rec[1:]).replace("\n", "")


# K-MER related methods
nucleotides = "ACGT"


def combinations(combi, n, instances=nucleotides):
    if n == 1:
        return combi
    else:
        return [f"{a}{n}" for a in combinations(combi, n-1) for n in instances]


def kmers_dic(n, choice=nucleotides):
    return {a:0 for a in combinations(choice, n)}


def window(fseq, window_size=4):
    for i in range(len(fseq) - window_size + 1):
        yield fseq[i:i+window_size]


def count_kmers(seq, kmer_template, n, w=100):
    """ Count all kmers, ignore kmers with N or other undecided nucleotides
        Return a list of dict for each window (w=100)
    """
    res = []
    current_split = 0
    next_split = current_split + w
    tmp_counts = kmer_template.copy()
    tmp_counts["start"] = current_split

    try:
        for i, kmer in enumerate(window(seq, n)):
            try:
                tmp_counts[kmer] += 1
            except:
                pass
            # To lower the computational need to split into windows
            if i == next_split:
                res.append(tmp_counts)
                current_split = next_split
                next_split += w
                tmp_counts = kmer_template.copy()
                tmp_counts["start"] = current_split

        return i + 1, res
    except Exception as e:
        print("type error: " + str(e))
        print(traceback.format_exc())
        return i, res


# todo : check this first, quite dirty
path_taxo_names = "/home/ubuntu/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/Kraken2_building/taxonomy/names.dmp"
path_kmer_freq  = "/home/ubuntu/Data/kmer_freq/"
taxo_table = pd.read_csv(path_taxo_names, sep="\t|\t")


def kmer_pkl_path(kmer_folder, fna_path, taxo_ext="gff"):
    """ Return a file name based on the taxonomy id instead of the file name.
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
        else:  return
        taxo_start = description.find(identificator)
        taxo = description[taxo_start + len(identificator):
                           taxo_start + description[taxo_start:].find('\n')]

        assert 1 <= len(taxo) <= 8, f"The taxo id search failed, found an id of length {len(taxo)}, \n" \
            f"for the file: {path_gbk} \n" \
            f"found string : {taxo[:min(50, len(taxo))]} ..."

        with open(path_taxon, "w") as f:
            f.write(taxo)

    query = taxo_table[(taxo_table.taxo == int(taxo)) & (taxo_table.class_name == "scientific name")]
    assert query.shape[0] == 1, f"Found {query.shape[0]} matches for the scientific name of taxo {taxo}." \
                                f" Display the taxo table: \n{taxo_table[taxo_table.taxo == int(taxo)]}"
    bacteria_name = query.name.iat[0]

    formatted_bacteria = re.sub('[^A-Za-z0-9]+', '_', bacteria_name)
    out_path = os.path.join(path_kmer_freq, kmer_folder, f"{taxo}__{fna_name}__{formatted_bacteria}.pd")
    return taxo, bacteria_name, fna_name, out_path


