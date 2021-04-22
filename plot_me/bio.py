#!/usr/bin/env python3
"""
#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
common resources for biology related functions and Classes
"""
import traceback

# todo: check if this logger works
import ete3.ncbi_taxonomy

import logging

logger = logging.getLogger(__name__)


# #############################################################################
# Methods for nucleotides manipulations
nucleotides = "ACGT"
nucl_dico = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
             'a': 0, 'c': 1, 'g': 2, 't': 3, }
conversion_table_rc = str.maketrans("ACTG", "TGAC")


def reverse_complement_string(seq):
    """ Return reverse complement string """
    return seq.translate(conversion_table_rc)[::-1]


def codon_addr(codon):
    """ Take a codon as str and return the address, given its nucleotides """
    length = len(codon)
    total = 0
    for i in range(length):
        codon_char = codon[i]
        total += nucl_dico[codon_char] * 4 ** (length-1 - i)
    return total


def n_dim_rc_combined(k):
    """ Return the number of dimensions, for a given k, for the unique k-mer counts (forward - reverse complement) """
    return 2**k + (4**k - 2**k)//2 if k % 2 == 0 else 4**k // 2


def codons_without_rev_comp(k):
    """ Forward codons only, without those which have been reverse complemented """
    l_codons_all = combinations(k)
    l_forward_only = []
    for cod in l_codons_all:
        rc = reverse_complement_string(cod)
        if rc not in l_forward_only:
            l_forward_only.append(cod)
    return l_forward_only


def table_forward_strand_to_rev_comp(k):
    """ Create a dict with the mapping of FORWARD strand to their REVERSE complements """
    l_codons_all = combinations(k)
    d_codons_orig_target = {}
    for index_codon, cod in enumerate(l_codons_all):
        rc = reverse_complement_string(cod)
        d_codons_orig_target[cod] = rc
    return d_codons_orig_target


def combine_counts_forward_w_rc(d_data, k):
    """ Combine forward and reverse codon into one. """
    combined = {}
    for forward, rc in table_forward_strand_to_rev_comp(k).items():
        if rc not in combined.keys():
            if forward == rc:
                combined[forward] = d_data[forward]
            else:
                combined[forward] = d_data[forward] + d_data[rc]
    return combined


def kmers_dic(n, choice=nucleotides):
    """ From a list, create a dict with these keys and value 0.0 """
    return {a: 0.0 for a in combinations(n, choice)}


def combinations(n=4, combi=nucleotides):
    """ Give all combinations of length n for the items in combi (str of list) """
    if n == 1:
        return combi
    else:
        return [f"{a}{b}" for a in combinations(n - 1, combi) for b in combi]


def seq_to_window(seq, window_size=4):
    """ Return a sliding window from a string """
    for i in range(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def seq_count_kmer(seq, kmer_count=None, k=4, ignore_N=True):
    """ Count all kmers, ignore kmers with N or other undecided nucleotides
        seq: string nucleotide input
        kmer_count: dict with all combinations of nucleotides, initialized with zeros (kmer_template["AAAA"] = 0, kmer_count["AAAC"] = 0...)
        the new string hashing behaviour, BiopythonWarning: Using str(seq) to use the new behaviour
    """
    # todo: numba. to speed the counting by a few folds
    # todo: add reverse complement into the same count
    # https://numba.pydata.org/numba-doc/dev/reference/pysupported.html#typed-dict
    # https://numba.pydata.org/numba-doc/dev/user/performance-tips.html#performance-tips
    # https://numba.pydata.org/numba-doc/latest/user/parallel.html
    if kmer_count is None:
        kmer_count = kmers_dic(k)
    logger.log(5, 'counting kmers')
    wrong_base = "N"*k
    kmer_count[wrong_base] = 0.0

    try:
        for kmer in seq_to_window(str(seq).upper(), k):
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


ncbi = ete3.ncbi_taxonomy.NCBITaxa()


# #############################################################################
# Related to taxonomy
def get_desired_ranks(taxid, desired_ranks, tolist=False):
    """ Get all taxonomy ids of desired ranks
        From stackoverflow
        https://stackoverflow.com/questions/36503042/how-to-get-taxonomic-specific-ids-for-kingdom-phylum-class-order-family-gen
    """
    try:
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        if tolist: return [ranks2lineage.get(rank, 0) for rank in desired_ranks]
        else:      return {f'{rank}_id': ranks2lineage.get(rank, 0) for rank in desired_ranks}
    except:
        print(f"retrieval of the lineage of {taxid} failed")
        if tolist: return [0 for rank in desired_ranks]
        else:      return {f'{rank}_id': 0 for rank in desired_ranks}


def get_list_rank(taxids, desired_rank="species"):
    """ Get the taxonomy id at the species level, from a list of id below species level """
    res = []
    for taxid in taxids:
        name = get_desired_ranks(taxid, [desired_rank], tolist=True)[0]
        res.append(name)
    return res











