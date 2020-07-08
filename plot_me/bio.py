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

from plot_me.tools import init_logger

logger = init_logger("bio")


# #############################################################################
# Methods for nucleotides manipulations
nucleotides = "ACGT"


def kmers_dic(n, choice=nucleotides):
    return {a: 0.0 for a in combinaisons(choice, n)}


def combinaisons(combi, n, instances=nucleotides):
    if n == 1:
        return combi
    else:
        return [f"{a}{n}" for a in combinaisons(combi, n-1) for n in instances]


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











