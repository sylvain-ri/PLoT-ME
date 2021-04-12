#! /usr/bin/env python3
# coding: utf-8
# cython: language_level=3, infer_types=True, boundscheck=False
# distutils: language=c++
# additional args: profile=True, wraparound=False, cdivision=True
# defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

"""
python3 setup_cython.py build_ext --inplace

TODO
 * check why codons don't reach the end of the loop / function call
 * make the reverse complement
 * make the binning
 * vectorize the distance calculation
 * save into fastq file


ADD LATER, TO MAKE IT COMPATIBLE FOR PURE PYTHON
import cython
@cython.locals(n=cython.int)
def fib_pure_python(n):
    cython.declare(a=cython.int,
                   b=cython.int,
                   i=cython.int)
"""

from plot_me.tools import init_logger

# Lib for Cython
cimport cython
import numpy as np
cimport numpy as np
from libcpp.unordered_map cimport unordered_map


logger = init_logger("cyt_ext")


# ##########################    CONSTANTS AND VARIABLES FOR FUNCTIONS   ##########################
# Related to DNA
cdef:
    unsigned int k = 4
    list nucleotides = [b"A", b"C", b"G", b"T"]
    float [::1] codon_k4 = np.zeros(256, dtype=np.float32)  # [float 0. for _ in range(256)]
    dict nucl_dico = {'A':0,'C':1,'G':2,'T':3,  'a':0,'c':1,'g':2,'t':3, }
    unordered_map [unsigned char, unsigned int] nucl_val
nucl_val[b"A"] = 0
nucl_val[b"a"] = 0
nucl_val[b"C"] = 1
nucl_val[b"c"] = 1
nucl_val[b"G"] = 2
nucl_val[b"g"] = 2
nucl_val[b"T"] = 3
nucl_val[b"t"] = 3

# Related to combining k-mer counts
cdef:
    list codons_list_all
    list codons_list_kept
    unsigned int dim_combined_codons
    unsigned int [:] codons_kept_indexes
    unsigned int [:] codons_orig_indexes
    dict codons_combined


# ##########################             FUNCTIONS             ##########################

# ##########################             UTILITIES             ##########################

cdef combinations(list combi, unsigned int n, list instances=nucleotides):
    """ Return combinations from the char in instances. Using for finding possible k-mers, for a given n/k """
    if n == 1:
        return combi
    else:
        return [f"{a}{k}" for a in combinations(combi, n - 1) for k in instances]

cdef conversion_table_rc = str.maketrans("ACTG", "TGAC")
cdef reverse_complement_string(str seq):
    """ Return reverse complement string """
    return seq.translate(conversion_table_rc)[::-1]

cdef unsigned int n_dim_rc_combined(unsigned int k):
    """ Return the number of dimensions, for a given k, for the unique k-mer counts (forward - reverse complement) """
    return 2**k + (4**k - 2**k)//2

cdef unsigned int codon_addr(str codon):
    """ Take a codon as char array / str and return the address, given its nucleotides """
    cdef unsigned int length = len(codon)
    cdef unsigned int i
    cdef unsigned int total = 0
    cdef unsigned char codon_char
    for i in range(length):
        codon_char = <unsigned char>codon[i]
        total += nucl_val[codon_char] * 4 ** (length-1 - i)
    return total

cdef combine_counts_with_reverse_complement(float[:] counts):
    """ Combine the forward and reverse complement in the k-mer profile  """
    cdef:
         res = np.empty(dim_combined_codons, dtype=np.float32)
         unsigned int i
    for i in range(dim_combined_codons):
        res[i] = counts[codons_orig_indexes[i]] + counts[codons_kept_indexes[i]]
    return res


 # ###################    INITIALIZATION OF VARIABLES    ########################
cdef init_variables():
    """ Initialize k and indexes for fast processing """
    # Build the mapping to convert fast
    logger.debug("Initializing Indexes for k-mer counting ")

    codons_list_all = combinations(nucleotides, k)
    codons_list_kept = []
    dim_combined_codons = n_dim_rc_combined(k)
    codons_kept_indexes = np.zeros(dim_combined_codons, dtype=np.uint32)
    codons_orig_indexes = np.zeros(dim_combined_codons, dtype=np.uint32)
    codons_combined = {}

    cdef unsigned int counter, index_codon, rc_address
    cdef str rc, cod
    counter = 0
    for index_codon, cod in enumerate(codons_list_all):
        rc = reverse_complement_string(cod)

        if index_codon not in codons_combined.keys():
            codons_list_kept.append(cod)
            codons_kept_indexes[counter] = index_codon
            rc_address = codon_addr(rc)
            codons_combined[rc_address] = index_codon
            codons_orig_indexes[counter] = rc_address
        counter += 1

init_variables()


# ##########################           MAIN  FUNCTIONS         ##########################

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# todo try with char [:] instead of bytes
cdef float [::1] kmer_counter_list(const unsigned char [:] stream):
    """
    Counting k-mers for one line/sequence/read. Return an array counts, alphabetically ordered
    :param stream: one line of a fastq/fast file
    :return: array of length n_dim_rc_combined(k)
    """
    # stream = str(stream)
    # codons = codon_template.copy()
    # if debug >= 2: print(len(stream), flush=True)
    logger.debug("comes to the kmer counter")
    cdef:
        float [::1] codons = codon_k4.copy()
        unsigned int last_failed = 0
        unsigned int recover_addr = 0
        unsigned int fails = 0
        unsigned long counter
        unsigned char lettre
    logger.debug(f"codons[0]{codons[0]}, stream[:10]{stream[:10]}")

    cdef unsigned int addr = nucl_val[stream[0]] + 4*nucl_val[stream[1]] + 16*nucl_val[stream[2]] + 64*nucl_val[stream[3]]
    logger.debug(f"addr for the kmer counter: {addr}")
    codons[addr] += 1
    counter = 4
    lettre = stream[counter]
    logger.debug(f"Starting loop")

    while lettre != b'\0':
        try:
            addr = addr // 4 + nucl_val[lettre] * 64
            codons[addr] += 1
        except:
            if last_failed == 0:
                addr = 1024       # force error for next step
                recover_addr = 0  # prepare the next accurate address
                last_failed += 1

            elif last_failed <= 2:
                addr = 1024
                recover_addr = recover_addr // 4 + nucl_val[lettre] * 64
                last_failed += 1

            else:  # last_failed reached 3
                addr = recover_addr
                last_failed = 0
                fails += 1
        counter += 1
        lettre = stream[counter]

    logger.info(f"stream length:{counter}, fails:{fails}, codons:{codons[0]}, {codons[1]}")
    return codons  #, fails

def pycy_kmer_counter_list(sequence):
    cdef unsigned char[:] chaine = sequence
    return kmer_counter_list(chaine)







