#! /usr/bin/env python3
# coding: utf-8
# cython: language_level=3, infer_types=True, boundscheck=True, profile=True, wraparound=False
# distutils: language=c, define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# unused args: wraparound=False, cdivision=True

"""
!!!! MUST run this init before using any methods !!!!
python3 setup.py build_ext --inplace
python3 -m pip install -e .

!!!!   CAREFUL   !!!!
While the _kmer_counter counts k-mer reverse way, the combined counts are alphabetically ordered
for _kmer_counter("AACCGGT", k=2)
 return is : 1,  0,  0,  0,  1,  1,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0
 (equiv for AA, CA, GA, TA, AC, CC, GC, TC, AG, CG, GG, TG, AT, CT, GT, TT
 instead of : 1,  1,  0,  0,  0,  1,  1,  0,  0,  0,  1,  1,  0,  0,  0,  0
(for k-mers  AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT

combine_counts_forward_w_rc(counts) will produce
 1,  2,  0,  0,  0,  2,  1,  0,  0,  0
AA, AC, AG, AT, CA, CC, CG, GA, GC, TA



TODO
 * count k-mers in fastq
 * count k-mers in genomes
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

import logging

# Lib for Cython
cimport cython     # For @cython.boundscheck(False)
import numpy as np
cimport numpy as np
from libc.stdio cimport FILE
from cython.parallel import prange

DEF ADDR_ERROR  = 8192  # Value above the possible address of the codon: if k=5, max addr is 4**(5+1)
DEF WARNING     = 30
DEF INFO        = 20
DEF DEBUG       = 10
DEF DEBUG_MORE  =  5

logger = logging.getLogger(__name__)


# ##########################    CONSTANTS AND VARIABLES FOR FUNCTIONS   ##########################
# Related to DNA
cdef:
    unsigned int verbosity   = 30
    unsigned int     k_val   = 0
    str          nucleotides = "ACGT"
    float [::1]  template_kmer_counts
    # dict nucl_dico = {'A':0,'C':1,'G':2,'T':3,  'a':0,'c':1,'g':2,'t':3, }

cdef inline unsigned int nucl_val(char c) nogil:
    """ Map of letters TO value (for addressing the k-mer array """
    if c == b"A" or c == b"a":
        return 0
    elif c == b"C" or c == b"c":
        return 1
    elif c == b"G" or c == b"g":
        return 2
    elif c == b"T" or c == b"t":
        return 3
    else:
        return ADDR_ERROR


# Related to combining k-mer counts.
cdef:
    unsigned int dim_combined_codons
    d_template_counts_all      = {}
    d_template_counts_combined = {}
    d_codons_orig_target       = {}
    list l_codons_all               = []
    list l_codons_combined          = []
    unsigned int [:] ar_codons_forward_addr
    unsigned int [:] ar_codons_rev_comp_addr

# ##########################              GETTERS              ##########################
def get_dim_combined_codons():
    return dim_combined_codons
def get_d_template_counts_all():
    return d_template_counts_all
def get_d_template_counts_combined():
    return d_template_counts_combined
def get_d_codons_orig_target():
    return d_codons_orig_target
def get_l_codons_all():
    return l_codons_all
def get_l_codons_combined():
    return l_codons_combined
def get_ar_codons_forward_addr():
    return ar_codons_forward_addr
def get_ar_codons_rev_comp_addr():
    return ar_codons_rev_comp_addr
def set_verbosity(v):
    global verbosity
    verbosity = v
def get_verbosity():
    return verbosity

# ##########################             FUNCTIONS             ##########################

# ##########################             UTILITIES             ##########################

cdef _combinations(int k, str combi=nucleotides):
    """ Return combinations from the char in instances. Using for finding possible k-mers, for a given n/k """
    if k == 1:
        return combi
    else:
        return [f"{a}{b}" for a in _combinations(k - 1, combi) for b in combi]

def combinations(k, combi=nucleotides):
    return _combinations(k, combi)

cdef conversion_table_rc = str.maketrans("ACGT", "TGCA")
cdef _reverse_complement_string(str seq):
    """ Return reverse complement string """
    return seq.translate(conversion_table_rc)[::-1]

def reverse_complement_string(seq):
    return _reverse_complement_string(seq)

cdef unsigned int _n_dim_rc_combined(unsigned int k):
    """ Return the number of dimensions, for a given k, for the unique k-mer counts (forward - reverse complement) 
        https://stackoverflow.com/questions/40952719/algorithm-to-collapse-forward-and-reverse-complement-of-a-dna-sequence-in-python
    """
    return 2**k + (4**k - 2**k)//2 if k % 2 == 0 else 4**k // 2

def n_dim_rc_combined(k):
    return _n_dim_rc_combined(k)


cdef unsigned int _codon_addr(str codon):
    """ Take a codon as char array / str and return the address, given its nucleotides """
    cdef:
        unsigned int length = len(codon)  # todo replace len() by class.k
        unsigned int i
        unsigned int total = 0
        char codon_char
    for i in range(length):
        codon_char = <char>codon[i]  # todo: needed ?
        total += nucl_val(codon_char) * 4 ** (length-1 - i)
    return total

def codon_addr(codon):
    return _codon_addr(codon)


#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)
cdef float[:] _combine_counts_forward_w_rc(float[:] counts):
    """ Combine the forward and reverse complement in the k-mer profile  """
    assert counts.shape[0] == 4 ** k_val, LookupError(f"k (={k_val}) hasn't been initialized properly, combining forward and RC can't be made ")
    cdef:
         float [:] res = np.empty(dim_combined_codons, dtype=np.float32)
         int i
    if verbosity <= DEBUG: logger.debug(f"Combining forward and reverse complement: {counts.base}")
    if verbosity <= DEBUG_MORE:
        logger.log(DEBUG_MORE, f"ar_codons_forward_addr={ar_codons_forward_addr.base}, ar_codons_rev_comp_addr={ar_codons_rev_comp_addr.base}")

#    for i in prange(dim_combined_codons, nogil=True):
    # todo buggy somewhere, the palindrome have all the same value
    for i in range(dim_combined_codons):
        if ar_codons_forward_addr[i] == ar_codons_rev_comp_addr[i]:
            res[i] = counts[ar_codons_forward_addr[i]]
        else:
            res[i] = counts[ar_codons_forward_addr[i]] + counts[ar_codons_rev_comp_addr[i]]
    if verbosity <= DEBUG: logger.debug(f"Combining forward and reverse complement: {res.base}")
    return res


def combine_counts_forward_w_rc(counts):
    """ Python API for cython method. combine forward and reverse complement on an array"""
    return _combine_counts_forward_w_rc(counts).base


 # ###################    INITIALIZATION OF VARIABLES    ########################
cdef _init_variables(unsigned int k):
    """ Initialize k and indexes for fast processing """
    # Build the mapping to convert fast
    if verbosity <= INFO: logger.info("Initializing Indexes for k-mer counting ")

    global k_val
    k_val = k
    global template_kmer_counts
    template_kmer_counts = np.zeros(4**k, dtype=np.float32)  # [float 0. for _ in range(256)]
    global l_codons_all
    l_codons_all = _combinations(k)
    global l_codons_combined
    l_codons_combined = []
    global dim_combined_codons
    dim_combined_codons = _n_dim_rc_combined(k)
    global ar_codons_forward_addr
    ar_codons_forward_addr = np.zeros(dim_combined_codons, dtype=np.uint32)
    global ar_codons_rev_comp_addr
    ar_codons_rev_comp_addr = np.zeros(dim_combined_codons, dtype=np.uint32)
    global d_template_counts_all
    d_template_counts_all = {}
    global d_template_counts_combined
    d_template_counts_combined = {}
    global d_codons_orig_target
    d_codons_orig_target = {}

    cdef unsigned int counter, index_codon, rc_address
    cdef str rc, forward
    counter = 0

    if verbosity <= DEBUG: logger.debug(f"Initializing variables with this list of k-mers: {l_codons_all}")
    for forward in l_codons_all:
        rc = _reverse_complement_string(forward)
        rc_address = _codon_addr(rc)
        fw_address = _codon_addr(forward)
        global d_codons_orig_target
        d_codons_orig_target[forward] = rc
        global d_template_counts_all
        d_template_counts_all[forward] = 0

        if verbosity <= DEBUG_MORE:
            logger.log(DEBUG_MORE, f"values: counter={counter:3}, fw_address={fw_address:3}, "
                                   f"rc_addr={_codon_addr(rc):3}, forward={forward:3}, rc={rc:3}, "
                                   f"rc {'IN' if rc in d_codons_orig_target.keys() else 'NOT in':^6} origin_target")

        if rc not in d_template_counts_combined.keys():
            global l_codons_combined
            l_codons_combined.append(forward)
            global d_template_counts_combined
            d_template_counts_combined[forward] = 0

            global ar_codons_forward_addr
            ar_codons_forward_addr[counter] = fw_address
            global ar_codons_rev_comp_addr
            ar_codons_rev_comp_addr[counter] = rc_address
            counter += 1
    if verbosity <= DEBUG: logger.debug(f"END of initialization. Final values:")
    if verbosity <= DEBUG: logger.debug(f"d_codons_orig_target={d_codons_orig_target}")
    if verbosity <= DEBUG: logger.debug(f"l_codons_combined={l_codons_combined[:10]} - {l_codons_combined[-10:]}")
    if verbosity <= DEBUG: logger.debug(f"d_template_counts_combined={d_template_counts_combined}")
    if verbosity <= DEBUG: logger.debug(f"ar_codons_forward_addr={ar_codons_forward_addr.base[:10]} - {ar_codons_forward_addr.base[-10:]}")
    if verbosity <= DEBUG: logger.debug(f"ar_codons_rev_comp_addr={ar_codons_rev_comp_addr.base[:10]} - {ar_codons_rev_comp_addr.base[-10:]}")


def init_variables(k):
    """ MUST run this init before using any methods """
    _init_variables(k)


# ##########################           MAIN  FUNCTIONS         ##########################

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# @cython.nonecheck(False)    # Check if receives None value
# @cython.profile(True)       # Allows profiling the function
# @cython.cdivision(False)    # No check for div by zero
# Bytes Faster than const unsigned char [:]
cdef float [::1] _kmer_counter(char *stream, unsigned int k_value=4):
    """
    Counting k-mers for one line/sequence/read. Return an array counts, alphabetically ordered
    :param stream: one line of a fastq/fast file
    :param k_value: value of k (k-mer)
    :return: array of length n_dim_rc_combined(k)
    """
    if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, "comes to the kmer counter")

    cdef:
        float [::1] kmer_counts = template_kmer_counts.copy()
        unsigned int stream_len = len(stream)
        unsigned int addr = 0
        unsigned int next_nucleotide_addr = 0
        unsigned int recover_addr = 0
        unsigned int k_minus_1 = k_value - 1
        unsigned int modulo_addr = 4 ** (k_value - 1)
        unsigned int last_failed = 0
        unsigned int fails = 0
        unsigned long long counter = 0

    if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"codons template.shape={kmer_counts.shape}, for k={k_value}, "
                                                       f"stream length={stream_len} and stream[:20]={stream[:min(20, stream_len)]}")
    if stream_len <= k_value:
        if verbosity <= WARNING: logger.warning(f"Sequence was shorter than the k used {stream}")
        return kmer_counts

    for counter in range(0, stream_len):
        # Value of the nucleotide
        next_nucleotide_addr = nucl_val(stream[counter])
        if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"counter={counter}, letter={stream[counter]}, "
           f"nucl_val={next_nucleotide_addr}, last_addr_mod={addr % modulo_addr}, next_addr={(addr % modulo_addr) * 4 + next_nucleotide_addr}")

        if next_nucleotide_addr != ADDR_ERROR:     # If the character was recognized
            addr = (addr % modulo_addr) * 4 + next_nucleotide_addr

            # before we have a full length k-mer, we just add the address, but skip counting the still forming k-mer
            if counter >= k_minus_1:
                kmer_counts[addr] += 1
            else:
                continue
        else:
            # Failed this time (an unknown character threw an address too big)
            if last_failed == 0:
                fails += 1
                addr = ADDR_ERROR       # force error for next step
                recover_addr = 0        # prepare the next accurate address
                last_failed  = 1        # We just failed, let's count until we can start again

            # wait until the wrong letter goes away
            elif last_failed < k_value-1:
                addr = ADDR_ERROR
                recover_addr = (recover_addr % modulo_addr) * 4 + next_nucleotide_addr
                last_failed += 1

            # last iteration, resetting the error variables
            else:  # last_failed reached 3
                addr = recover_addr
                last_failed = 0

    if verbosity <= DEBUG: logger.debug(f"stream length:{counter}, fails:{fails}")
    return kmer_counts  #, fails

def kmer_counter(sequence, k=4, dictionary=True, combine=True):
    """ Python interface for the Cython k-mer counter """
    cdef:
        float [:] kmer_counts
        dict dict_kmer_counts
        unsigned int i
        str key
        char* seq

    seq = <char*>sequence
    if dictionary:
        if combine:
            kmer_counts = _combine_counts_forward_w_rc(_kmer_counter(seq, k))
            dict_kmer_counts = d_template_counts_combined.copy()
        else:
            kmer_counts = _kmer_counter(seq, k)
            dict_kmer_counts = d_template_counts_all.copy()

        if verbosity <= DEBUG:
            dict_kmers_keys = list(dict_kmer_counts.keys())
            logger.debug(f"MemoryView (len={kmer_counts.shape}, first 10={kmer_counts[0]}) "
                         f"to Dict (template keys={dict_kmers_keys[:5]} - {dict_kmers_keys[-5:]}")
        for i, key in enumerate(dict_kmer_counts.keys()):
            dict_kmer_counts[key] = kmer_counts[i]
        return dict_kmer_counts

    else:  # raw data, no dictionary
        if combine:
            return _combine_counts_forward_w_rc(_kmer_counter(seq, k)).base
        else:
            return _kmer_counter(seq, k).base



# Related to the file reader. Can be replaced by from libc.stdio cimport fopen, fclose, getline ; +10% time
# from https://gist.github.com/pydemo/0b85bd5d1c017f6873422e02aeb9618a
cdef extern from "stdio.h":
    # FILE * fopen ( const char * filename, const char * mode )
    FILE *fopen(const char *, const char *)
    # int fclose ( FILE * stream )
    int fclose(FILE *)
    # ssize_t getline(char **lineptr, size_t *n, FILE *stream);
    ssize_t getline(char **, size_t *, FILE *)

def read_file(filename):
    """ Fast Cython file reader 
        from https://gist.github.com/pydemo/0b85bd5d1c017f6873422e02aeb9618a
        For even faster, could look at : https://github.com/EveryTimeIWill18/Cython_Repo/blob/master/FastFileProcessingWithCython.ipynb
    """
    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef char * line = NULL
    cdef size_t seed = 0
    cdef size_t last_seed = 0
    cdef size_t length = 0
    cdef ssize_t read
    while True:
        read = getline(&line, &seed, cfile)
        if read == -1: break
        # length = seed - last_seed
        # last_seed = seed
        if verbosity <= DEBUG: logger.info(f"Line read: given length={length}, measured length:{len(line)}, {line[:10]}")
        # if length > 1:
        #     if line[length-1] == 10:
        #         print(f"(-1) length is {length} ")
        #         length -= 1
        #     if line[length-2] == 10:
        #         print(f"(-2) length is {length} ")
        #         length -= 2
        yield line

    fclose(cfile)
    return (b"", 0)


cdef _process_file(str filename, str file_format="fastq"):
    """ Read a file and return the k-mer count """
    # todo: buggy, doesn't return the whole list of arrays, some are empty after index ~10
    cdef:
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        size_t length
        char * line = NULL
        kmer_counts = []

    for line in read_file(filename):
        if line_nb % modulo == 1:
            kmer_counts.append(_kmer_counter(line))
        line_nb+= 1
    return kmer_counts

def process_file(filename, file_format="fastq"):
    return _process_file(filename, file_format)

#
cdef _classify_reads(str filename, unsigned int window, str file_format="fastq"):
    """ Fast Cython file reader
        from https://gist.github.com/pydemo/0b85bd5d1c017f6873422e02aeb9618a
    """
    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef:
        char * line = NULL
        size_t length = 0
        ssize_t read
        float [::1] counts
        list list_counts = []
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        # todo: add variables to hold the whole read (id, seq, +, quality)

    while True:
        read = getline(&line, &length, cfile)
        if read == -1: break
        if line_nb % modulo == 1:
            # Count all kmers
            counts = _kmer_counter(line)
            # todo: fold forward and reverse strand

            # todo: find cluster

            # todo: copy read to bin

        line_nb += 1

    fclose(cfile)
    return list_counts

def classify_reads(filename, file_format="fastq"):
    return _classify_reads(filename, file_format)



def python_preprocess(filename, file_format="fastq"):

    cdef long long modulo = 4 if file_format.lower() == "fastq" else 2
    cdef long long line_nb = 0
    cdef char * line
    cdef float [:] counts
    list_counts = []
    with open(filename, "rb") as f:
        for line in f:
            if line_nb % modulo == 1:
                counts = _kmer_counter(line)
                list_counts.append(counts)
            line_nb += 1
    return list_counts