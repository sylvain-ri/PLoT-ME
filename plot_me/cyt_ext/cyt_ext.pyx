#! /usr/bin/env python3
# coding: utf-8
# cython: language_level=3, infer_types=True, boundscheck=True, profile=True, wraparound=False
# distutils: language=c
# unused args: wraparound=False, cdivision=True
# defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

"""
python3 setup.py build_ext --inplace

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

from plot_me.tools import init_logger

# Lib for Cython
cimport cython     # For @cython.boundscheck(False)
import numpy as np
cimport numpy as np
from libc.stdio cimport FILE
from cython.parallel import prange

DEF ADDR_ERROR  = 8192  # Value above the possible address of the codon: if k=5, max addr is 4**(5+1)
DEF INFO        = 20
DEF DEBUG       = 10

logger = init_logger("cyt_ext")


# ##########################    CONSTANTS AND VARIABLES FOR FUNCTIONS   ##########################
# Related to DNA
cdef:
    unsigned int verbosity   = 30
    str          nucleotides = "ACGT"
    float [::1]  template_kmer_counts
    # dict nucl_dico = {'A':0,'C':1,'G':2,'T':3,  'a':0,'c':1,'g':2,'t':3, }

cdef unsigned int nucl_val(char c) nogil:
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
    list codons_list_all = []
    list codons_list_kept = []
    unsigned int dim_combined_codons
    unsigned int [:] codons_kept_indexes
    unsigned int [:] codons_orig_indexes
    dict codons_combined = {}


# ##########################             FUNCTIONS             ##########################

# ##########################             UTILITIES             ##########################

cdef combinations(str combi, unsigned int k):
    """ Return combinations from the char in instances. Using for finding possible k-mers, for a given n/k """
    if k == 1:
        return combi
    else:
        return [f"{a}{b}" for a in combinations(combi, k - 1) for b in combi]

cdef conversion_table_rc = str.maketrans("ACTG", "TGAC")
cdef reverse_complement_string(str seq):
    """ Return reverse complement string """
    return seq.translate(conversion_table_rc)[::-1]

cdef unsigned int n_dim_rc_combined(unsigned int k):
    """ Return the number of dimensions, for a given k, for the unique k-mer counts (forward - reverse complement) """
    return 2**k + (4**k - 2**k)//2

cdef unsigned int codon_addr(str codon):
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

#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)
cdef combine_counts_with_reverse_complement(float[:] counts):
    """ Combine the forward and reverse complement in the k-mer profile  """
    cdef:
         float [:] res = np.empty(dim_combined_codons, dtype=np.float32)
         int i
#    for i in prange(dim_combined_codons, nogil=True):
    for i in range(dim_combined_codons):
        res[i] = counts[codons_orig_indexes[i]] + counts[codons_kept_indexes[i]]
    return res


 # ###################    INITIALIZATION OF VARIABLES    ########################
cdef _init_variables(unsigned int k, unsigned int logging_level=30):
    """ Initialize k and indexes for fast processing """
    # Build the mapping to convert fast
    if verbosity <= INFO: logger.debug("Initializing Indexes for k-mer counting ")
    global verbosity
    verbosity = logging_level
    logger.setLevel(logging_level)

    global template_kmer_counts
    template_kmer_counts = np.zeros(4**k, dtype=np.float32)  # [float 0. for _ in range(256)]
    global codons_list_all
    codons_list_all = combinations(nucleotides, k)
    global dim_combined_codons
    dim_combined_codons = n_dim_rc_combined(k)
    global codons_kept_indexes
    codons_kept_indexes = np.zeros(dim_combined_codons, dtype=np.uint32)
    global codons_orig_indexes
    codons_orig_indexes = np.zeros(dim_combined_codons, dtype=np.uint32)

    cdef unsigned int counter, index_codon, rc_address
    cdef str rc, cod
    counter = 0
    for index_codon, cod in enumerate(codons_list_all):
        rc = reverse_complement_string(cod)

        if index_codon not in codons_combined.keys():
            global codons_list_kept
            codons_list_kept.append(cod)
            global codons_kept_indexes
            codons_kept_indexes[counter] = index_codon
            rc_address = codon_addr(rc)
            global codons_combined
            codons_combined[rc_address] = index_codon
            global codons_orig_indexes
            codons_orig_indexes[counter] = rc_address
            counter += 1

def init_variables(unsigned int k=4, unsigned int logging_level=30):
    _init_variables(k, logging_level)

init_variables()


# ##########################           MAIN  FUNCTIONS         ##########################

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# Bytes Faster than const unsigned char [:]
cdef float [::1] kmer_counter(char *stream, unsigned int k_value=4):
    """
    Counting k-mers for one line/sequence/read. Return an array counts, alphabetically ordered
    :param stream: one line of a fastq/fast file
    :param k_value: value of k (k-mer)
    :return: array of length n_dim_rc_combined(k)
    """
    # stream = str(stream)
    # codons = codon_template.copy()
    # if debug >= 2: print(len(stream), flush=True)
    if verbosity <= DEBUG: logger.debug("comes to the kmer counter")

    cdef:
        float [::1] kmer_counts = template_kmer_counts.copy()
        unsigned int stream_len = len(stream)
        unsigned int last_failed = 0
        unsigned int addr = 0
        unsigned int max_addr = 4**k_value - 1
        unsigned int new_addr_mul = 4**(k_value - 1)
        unsigned int recover_addr = 0
        unsigned int fails = 0
        unsigned long long counter = k_value
        char letter = b'N'  # initializing the while loop

    if verbosity <= DEBUG: logger.debug(f"empty codons counts[0]{kmer_counts[0]}, stream[:10]{stream[:10]}")
    if stream_len <= k_value:
        if verbosity <= 30: logger.warning("Sequence was shorter than the k used")
        if verbosity <= INFO: logger.warning(stream)
        return kmer_counts

    addr = nucl_val(stream[0]) + 4*nucl_val(stream[1]) + 16*nucl_val(stream[2]) + 64*nucl_val(stream[3])
    if addr > max_addr:
        addr = ADDR_ERROR       # force error for next step
        recover_addr = 0  # prepare the next accurate address
        last_failed += 1
    else:
        kmer_counts[addr] += 1

    if verbosity <= DEBUG: logger.debug(f"Starting loop")
    for counter in range(4, stream_len):
        letter = stream[counter]
        addr = addr // 4 + nucl_val(letter) * new_addr_mul
        if addr < max_addr:
            kmer_counts[addr] += 1
        else:
            if last_failed == 0:
                addr = ADDR_ERROR       # force error for next step
                recover_addr = 0  # prepare the next accurate address
                last_failed += 1

            elif last_failed <= 2:
                addr = ADDR_ERROR
                recover_addr = recover_addr // 4 + nucl_val(letter) * new_addr_mul
                last_failed += 1

            else:  # last_failed reached 3
                addr = recover_addr
                last_failed = 0
                fails += 1

    if verbosity <= INFO: logger.info(f"stream length:{counter}, fails:{fails}")
    return kmer_counts  #, fails

def py_kmer_counter(sequence, k=4):
    """ Python interface for the Cython k-mer counter """
    return kmer_counter(sequence, k)


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
        yields the line AND its length !
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


cdef process_file(str filename, str file_format="fastq"):
    """ Read a file and return the k-mer count """
    cdef:
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        size_t length
        line
        kmer_counts = []

    for line in read_file(filename):
        if line_nb % modulo == 1:
            kmer_counts.append(kmer_counter(line))
        line_nb+= 1
    return kmer_counts

def py_process_file(filename, file_format="fastq"):
    return process_file(filename, file_format)

#
cdef process_file3(str filename, str file_format="fastq"):
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

    while True:
        read = getline(&line, &length, cfile)
        if read == -1: break
        if line_nb % modulo == 1:
            counts = kmer_counter(line)
            list_counts.append(counts)
        line_nb += 1

    fclose(cfile)
    return list_counts

def py_process_file3(filename, file_format="fastq"):
    return process_file3(filename, file_format)



def python_process(filename, file_format="fastq"):

    cdef long long modulo = 4 if file_format.lower() == "fastq" else 2
    cdef long long line_nb = 0
    cdef char * line
    cdef float [:] counts
    list_counts = []
    with open(filename, "rb") as f:
        for line in f:
            if line_nb % modulo == 1:
                counts = kmer_counter(line)
                list_counts.append(counts)
            line_nb += 1
    return list_counts