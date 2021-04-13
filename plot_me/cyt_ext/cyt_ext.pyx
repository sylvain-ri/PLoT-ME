#! /usr/bin/env python3
# coding: utf-8
# cython: language_level=3, infer_types=True, boundscheck=True, profile=True, cdivision=True, wraparound=False
# distutils: language=c++
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
import cython
cimport cython
import numpy as np
cimport numpy as np
from libcpp.unordered_map cimport unordered_map
from libc.stdio cimport FILE


logger = init_logger("cyt_ext")


# ##########################    CONSTANTS AND VARIABLES FOR FUNCTIONS   ##########################
# Related to DNA
cdef:
    unsigned int verbosity = 30
    str nucleotides = "ACGT"
    float [::1] template_kmer_counts
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
    cdef unsigned int length = len(codon)  # todo replace len() by class.k
    cdef unsigned int i
    cdef unsigned int total = 0
    cdef unsigned char codon_char
    for i in range(length):
        codon_char = <unsigned char>codon[i]  # todo: needed ?
        total += nucl_val[codon_char] * 4 ** (length-1 - i)
    return total

cdef combine_counts_with_reverse_complement(float[:] counts):
    """ Combine the forward and reverse complement in the k-mer profile  """
    cdef:
         cdef float [:] res = np.empty(dim_combined_codons, dtype=np.float32)
         unsigned int i
    for i in range(dim_combined_codons):
        res[i] = counts[codons_orig_indexes[i]] + counts[codons_kept_indexes[i]]
    return res


 # ###################    INITIALIZATION OF VARIABLES    ########################
cdef _init_variables(unsigned int k):
    """ Initialize k and indexes for fast processing """
    # Build the mapping to convert fast
    if verbosity <= 20: logger.debug("Initializing Indexes for k-mer counting ")

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

def init_variables(unsigned int k=4):
    _init_variables(k)

init_variables()


# ##########################           MAIN  FUNCTIONS         ##########################

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# Bytes Faster than char [:]
cdef float [::1] kmer_counter(bytes stream):
    """
    Counting k-mers for one line/sequence/read. Return an array counts, alphabetically ordered
    :param stream: one line of a fastq/fast file
    :return: array of length n_dim_rc_combined(k)
    """
    # stream = str(stream)
    # codons = codon_template.copy()
    # if debug >= 2: print(len(stream), flush=True)
    if verbosity <= 20: logger.debug("comes to the kmer counter")
    cdef:
        float [::1] kmer_counts = template_kmer_counts.copy()
        unsigned int last_failed = 0
        unsigned int recover_addr = 0
        unsigned int fails = 0
        unsigned long counter
        unsigned char lettre
    if verbosity <= 20: logger.debug(f"codons[0]{kmer_counts[0]}, stream[:10]{stream[:10]}")

    cdef unsigned int addr = nucl_val[stream[0]] + 4*nucl_val[stream[1]] + 16*nucl_val[stream[2]] + 64*nucl_val[stream[3]]
    if verbosity <= 20: logger.debug(f"addr for the kmer counter: {addr}")
    kmer_counts[addr] += 1
    counter = 4
    if verbosity <= 20: logger.debug(f"Starting loop")

    while lettre != b'\0':
        lettre = stream[counter]
        try:
            addr = addr // 4 + nucl_val[lettre] * 64
            kmer_counts[addr] += 1
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

    if verbosity <= 20: logger.info(f"stream length:{counter}, fails:{fails}")
    return kmer_counts  #, fails

def py_kmer_counter(sequence):
    """ Python interface for the Cython k-mer counter """
    return kmer_counter(sequence)


# Related to the file reader. Can be replaced by from libc.stdio cimport fopen, flcose, getline ; +10% time
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
    """
    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef char * line = NULL
    cdef size_t l = 0
    cdef ssize_t read

    while True:
        read = getline(&line, &l, cfile)
        if read == -1: break
        yield line

    fclose(cfile)
    return b""


cdef process_file(str filename, str file_format="fastq"):
    """ Read a file and return the k-mer count """
    cdef:
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        list kmer_counts = []

    for line in read_file(filename):
        if line_nb % modulo == 1:
            kmer_counts.append(kmer_counter(line))
        line_nb+= 1
    return kmer_counts

def py_process_file(filename, file_format="fastq"):
    return process_file(filename, file_format)

cdef process_file2(str filename, str file_format="fastq"):
    """ Read a file and return the k-mer count """
    cdef:
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        float [::1] counts

    for line in read_file(filename):
        if line_nb % modulo == 1:
            counts = kmer_counter(line)
        line_nb+= 1
    return counts

def py_process_file2(filename, file_format="fastq"):
    return process_file(filename, file_format)

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
        size_t l = 0
        ssize_t read
        float [::1] counts
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0

    while True:
        read = getline(&line, &l, cfile)
        if read == -1: break
        if line_nb % modulo == 1:
            counts = kmer_counter(line)
        line_nb += 1

    fclose(cfile)
    return counts

def py_process_file3(filename, file_format="fastq"):
    return process_file3(filename, file_format)

