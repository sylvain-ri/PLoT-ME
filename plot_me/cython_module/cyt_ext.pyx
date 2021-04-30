#! /usr/bin/env python3
# coding: utf-8
# cython: language_level=3, infer_types=True, boundscheck=True, profile=False, wraparound=False, cdivision=True
# distutils: language=c, define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# unused args: DISABLED_cdivision=True, DISABLED_define_macros=CYTHON_TRACE_NOGIL=1,
# todo: set to False the cython flags, one by one
# todo: use a class
#       https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html
#       https://cython.readthedocs.io/en/latest/src/tutorial/cdef_classes.html
"""
!!!! MUST run this init before using any methods !!!!
python3 -m build
python3 setup.py build_ext --inplace
python3 -m pip install -e .

TODO
 * count k-mers in fastq
 * make the binning
 * vectorize the distance calculation
 * save into fastq file
"""

import logging


# Lib for Cython
cimport cython     # For @cython.boundscheck(False)
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from cython.parallel import prange
# todo: try with cpython.array : https://cython.readthedocs.io/en/latest/src/tutorial/array.html#array-array
# from cpython cimport array
# import array

# ######################################      FILE PROCESSING      ########################################
# Related to the file reader. Can be replaced by from libc.stdio cimport fopen, fclose, getline ; +10% time
# from https://gist.github.com/pydemo/0b85bd5d1c017f6873422e02aeb9618a
from libc.stdio cimport FILE, fopen, fclose, getline, fwrite, fputs, fputc, fprintf

logger = logging.getLogger(__name__)


DEF ADDR_ERROR  = 8192  # Value above the possible address of the codon: if k=5, max addr is 4**(5+1)
DEF WARNING     = 30
DEF INFO        = 20
DEF DEBUG       = 10
DEF DEBUG_MORE  =  5


# ##########################    CONSTANTS AND VARIABLES FOR FUNCTIONS   ##########################
# Related to DNA
cdef:
    unsigned int verbosity          = 30
    unsigned int     k_val          = 0
    unsigned int dim_combined_codons
    str          nucleotides        = "ACGT"

    # Related to combining k-mer counts.
    d_template_counts_all           = {}
    d_template_counts_combined      = {}
    d_codons_orig_target            = {}
    list         l_codons_all       = []
    list         l_codons_combined  = []
    unsigned int [:] ar_codons_forward_addr
    unsigned int [:] ar_codons_rev_comp_addr
    float [:]        template_kmer_counts
    float [:]        template_distances

# ##########################              GETTERS / SETTERS              ##########################
# Getters in Python only. Setter for Verbosity
def get_verbosity():
    return verbosity
def set_verbosity(v):
    global verbosity
    verbosity = v
def get_k_val():
    return k_val
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

def get_template_kmer_counts():
    return template_kmer_counts


# ##########################             FUNCTIONS             ##########################

# ##########################             UTILITIES             ##########################

# https://stackoverflow.com/a/54081075/4767645
cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8(object unicode)

cdef const char** to_cstring_array(list_str):
    cdef const char** ret = <const char**>malloc(len(list_str) * sizeof(char *))
    for i in range(len(list_str)):
        ret[i] = PyUnicode_AsUTF8(list_str[i])
    return ret


# ##########################               BIO               ##########################

@cython.profile(False)       # Allows profiling the function
cdef inline unsigned int nucl_val(char c) nogil:
    """ Map of letters To Value (for addressing the k-mer array """
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

cdef _combinations(int k, str combi=nucleotides):
    """ Return combinations from the char in instances. Using for finding possible k-mers, for a given n/k """
    if k == 1:
        return combi
    else:
        return [f"{a}{b}" for a in _combinations(k - 1, combi) for b in combi]

def combinations(k, combi=nucleotides):
    return _combinations(k, combi)

cdef conversion_table_rc = str.maketrans("ACGT", "TGCA")
cdef str _reverse_complement_string(str seq):
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
        unsigned int length = len(codon)
        unsigned int i
        unsigned int total = 0
        char codon_char
    for i in range(length):
        codon_char = codon[i]
        total += nucl_val(codon_char) * 4 ** (length-1 - i)
    return total

def codon_addr(codon):
    return _codon_addr(codon)


#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)
cdef float[:] _combine_counts_forward_w_rc(float[:] counts):
    """ Return the combined forward and reverse complement from the k-mer profile  """
    if counts.shape[0] != 4 ** k_val:
        LookupError(f"k (={k_val}) hasn't been initialized properly, combining forward and RC can't be made ")
    cdef:
         float [:] res = np.empty(dim_combined_codons, dtype=np.float32)
         unsigned int i
    if verbosity <= DEBUG_MORE: logger.debug(f"Combining forward and reverse complement. count[0]={counts[0]}")
    if verbosity <= DEBUG_MORE:
        logger.log(DEBUG_MORE, f"ar_codons_forward_addr={ar_codons_forward_addr.base}, ar_codons_rev_comp_addr={ar_codons_rev_comp_addr.base}")

    # todo: use prange ? then `i` must be signed
#    for i in prange(dim_combined_codons, nogil=True):
    for i in range(dim_combined_codons):
        if ar_codons_forward_addr[i] == ar_codons_rev_comp_addr[i]:
            res[i] = counts[ar_codons_forward_addr[i]]
        else:
            res[i] = counts[ar_codons_forward_addr[i]] + counts[ar_codons_rev_comp_addr[i]]
    if verbosity <= DEBUG_MORE: logger.debug(f"Combined forward and reverse complement={res.base}")
    return res


def combine_counts_forward_w_rc(counts):
    """ Python API for cython method. combine forward and reverse complement on an array"""
    return _combine_counts_forward_w_rc(counts).base


cdef void _scale_counts(float[:] counts, unsigned int k, ssize_t length):
    """ Scale the counts by their length, in place """
    # todo: should scale by the actual number of columns (palindromes and reverse complemented k-mers)
    #       then do it in parse_DB tools.scale_df_by_length() as well
    cdef:
        unsigned int i
        float divisor, factor

    divisor = <float>length - <float>k + 1.
    factor  = 4.**k / divisor

    for i in range(0, counts.shape[0]):
        counts[i] = counts[i] * factor

def scale_counts(counts, k, length):
    """ Scale the counts by their length, in place """
    return _scale_counts(counts, k, length)

cdef unsigned int _find_cluster(float[:] counts, const float[:,::1] centers):
    """ Compute the distance to each centroid, given the centers for each centroid, for all dimensions 
        Return the cluster number (unsigned int)
    """
    cdef unsigned int cluster_nb = centers.shape[0]
    cdef float [:] distances = template_distances    # copy template with zeros, faster than initializing each time
    cdef float tmp_distance
    cdef unsigned int cluster_choice = 0
    cdef unsigned int i, j
    if verbosity <= DEBUG_MORE: logger.debug(f"Find_clusters: Centroids of shape={centers.shape[0]},{centers.shape[1]}, counts of shape={counts.shape[0]}")

    # Compute the distance to each centroid
    for i in range(cluster_nb):
        for j in range(centers.shape[1]):
            tmp_distance = counts[j] - centers[i][j]
            distances[i] += tmp_distance * tmp_distance

        # Find the cluster with the minimum distance
        if i != 0 and distances[cluster_choice] > distances[i]:
            cluster_choice = i
    return cluster_choice

def find_cluster(counts, centroid_centers):
    return _find_cluster(counts, centroid_centers)


cdef _copy_read_to_bin(char** outputs, unsigned int cluster, char* line_0, char* line_1, char* line_2, char* line_3):
    """" Write the 4 lines of a read to its designated cluster. """

    cdef FILE* file_write_ptr
    file_write_ptr = fopen(outputs[cluster], "ab")
    if file_write_ptr == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % outputs[2])

    # fprintf(file_write_ptr, "%d", cluster)
    fputs(line_0, file_write_ptr)
    fputs(line_1, file_write_ptr)
    fputs(line_2, file_write_ptr)
    fputs(line_3, file_write_ptr)

    fclose(file_write_ptr)
    return


 # ###################    INITIALIZATION OF VARIABLES    ########################
cdef void _init_variables(unsigned int k):
    """ Initialize k and indexes for fast processing """
    # Build the mapping to convert fast
    if verbosity <= INFO: logger.info("Initializing Indexes for k-mer counting ")

    global dim_combined_codons
    dim_combined_codons = _n_dim_rc_combined(k)
    cdef unsigned int dim_k = 4 ** k

    cdef:
        template_counts_all       = {}
        template_counts_combined  = {}
        codons_orig_target        = {}
        codons_all                = _combinations(k)
        codons_combined           = []
        unsigned int [:] codons_forward_addr = np.zeros(dim_combined_codons, dtype=np.uint32)
        unsigned int [:] codons_rev_comp_addr = np.zeros(dim_combined_codons, dtype=np.uint32)

    cdef unsigned int counter, index_codon, rc_address
    cdef str rc, forward
    counter = 0

    if verbosity <= DEBUG: logger.debug(f"Initializing variables with this list of k-mers: {codons_all}")
    for forward in codons_all:
        rc = _reverse_complement_string(forward)
        rc_address = _codon_addr(rc)
        fw_address = _codon_addr(forward)
        codons_orig_target[forward] = rc
        template_counts_all[forward] = 0

        if verbosity <= DEBUG_MORE:
            logger.log(DEBUG_MORE, f"values: counter={counter:3}, fw_address={fw_address:3}, "
                                   f"rc_addr={_codon_addr(rc):3}, forward={forward:3}, rc={rc:3}, "
                                   f"rc {'IN' if rc in codons_orig_target.keys() else 'NOT in':^6} origin_target")

        if rc not in template_counts_combined.keys():
            codons_combined.append(forward)
            template_counts_combined[forward] = 0

            codons_forward_addr[counter] = fw_address
            codons_rev_comp_addr[counter] = rc_address
            counter += 1
    if verbosity <= DEBUG: logger.debug(f"END of initialization. Final values:")
    if verbosity <= DEBUG: logger.debug(f"d_codons_orig_target={codons_orig_target}")
    if verbosity <= DEBUG: logger.debug(f"l_codons_combined={codons_combined[:min(20, dim_combined_codons)]}")
    if verbosity <= DEBUG: logger.debug(f"d_template_counts_combined={template_counts_combined}")
    if verbosity <= DEBUG: logger.debug(f"codons_forward_addr={codons_forward_addr.base[:min(20, dim_k)]}")
    if verbosity <= DEBUG: logger.debug(f"codons_rev_comp_addr={codons_rev_comp_addr.base[:min(20, dim_k)]}")

    # Setting the values
    global k_val
    k_val = k
    global template_kmer_counts
    template_kmer_counts = np.zeros(4**k, dtype=np.float32)
    global template_distances
    template_distances = np.zeros(dim_combined_codons, dtype=np.float32)
    global l_codons_all
    l_codons_all = codons_all
    global l_codons_combined
    l_codons_combined = codons_combined
    global d_template_counts_all
    d_template_counts_all = template_counts_all
    global d_template_counts_combined
    d_template_counts_combined = template_counts_combined
    global d_codons_orig_target
    d_codons_orig_target = codons_orig_target
    global ar_codons_forward_addr
    ar_codons_forward_addr = codons_forward_addr
    global ar_codons_rev_comp_addr
    ar_codons_rev_comp_addr = codons_rev_comp_addr


def init_variables(k):
    """ MUST run this init before using any methods """
    _init_variables(k)


def init_binner():
    # todo: load the model
    # todo: save weights in an array
    raise NotImplementedError

# ##########################           MAIN  FUNCTIONS         ##########################

# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# @cython.nonecheck(False)    # Check if receives None value
# @cython.cdivision(True)    # No check for div by zero
@cython.profile(False)       # Allows profiling the function
cdef float [:] _kmer_counter(const char *stream, unsigned int k_value=4, ssize_t length=-1):
    """
    Counting k-mers for one line/sequence/read. Return an array counts, alphabetically ordered
    :param length: 
    :param stream: one line of a fastq/fast file
    :param k_value: value of k (k-mer)
    :return: array of length n_dim_rc_combined(k)
    """
    cdef:
        float [:] kmer_counts = template_kmer_counts.copy()
        unsigned int stream_len
        unsigned int addr = 0
        unsigned int next_nucleotide_addr = 0
        unsigned int recover_addr = 0
        unsigned int k_minus_1 = k_value - 1
        unsigned int modulo_addr = 4 ** (k_value - 1)
        unsigned int last_failed = 0
        unsigned int fails = 0
        unsigned long long counter = 0

    if length >= 0:
        stream_len = length
    else:
        stream_len = len(stream)

    if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"codons template.shape={kmer_counts.shape}, for k={k_value}, "
                                                       f"stream length={stream_len} and stream[:20]={stream[:min(20, stream_len)]}")
    if stream_len <= k_value:
        if verbosity <= WARNING: logger.warning(f"Sequence was shorter than the k used {stream}")
        return kmer_counts

    for counter in range(0, stream_len):
        # Value of the nucleotide
        next_nucleotide_addr = nucl_val(stream[counter])
        # todo: might get the linefeed character
        # if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"counter={counter}, letter={stream[counter]}, "
        #    f"nucl_val={next_nucleotide_addr}, last_addr_mod={addr % modulo_addr}, next_addr={(addr % modulo_addr) * 4 + next_nucleotide_addr}")

        if next_nucleotide_addr != ADDR_ERROR:     # If the character was recognized
            addr = (addr % modulo_addr) * 4 + next_nucleotide_addr

            # before we have a full length k-mer, we just add the address, but skip counting the still forming k-mer
            if counter >= k_minus_1:
                kmer_counts[addr] += 1.
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
            elif last_failed < k_minus_1:
                addr = ADDR_ERROR
                recover_addr = (recover_addr % modulo_addr) * 4 + next_nucleotide_addr
                last_failed += 1

            # last iteration, resetting the error variables
            else:  # last_failed reached 3
                addr = recover_addr
                last_failed = 0

    if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"stream length:{counter}, fails:{fails}")
    return kmer_counts  #, fails

def kmer_counter(sequence, k=4, dictionary=True, combine=True, ssize_t length=-1):
    """ Python interface for the Cython k-mer counter
        sequence   : bytes (str.encode(dna_sequence)
        k          : int, k-mer value
        dictionary : flag for the return format (else array)
        combine    : flag for the return format (else all k-mers)
        length     : length of the sequence. If provided, assumes that the sequence is a binary sequence (bytes)
    """
    cdef:
        float [:] kmer_counts
        dict dict_kmer_counts
        unsigned int i
        str key
        char* seq             # Do NOT free it. Get munmap_chunk(): invalid pointer | Fatal Python error: Aborted | Aborted (core dumped)

    if length == -1:
        length = len(sequence)
    if isinstance(sequence, str):
        sequence = str.encode(sequence)
    seq = <char*>sequence

    if dictionary:
        if combine:
            kmer_counts = _combine_counts_forward_w_rc(_kmer_counter(seq, k_value=k, length=length))
            dict_kmer_counts = d_template_counts_combined.copy()
        else:
            kmer_counts = _kmer_counter(seq, k_value=k, length=length)
            dict_kmer_counts = d_template_counts_all.copy()

        if verbosity <= DEBUG_MORE:
            dict_kmers_keys = list(dict_kmer_counts.keys())
            logger.debug(f"MemoryView (len={kmer_counts.shape}, first 10={kmer_counts[0]}) "
                         f"to Dict (template keys={dict_kmers_keys[:5]}")
        for i, key in enumerate(dict_kmer_counts.keys()):
            dict_kmer_counts[key] = kmer_counts[i]
        return dict_kmer_counts

    else:  # raw data, no dictionary, but as numpy arrays
        if combine:
            return _combine_counts_forward_w_rc(_kmer_counter(seq, k_value=k, length=length)).base
        else:
            return np.asarray(_kmer_counter(seq, k_value=k, length=length))


# ######################################      FILE PROCESSING      ########################################
def read_file(filename):
    """ Fast Cython file reader
        return each line as char array and its length
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
    cdef ssize_t length
    while True:
        length = getline(&line, &seed, cfile)
        if length <= -1: break
        if verbosity <= DEBUG_MORE: logger.log(DEBUG_MORE, f"Line read: measured length:{len(line)}, returned length={length}, {line[:10]}")
        yield line, length

    fclose(cfile)
    return (b"", 0)


#
cdef unsigned long long _classify_reads(char* fastq_file, unsigned int k, const float[:,::1] centroid_centers,
                                        const char** outputs, unsigned int modulo=4, unsigned int dev=12):
    """ Fast Cython file reader
        from https://gist.github.com/pydemo/0b85bd5d1c017f6873422e02aeb9618a
        
        For whole file reading:
        in C: https://stackoverflow.com/questions/2029103/correct-way-to-read-a-text-file-into-a-buffer-in-c
        in cython: https://stackoverflow.com/questions/44721835/cython-read-binary-file-into-string-and-return
    """
    _init_variables(k)
    if dim_combined_codons != centroid_centers.shape[1]:
        logger.error(f"The number of dimension, based on the provided k={k_val}, should be {dim_combined_codons}. "
                     f"Nevertheless the centroids have a shape of ({centroid_centers.shape[0]},{centroid_centers.shape[1]}).")

    cdef:
        char * line_0        = NULL
        char * line_1        = NULL
        char * line_2        = NULL
        char * line_3        = NULL
        size_t buffer_size_0 = 0
        size_t buffer_size_1 = 0
        size_t buffer_size_2 = 0
        size_t buffer_size_3 = 0
        ssize_t length_line
        ssize_t length_sequence
        float [:] counts     # = np.empty(4**k, dtype=np.float32)
        float [:] combined   # = np.empty(dim_combined_codons, dtype=np.float32)
        unsigned long long number_of_reads = 0
        unsigned int cluster
    results = []

    cdef FILE* cfile
    cfile = fopen(fastq_file, "rb")

    if verbosity <= INFO: logger.info("File opened, Cython initialized, counting k-mers and splitting fastq file")
    while True:
        # Read line 4 by 4 if fastq, otherwise 2 by 2. Save all
        length_line = getline(&line_0, &buffer_size_0, cfile)
        if length_line < 0: break
        length_sequence = getline(&line_1, &buffer_size_1, cfile)
        if length_line < 0: break
        if modulo == 4:
            length_line = getline(&line_2, &buffer_size_2, cfile)
            if length_line < 0: break
            length_line = getline(&line_3, &buffer_size_3, cfile)
            if length_line < 0: break

        # Count k-mers
        if verbosity <= DEBUG_MORE: logger.debug(f"Lines read, sequence len={length_sequence-1}, counting k-mers now")
        counts = _kmer_counter(line_1, k_value=k, length=length_sequence - 1)     # because getline() keeps \n
        combined = _combine_counts_forward_w_rc(counts)
        if verbosity <= DEBUG_MORE: logger.debug(f"combined counts, shape={combined.shape[0]}, counts[0]={combined[0]}, scaling now")
        # scale
        _scale_counts(combined, k, length_sequence - 1)                           # because getline() keeps \n
        if verbosity <= DEBUG_MORE: logger.debug(f"scaled value={combined[0]}")
        # find cluster
        cluster = _find_cluster(combined, centroid_centers)
        if verbosity <= DEBUG_MORE: logger.debug(f"assigned to cluster={cluster}")
        results.append(cluster)
        # todo: copy read to bin
        _copy_read_to_bin(outputs, cluster, line_0, line_1, line_2, line_3)

        number_of_reads += 1
        if number_of_reads > dev:
            break

    free(line_0)
    free(line_1)
    free(line_2)
    free(line_3)
    fclose(cfile)
    return number_of_reads


def classify_reads(p_fastq, k, centroids, list outputs, file_format="fastq", dev=12):
    """ Interface for Cython's function.
        Read fastq file, count k-mers for each read, find its cluster, copy to read to its bin
    """
    # todo: use ramfs to preload the file ()
    # todo: pre process python variables

    # ouputs to char array
    cdef const char* filename
    cdef const char** p_output_parts = to_cstring_array(outputs)
    cdef unsigned long long number_of_reads

    filename = PyUnicode_AsUTF8(p_fastq)
    cdef float [:,::1] kmeans_centroids = centroids
    cdef unsigned int modulo = 4 if file_format.lower() == "fastq" else 2

    number_of_reads = _classify_reads(filename, k=k, centroid_centers=kmeans_centroids, outputs=p_output_parts,
                                      modulo=modulo, dev=dev)
    # Free memory
    # crashes if trying to free the char** or PyUnicode_AsUTF8()
    # PyMem_Free(filename)
    # for i in range(len(outputs)):
    #     PyMem_Free(p_output_parts[i])
    free(p_output_parts)
    return number_of_reads





# ###################################         UNUSED and DISCONTINUED METHODS        ##################################
cdef _process_file(str filename, unsigned int k=4, file_format="fastq"):
    """ Read a file and return the k-mer count """
    # todo: buggy, doesn't return the whole list of arrays, some are empty after index ~10
    cdef:
        unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
        unsigned long long line_nb = 0
        size_t length
        char * line = NULL
        kmer_counts = []

    for line, length in read_file(filename):
        if line_nb % modulo == 1:
            kmer_counts.append(_kmer_counter(line, k_value=k, length=length))
        line_nb+= 1
    free(line)
    return kmer_counts

def process_file(filename, k, file_format="fastq"):
    return _process_file(filename, k=k, file_format=file_format)


def python_preprocess(filename, k, file_format="fastq"):

    cdef unsigned int modulo = 4 if file_format.lower() == "fastq" else 2
    cdef long long line_nb = 0
    cdef char * line
    cdef float [:] counts
    list_counts = []

    _init_variables(k)

    with open(filename, "rb") as f:
        for line in f:
            if verbosity <= DEBUG: logger.info(f"Line length={len(line)}, seq={line[:20]}")
            if line_nb % modulo == 1:
                counts = _kmer_counter(line, k_value=k, length=len(line))
                if verbosity <= DEBUG: logger.info(f"counts length={len(line)}, counts={counts}")
                list_counts.append(counts)
            line_nb += 1
    free(line)
    return list_counts