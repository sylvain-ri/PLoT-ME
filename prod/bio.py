#!/usr/bin/env python3
# #############################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-12-11
# Reads Binning Project
#
# #############################################################################
#
#
# common resources for biology related functions and Classes
#
import logging
import numpy as np
import os.path as osp
import pickle
import traceback

# todo: check if this logger works
from Bio import SeqRecord, SeqIO
import ete3.ncbi_taxonomy
from tqdm import tqdm

from prod.util import PATHS, util_logger


# #############################################################################
# Methods for nucleotides manipulations
nucleotides = "ACGT"


def kmers_dic(n, choice=nucleotides):
    return {a:0 for a in combinaisons(choice, n)}


def combinaisons(combi, n, instances=nucleotides):
    if n == 1:
        return combi
    else:
        return [f"{a}{n}" for a in combinaisons(combi, n-1) for n in instances]


def seq_to_window(seq, window_size=4):
    """ Return a sliding window from a string """
    for i in range(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def seq_count_kmer(seq, kmer_count, k=4, ignore_N=True):
    """ Count all kmers, ignore kmers with N or other undecided nucleotides
        seq: string nucleotide input
        kmer_count: dict with all combinations of nucleotides, initialized with zeros (kmer_template["AAAA"] = 0, kmer_count["AAAC"] = 0...)
        the new string hashing behaviour, BiopythonWarning: Using str(seq) to use the new behaviour
    """
    util_logger.info('counting kmers')
    wrong_base = "N"*k
    kmer_count[wrong_base] = 0

    try:
        for kmer in seq_to_window(str(seq), k):
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


# #############################################################################
class CustomRead(SeqRecord.SeqRecord):
    """ Customized Read Sequence. Wrapping SeqIO.Record """
    KMER4 = kmers_dic(4)
    FASTQ_PATH = None
    BASE_PATH = None

    # Load the models to be able to apply them on each read
    LDA = pickle.load(open(PATHS.lda_model, 'rb'))
    KMEANS = pickle.load(open(PATHS.kmeans_model, 'rb'))

    def __init__(self, obj, k=4):
        self.logger = logging.getLogger('classify.CustomRead')
        self.logger.debug('Creating new instance')
        # wrap the object
        self._wrapped_obj = obj
        # Additional attributes
        self.k = k
        self.bin = None
        self._kmer_count = None
        self.lda_feat = []
        self.path_out = None

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self._wrapped_obj, attr)

    @property
    def kmer_count(self, ignore_N=True):
        """ common method """
        if self._kmer_count is None:
            self._kmer_count = seq_count_kmer(self.seq, self.KMER4.copy(), self.k, ignore_N=ignore_N)
        return self._kmer_count

    def lda_reduce(self):
        self.logger.info('reducing dimension of kmer frequency to lda representation')
        self.lda_feat = self.LDA.transform(np.fromiter(
            self.kmer_count.values(), dtype=int).reshape(-1, 256))  # Put into 2D one row

    def find_bin(self):
        self.logger.info('finding bins for each read')
        self.bin = self.KMEANS.predict(self.lda_feat)[0]
        self.description = f"{self.description}, bin_id={self.bin}"
        self.path_out = f"{self.BASE_PATH}.bin-{self.bin}.fastq"
        return self.bin

    def to_fastq(self):
        assert self.FASTQ_PATH is not None, AttributeError("Path of the fastq file must first be defined")
        with open(self.path_out, "a") as f:
            SeqIO.write(self, f, "fasta")

    @classmethod
    def set_fastq_path(cls, path_fastq):
        assert osp.isfile(path_fastq), FileNotFoundError(f"{path_fastq} cannot be found")
        cls.FASTQ_PATH = path_fastq
        cls.BASE_PATH = osp.splitext(path_fastq)[0]

    @classmethod
    def bin_reads(cls, fastq_file):
        """ Bin all reads from provide file """
        counter = 0
        cls.set_fastq_path(fastq_file)

        for record in tqdm(SeqIO.parse(fastq_file, "fasta")):
            custom_read = CustomRead(record)
            custom_read.count_kmer()
            custom_read.lda_reduce()
            custom_read.find_bin()
            custom_read.to_fastq()
            counter += 1
        print(counter)

    @staticmethod
    def split_genome(record=SeqRecord.SeqRecord, window=10000):
        """ Split a genome/plasmid into multiple windows, to count the k-mer
            or to create the .fna files for kraken2-build
        """
        can_be = ["plasmid", "chloroplaste", "scaffold", "contig",
                  "chromosome", "complete genome", "whole genome shotgun sequence", ]

        full_seq = record.seq
        len_genome = len(full_seq)
        segments = []
        for i in range(0, len_genome - window, window):
            segment = SeqRecord.SeqRecord(full_seq[i: min(i + window - 1, len_genome - 1)],
                                          record.id, record.name, record.description, record.dbxrefs, record.features,
                                          record.annotations, record.letter_annotations)
            segments.append(segment)

        if "write_fna":
            raise NotImplementedError
        if "count_kmer":
            raise NotImplementedError














