#!/usr/bin/env python3
"""
First attempt to add test to project, using pytest
python3 -m pytest -v
python3 -m pytest -v --color=yes --log-level=5 | less -r

Testing both plot_me.bio and plot_me.cython_module.cyt_ext
"""
from plot_me import bio
from plot_me.cython_module import cyt_ext
from plot_me.tools import init_logger

import logging
import numpy as np
import pytest


_ = init_logger(__package__)
logger = logging.getLogger(__name__)
cyt_ext.set_verbosity(5)

# ######################    TESTING COMBINATIONS    ######################
combinations = [
    (2, ["AA", "CA", "GA", "TA", "AC", "CC", "GC", "TC",
         "AG", "CG", "GG", "TG", "AT", "CT", "GT", "TT", ]),
    (3, ["AAA", "CAA", "GAA", "TAA", "ACA", "CCA", "GCA", "TCA", "AGA", "CGA", "GGA", "TGA", "ATA", "CTA", "GTA", "TTA",
         "AAC", "CAC", "GAC", "TAC", "ACC", "CCC", "GCC", "TCC", "AGC", "CGC", "GGC", "TGC", "ATC", "CTC", "GTC", "TTC",
         "AAG", "CAG", "GAG", "TAG", "ACG", "CCG", "GCG", "TCG", "AGG", "CGG", "GGG", "TGG", "ATG", "CTG", "GTG", "TTG",
         "AAT", "CAT", "GAT", "TAT", "ACT", "CCT", "GCT", "TCT", "AGT", "CGT", "GGT", "TGT", "ATT", "CTT", "GTT", "TTT", ]),
]
@pytest.mark.parametrize("k, list_combinations", combinations)
def test_bio_combinations(k, list_combinations):
    assert bio.combinations(k) == list_combinations

@pytest.mark.parametrize("k, list_combinations", combinations)
def test_cyt_combinations(k, list_combinations):
    assert cyt_ext.combinations(k) == list_combinations

# ######################    TESTING REVERSE COMPLEMENT    ######################
codon_reversecomplement = [
    ("AACG", "CGTT"),
    ("TGCA", "TGCA"),
    ("TGC", "GCA"),
    ("TTCCA", "TGGAA")]

@pytest.mark.parametrize("codon, reverse_complement", codon_reversecomplement)
def test_bio_reverse_complement_string(codon, reverse_complement):
    assert bio.reverse_complement_string(codon) == reverse_complement

@pytest.mark.parametrize("codon, reverse_complement", codon_reversecomplement)
def test_cyt_reverse_complement_string(codon, reverse_complement):
    assert cyt_ext.reverse_complement_string(codon) == reverse_complement


# ######################    TESTING K-MER ADDRESS    ######################
kmer_addresses = [
    ("AAAA",    0),
    ("TT",     15),
    ("AACG",    6),
    ("TGCA",  228),
    ("TGC",    57),
    ("TTCCA", 980)]

@pytest.mark.parametrize("k_mer, address", kmer_addresses)
def test_bio_codon_addr(k_mer, address):
    assert bio.codon_addr(k_mer) == address

@pytest.mark.parametrize("k_mer, address", kmer_addresses)
def test_cyt_codon_addr(k_mer, address):
    assert cyt_ext.codon_addr(k_mer) == address


# ######################    TESTING DIMENSION OF AN ARRAY OF COMBINED K-MERS    ######################
expected_array_dimension_for_a_given_k = [
    (2, 10),
    (3, 32),
    (4, 136),
    (5, 512), ]

@pytest.mark.parametrize("k, combined_dim", expected_array_dimension_for_a_given_k)
def test_bio_dim_combined_codons(k, combined_dim):
    assert bio.n_dim_rc_combined(k) == combined_dim

@pytest.mark.parametrize("k, combined_dim", expected_array_dimension_for_a_given_k)
def test_cyt_dim_combined_codons(k, combined_dim):
    assert cyt_ext.n_dim_rc_combined(k) == combined_dim
    cyt_ext.init_variables(k)
    assert cyt_ext.get_dim_combined_codons() == combined_dim


# ######################    TESTING INITIALIZATION OF CYTHON'S VARIABLES    ######################
# Testing if the initialization of values work in cyt_ext
# The arrays are initialized to speed up methods
initialized_data = [
    # k, dict mapping of k-mers to their reverse complement if they Should be combined
    (2,  {"AA": "TT", "CA": "TG", "GA": "TC", "TA": "TA", "AC": "GT", "CC": "GG", "GC": "GC", "TC": "GA",
          "AG": "CT", "CG": "CG", "GG": "CC", "TG": "CA", "AT": "AT", "CT": "AG", "GT": "AC", "TT": "AA", },
         # dic with only the kept k-mers
         {"AA": 0., "CA": 0., "GA": 0., "TA": 0., "AC": 0., "CC": 0., "GC": 0.,
          "AG": 0., "CG": 0., "AT": 0., },
         # mapping of array indexes, for the forward strand, FROM where it should be added
         # AA, CA, GA, TA, AC, CC, GC, AG, CG, AT
         [  0,  4,  8, 12,  1,  5,  9,  2,  6,  3, ],
         # mapping of array indexes, for the reverse complement, FROM where it should be added
         # TT, TG, TC, TA, GT, GG, GC, CT, CG, AT
         [ 15, 14, 13, 12, 11, 10,  9,  7,  6,  3], ),
]

@pytest.mark.parametrize("k, origin_target, template_combined, map_forward_adr, map_rc_adr", initialized_data)
def test_bio_table_rev_comp_to_forward_strand(k, origin_target, template_combined, map_forward_adr, map_rc_adr):
    assert bio.table_forward_strand_to_rev_comp(k) == origin_target

@pytest.mark.parametrize("k, origin_target, template_combined, map_forward_adr, map_rc_adr", initialized_data)
def test_cyt_init_variables(k, origin_target, template_combined, map_forward_adr, map_rc_adr):
    cyt_ext.init_variables(k)
    assert cyt_ext.get_d_codons_orig_target()       == origin_target
    assert cyt_ext.get_d_template_counts_combined() == template_combined
    np.testing.assert_array_equal(cyt_ext.get_ar_codons_forward_addr(),
                                         np.array(map_forward_adr, dtype=np.float32))
    np.testing.assert_array_equal(cyt_ext.get_ar_codons_rev_comp_addr(),
                                         np.array(map_rc_adr, dtype=np.float32))


def test_codons_without_rev_comp():
    assert bio.codons_without_rev_comp(2) == ["AA", "CA", "GA", "TA", "AC", "CC", "GC", "AG", "CG", "AT"]


# #############   KMER COUNTER   #######################
seq_w_expected_counts = [
    {'k':2, 'seq': "AAAACCCCGGGGTTTT",
     'counts': {"AA": 3., "CA": 0., "GA": 0., "TA": 0., "AC": 1., "CC": 3., "GC": 0., "TC": 0.,
                "AG": 0., "CG": 1., "GG": 3., "TG": 0., "AT": 0., "CT": 0., "GT": 1., "TT": 3., }, },
    {'k':3, 'seq': "AAAACCCCGGGGTTTT",
     'counts': {"AAA":2., "CAA":0., "GAA":0., "TAA":0., "ACA":0., "CCA":0., "GCA":0., "TCA":0., "AGA":0., "CGA":0., "GGA":0., "TGA":0., "ATA":0., "CTA":0., "GTA":0., "TTA":0.,
                "AAC":1., "CAC":0., "GAC":0., "TAC":0., "ACC":1., "CCC":2., "GCC":0., "TCC":0., "AGC":0., "CGC":0., "GGC":0., "TGC":0., "ATC":0., "CTC":0., "GTC":0., "TTC":0.,
                "AAG":0., "CAG":0., "GAG":0., "TAG":0., "ACG":0., "CCG":1., "GCG":0., "TCG":0., "AGG":0., "CGG":1., "GGG":2., "TGG":0., "ATG":0., "CTG":0., "GTG":0., "TTG":0.,
                "AAT":0., "CAT":0., "GAT":0., "TAT":0., "ACT":0., "CCT":0., "GCT":0., "TCT":0., "AGT":0., "CGT":0., "GGT":1., "TGT":0., "ATT":0., "CTT":0., "GTT":1., "TTT":2.,
                }, },
]
seq_w_expected_counts_tuples = []
for row in seq_w_expected_counts:
    as_np_array = np.fromiter(row['counts'].values(), dtype=np.float32)
    seq_w_expected_counts_tuples.append((row['k'], row['seq'], row['counts'], as_np_array))

@pytest.mark.parametrize("k, seq, counts, np_counts", seq_w_expected_counts_tuples)
def test_bio_seq_count_kmer(k, seq, counts, np_counts):
    assert bio.seq_count_kmer(seq, k=k) == counts

@pytest.mark.parametrize("k, seq, counts, np_counts", seq_w_expected_counts_tuples)
def test_cyt_seq_count_kmer_return_array(k, seq, counts, np_counts):
    cyt_ext.init_variables(k)
    with pytest.raises(TypeError):
        assert cyt_ext.kmer_counter(4568, k=k, dictionary=True)  == counts
    cyt_counts = cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=False, combine=False)
    np.testing.assert_array_equal(cyt_counts, np_counts)

@pytest.mark.parametrize("k, seq, counts, np_counts", seq_w_expected_counts_tuples)
def test_cyt_seq_count_kmer_return_dict(k, seq, counts, np_counts):
    cyt_ext.init_variables(k)
    assert cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=True, combine=False)  == counts


# ######################    TESTING COMBINATION FOR FORWARD AND REVERSE COMPLEMENT    ######################
# Combining forward and reverse into one array / list
counts_all_to_combined = [
    ({"AA": 2., "CA": 2., "GA": 5., "TA": 9., "AC": 4., "CC": 3., "GC": 6., "TC": 1.,
      "AG": 6., "CG": 4., "GG": 7., "TG": 2., "AT": 1., "CT": 5., "GT": 8., "TT": 3., },
     {"AA": 5., "AC":12., "AG":11., "AT": 1., "CA": 4., "CC":10., "CG": 4.,
      "GA": 6., "GC": 6., "TA": 9., },
     2),
]
np_counts_all_to_combined = []
for item in counts_all_to_combined:
    entry = (np.fromiter(item[0].values(), dtype=np.float32),
             np.fromiter(item[1].values(), dtype=np.float32),
             item[2])
    np_counts_all_to_combined.append(entry)

@pytest.mark.parametrize("dict_counts, dict_combined, k", counts_all_to_combined)
def test_bio_combine_forward_rv(dict_counts, dict_combined, k):
    assert bio.combine_counts_forward_w_rc(dict_counts, k) == dict_combined

@pytest.mark.parametrize("array_counts, array_combined, k", np_counts_all_to_combined)
def test_cyt_combine_counts_with_reverse_complement(array_counts, array_combined, k):
    cyt_ext.init_variables(k)
    np.testing.assert_array_equal(cyt_ext.combine_counts_forward_w_rc(array_counts), array_combined)


# ######################    TESTING WHOLE COUNTING: KMER AND COMBINATION    ######################
# Count k-mers in sequence And combine them
seq_combined_k_mers = [
    (2, "GGGGGGGGGGGGGGGGGGGGG",
     {"AA": 0., "AC": 0., "AG": 0., "AT": 0., "CA": 0., "CC": 20., "CG": 0.,
      "GA": 0., "GC": 0., "TA": 0., }),
    (2, "TTTTTTTTTTTGGGGCCC",
     {"AA":10., "AC": 0., "AG": 0., "AT": 0., "CA": 1., "CC": 5., "CG": 0.,
      "GA": 0., "GC": 1., "TA": 0., }),
    (3, "TTTTTTTTTTTGGGGCCCGAAT",
    {"AAA": 9., "AAC": 0., "AAG": 0., "AAT": 1., "ACA": 0., "ACC": 0., "ACG": 0., "ACT": 0., "AGA": 0., "AGC": 0., "AGG": 0., "AGT": 0., "ATA": 0., "ATC": 0., "ATG": 0., "ATT": 0.,
     "CAA": 1., "CAC": 0., "CAG": 0., "CAT": 0., "CCA": 1., "CCC": 3., "CCG": 1., "CCT": 0., "CGA": 1., "CGC": 0., "CGG": 0., "CGT": 0., "CTA": 0., "CTC": 0., "CTG": 0., "CTT": 0.,
     "GAA": 1., "GAC": 0., "GAG": 0., "GAT": 0., "GCA": 0., "GCC": 2., "GCG": 0., "GCT": 0., "GGA": 0., "GGC": 0., "GGG": 0., "GGT": 0., "GTA": 0., "GTC": 0., "GTG": 0., "GTT": 0.,
     "TAA": 0., "TAC": 0., "TAG": 0., "TAT": 0., "TCA": 0., "TCC": 0., "TCG": 0., "TCT": 0., "TGA": 0., "TGC": 0., "TGG": 0., "TGT": 0., "TTA": 0., "TTC": 0., "TTG": 0., "TTT": 0., })
]
seq_combined_k_mers_array = []
for item in seq_combined_k_mers:
    entry = (item[0], item[1], np.fromiter(item[2].values(), dtype=np.float32), )
    seq_combined_k_mers_array.append(entry)

@pytest.mark.parametrize("k, seq, counts", seq_combined_k_mers)
def test_bio_count_and_combine_seq(k, seq, counts):
    assert bio.combine_counts_forward_w_rc(bio.seq_count_kmer(seq, k=k), k) == counts

@pytest.mark.parametrize("k, seq, counts", seq_combined_k_mers)
def test_cyt_count_and_combine_seq_dicts(k, seq, counts):
    np.testing.assert_array_equal(counts, cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=True, combine=True))

@pytest.mark.parametrize("k, seq, counts", seq_combined_k_mers_array)
def test_cyt_count_and_combine_seq_arrays(k, seq, counts):
    np.testing.assert_array_equal(counts, cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=False, combine=True))








