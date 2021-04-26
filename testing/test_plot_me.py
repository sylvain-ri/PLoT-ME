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
    (2, [ "AA",  "AC",  "AG",  "AT",  "CA",  "CC",  "CG",  "CT",  "GA",  "GC",  "GG",  "GT",  "TA",  "TC",  "TG",  "TT", ]),
    (3, ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
         "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
         "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
         "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]),
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
    (2,  {"AA": "TT", "AC": "GT", "AG": "CT", "AT": "AT", "CA": "TG", "CC": "GG", "CG": "CG", "CT": "AG",
          "GA": "TC", "GC": "GC", "GG": "CC", "GT": "AC", "TA": "TA", "TC": "GA", "TG": "CA", "TT": "AA", },
         # dic with only the kept k-mers
         {"AA": 0., "AC": 0., "AG": 0., "AT": 0., "CA": 0., "CC": 0., "CG": 0.,
          "GA": 0., "GC": 0., "TA": 0., },
         # mapping of array indexes, for the forward strand, FROM where it should be added
         # mapping of array indexes, for the forward strand, to where it should be added
         # AA, AC, AG, AT, CA, CC, CG, GA, GC, TA
         [  0,  1,  2,  3,  4,  5,  6,  8,  9, 12, ],
         # mapping of array indexes, for the reverse complement, to where it should be added
         # TT, GT, CT, AT, TG, GG, CG, TC, GC, TA
         [ 15, 11,  7,  3, 14, 10,  6, 13,  9, 12], ),
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
    assert bio.codons_without_rev_comp(2) == ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA", ]


# #############   KMER COUNTER   #######################
seq_w_expected_counts = [
    {'k':2, 'seq': "AAAACCCCGGGGTTTT",
     'counts': {"AA": 3., "AC": 1., "AG": 0., "AT": 0., "CA": 0., "CC": 3., "CG": 1., "CT": 0.,
                "GA": 0., "GC": 0., "GG": 3., "GT": 1., "TA": 0., "TC": 0., "TG": 0., "TT": 3., }, },
    {'k':3, 'seq': "AAAACCCCGGGGTTTT",
     'counts': {"AAA":2., "AAC":1., "AAG":0., "AAT":0., "ACA":0., "ACC":1., "ACG":0., "ACT":0., "AGA":0., "AGC":0., "AGG":0., "AGT":0., "ATA":0., "ATC":0., "ATG":0., "ATT":0.,
                "CAA":0., "CAC":0., "CAG":0., "CAT":0., "CCA":0., "CCC":2., "CCG":1., "CCT":0., "CGA":0., "CGC":0., "CGG":1., "CGT":0., "CTA":0., "CTC":0., "CTG":0., "CTT":0.,
                "GAA":0., "GAC":0., "GAG":0., "GAT":0., "GCA":0., "GCC":0., "GCG":0., "GCT":0., "GGA":0., "GGC":0., "GGG":2., "GGT":1., "GTA":0., "GTC":0., "GTG":0., "GTT":1.,
                "TAA":0., "TAC":0., "TAG":0., "TAT":0., "TCA":0., "TCC":0., "TCG":0., "TCT":0., "TGA":0., "TGC":0., "TGG":0., "TGT":0., "TTA":0., "TTC":0., "TTG":0., "TTT":2.,
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
    ({"AA":2., "AC": 4., "AG": 6., "AT":1., "CA":2., "CC": 3., "CG":4., "CT":5.,
      "GA":5., "GC": 6., "GG": 7., "GT":8., "TA":9., "TC": 1., "TG":2., "TT":3., },
     {"AA":5., "AC":12., "AG":11., "AT":1., "CA":4., "CC":10., "CG":4.,
      "GA":6., "GC": 6., "TA": 9., },
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
    (3, "AAAAAAACCCCCCAAAGGGGGTTTTTT",
     {"AAA":10., "AAC": 2., "AAG": 1., "AAT": 0., "ACA": 0., "ACC": 2., "ACG": 0., "ACT": 0., "AGA": 0., "AGC": 0., "AGG": 1., "ATA": 0., "ATC": 0., "ATG": 0.,
      "CAA": 1., "CAC": 0., "CAG": 0., "CCA": 1., "CCC": 7., "CCG": 0., "CGA": 0., "CGC": 0., "CTA": 0., "CTC": 0.,
      "GAA": 0., "GAC": 0., "GCA": 0., "GCC": 0., "GGA": 0., "GTA": 0., "TAA": 0., "TCA": 0., }),
    (3, "TTTTTTTTTTTGGGGCCCGAAT",
     {"AAA": 9., "AAC": 0., "AAG": 0., "AAT": 1., "ACA": 0., "ACC": 0., "ACG": 0., "ACT": 0., "AGA": 0., "AGC": 0., "AGG": 0., "ATA": 0., "ATC": 0., "ATG": 0.,
      "CAA": 1., "CAC": 0., "CAG": 0., "CCA": 1., "CCC": 3., "CCG": 1., "CGA": 1., "CGC": 0., "CTA": 0., "CTC": 0.,
      "GAA": 1., "GAC": 0., "GCA": 0., "GCC": 2., "GGA": 0., "GTA": 0., "TAA": 0., "TCA": 0., })
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
    cyt_ext.init_variables(k)
    np.testing.assert_array_equal(counts, cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=True, combine=True))

@pytest.mark.parametrize("k, seq, counts", seq_combined_k_mers_array)
def test_cyt_count_and_combine_seq_arrays(k, seq, counts):
    cyt_ext.init_variables(k)
    np.testing.assert_array_equal(counts, cyt_ext.kmer_counter(str.encode(seq), k=k, dictionary=False, combine=True))


scaling_counts = [
    (2, 11,
     #        AA, AC, AG, AT, CA, CC, CG, GA, GC, TA
     np.array([0,  0,  0,  0,        10    ,  0,  0,  0,  0,  0], dtype=np.float32),
     np.array([0,  0,  0,  0, 4**2/(10-2+1),  0,  0,  0,  0,  0], dtype=np.float32), ),  # 1.777777778
    (2, 30,
     #        AA, AC, AG, AT, CA, CC, CG, GA, GC, TA
     np.array([2,  2,  3,  1,  4,  5,  2,  1,  3,  6], dtype=np.float32),
     np.array([0.068965517, 0.068965517, 0.103448276, 0.034482759, 0.137931034, 0.172413793, 0.068965517, 0.034482759, 0.103448276, 0.206896552], dtype=np.float32), ),
    (2, 102,
     #        AA, AC, AG, AT, CA, CC, CG, GA, GC, TA
     np.array([2,  8,  3,  7, 10,  1,  9, 20, 40,  1], dtype=np.float32),
     np.array([0.316831683, 1.267326733, 0.475247525, 1.425742574, 1.584158416, 0.158415842, 1.425742574, 3.168316832, 6.336633663, 0.158415842], dtype=np.float32), ),
    (3, 102,
     #
     #        AAA,AAC,AAG,        AAT,ACA,ACC,ACG,ACT,AGA,        AGC,AGG,ATA,ATC,ATG,CAA,CAC,CAG,CCA,CCC,        CCG,CGA,CGC,CTA,CTC,GAA,        GAC,GCA,GCC,GGA,        GTA,TAA,TCA
     np.array([ 0,  0,  0,         10,  0,  0,  0,  0,  0,          4,  0,  0,  0,  0,  0,  0,  0,  0,  0,          7,  0,  0,  0,  0,  0,          2,  0,  0,  0,          1,  0,  0], dtype=np.float32),
     np.array([ 0,  0,  0,0.434782609,  0,  0,  0,  0,  0,0.173913043,  0,  0,  0,  0,  0,  0,  0,  0,  0,0.304347826,  0,  0,  0,  0,  0,0.086956522,  0,  0,  0,0.043478261,  0,  0], dtype=np.float32), ),
]

@pytest.mark.parametrize("k, length, counts_unscaled, scaled", scaling_counts)
def test_cyt_scale_counts(k, length, counts_unscaled, scaled):
    cyt_ext.init_variables(k)
    cyt_ext.scale_counts(counts_unscaled, k, length)  # Scaling in place
    np.testing.assert_array_almost_equal(scaled, counts_unscaled)


cluster_for_kmer_profile = [
    (2, 8,
     np.array([1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778,  1.777777778, ], dtype=np.float32),
     np.array([[  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, ],
               [  3,  3,  3,  3,  3,  3,  3,  3,  3,  3, ],
               [ -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, ],
               [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ],
               [3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5, ],
               [  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, ],
               [  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, ],
               [  4,  4,  4,  4,  4,  4,  4,  4,  4,  4, ],
               [  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, ],
               [ -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, ]], dtype=np.float32), ),
    #           AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,ATA,ATC,ATG,CAA,CAC,CAG,CCA,CCC,CCG,CGA,CGC,CTA,CTC,GAA,GAC,GCA,GCC,GGA,GTA,TAA,TCA
    (3, 2,
     np.array( [ -1, -2, -1, -1,  0, -1, -1, -1, -1, -3, -1, -1,1.5, -1, -1, -1, -1,-0.5,-1, -1, -1, -1, -1, -1, -1, -1, -1,0.3, -1,-1.7,-1, -1, ], dtype=np.float32),
     np.array([[  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3, ],
               [  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5, ],
               [ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ],
               [-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3.,-3., ],
               [7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2,7.2, ],
               [  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, ],
               [ 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, ],
               [-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14,-14, ],
               [3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9,3.9, ],
               [ -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, ]], dtype=np.float32), ),
]

@pytest.mark.parametrize("k, cluster, counts, centroids", cluster_for_kmer_profile)
def test_cyt_find_cluster(k, cluster, counts, centroids):
    cyt_ext.init_variables(k)
    assert cluster == cyt_ext.find_cluster(counts, centroids)  # Scaling in place



