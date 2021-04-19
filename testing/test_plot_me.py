#!/usr/bin/env python3
"""
First attempt to add test to project, using pytest
python3 -m pytest -v

Testing both plot_me.bio and plot_me.cython_module.cyt_ext
"""
import pytest
from plot_me import bio
from plot_me.cython_module import cyt_ext
import numpy as np


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


codon_addresses = [
    ("AAAA",    0),
    ("TT",     15),
    ("AACG",    6),
    ("TGCA",  228),
    ("TGC",    57),
    ("TTCCA", 980)]

@pytest.mark.parametrize("codon, address", codon_addresses)
def test_bio_codon_addr(codon, address):
    assert bio.codon_addr(codon) == address

@pytest.mark.parametrize("codon, address", codon_addresses)
def test_cyt_codon_addr(codon, address):
    assert cyt_ext.codon_addr(codon) == address


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


# Combining forward and reverse into one array / list
counts_all_to_combined = [
    ({"AA":2, "AC": 4, "AG": 6, "AT":1, "CA":2, "CC": 3, "CG":4, "CT":5,
      "GA":5, "GC": 6, "GG": 7, "GT":8, "TA":9, "TC": 1, "TG":2, "TT":3, },
     {"AA":5, "AC":12, "AG":11, "AT":1, "CA":4, "CC":10, "CG":4,
      "GA":6, "GC": 6, "TA": 9, },
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

@pytest.mark.parametrize("dict_counts, dict_combined, k", np_counts_all_to_combined)
def test_cyt_combine_counts_with_reverse_complement(dict_counts, dict_combined, k):
    cyt_ext.init_variables(k)
    res = [np.float32(i) for i in cyt_ext.combine_counts_forward_w_rc(dict_counts)]
    print(res)
    np.testing.assert_array_equal(cyt_ext.combine_counts_forward_w_rc(dict_counts), dict_combined)


# Testing if the initialization of values work in cyt_ext
# The arrays are initialized to speed up methods
initialized_arrays = [
    # k, dict mapping of k-mers to their reverse complement if they Should be combined
    (2,  {"AA": "AA", "AC": "AC", "AG": "AG", "AT": "AT", "CA": "CA", "CC": "CC", "CG": "CG", "CT": "AG",
          "GA": "GA", "GC": "GC", "GG": "CC", "GT": "AC", "TA": "TA", "TC": "GA", "TG": "CA", "TT": "AA", },
         # dic with only the kept k-mers
         {"AA": 0, "AC": 0, "AG": 0, "AT": 0, "CA": 0, "CC": 0, "CG": 0,
          "GA": 0, "GC": 0, "TA": 0, },
         # mapping of array indexes, for the forward strand, to where it should be added
         [0, 1, 2, 3, 4, 5, 6,  ],
         # mapping of array indexes, for the reverse complement, to where it should be added
         [], ),
]

@pytest.mark.parametrize("k, d_codons_orig_target, d_template_counts_combined, ar_codons_forward_addr, ar_codons_rev_comp_addr", initialized_arrays)
def test_table_rev_comp_to_forward_strand(k, d_codons_orig_target, *_):
    assert bio.table_rev_comp_to_forward_strand(k) == d_codons_orig_target

@pytest.mark.parametrize("k, d_codons_orig_target, d_template_counts_combined, ar_codons_forward_addr, ar_codons_rev_comp_addr", initialized_arrays)
def test_init_variables(k, d_codons_orig_target, d_template_counts_combined, ar_codons_forward_addr, ar_codons_rev_comp_addr):
    cyt_ext.init_variables(k)
    assert cyt_ext.get_d_codons_orig_target()       == d_codons_orig_target
    assert cyt_ext.get_d_template_counts_combined() == d_template_counts_combined
    assert cyt_ext.get_ar_codons_forward_addr()     == ar_codons_forward_addr
    assert cyt_ext.get_ar_codons_rev_comp_addr()    == ar_codons_rev_comp_addr
