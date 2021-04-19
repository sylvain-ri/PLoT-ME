#!/usr/bin/env python3
"""
First attempt to add test to project, using pytest

"""
# import pytest
from plot_me import bio


# @pytest.mark.parametrize()
def test_bio_reverse_complement_string():
    assert bio.reverse_complement_string("AACG") == "CGTT"


def test_bio_codon_addr():
    assert bio.codon_addr("AAAA") == 0
    assert bio.codon_addr("TT") == 15
