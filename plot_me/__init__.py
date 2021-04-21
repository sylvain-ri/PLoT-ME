#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#############################################################################
Package script, testing who python works

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

from datetime import datetime
from pathlib import Path

PLOT_ME_ROOT = Path.home().joinpath("tmp/PLoT-ME")
LOGS = PLOT_ME_ROOT.joinpath(f"logs/{datetime.now():%Y-%m-%d_%H-%M}.log")
LOGS.parent.mkdir(parents=True, exist_ok=True)
RECORDS = PLOT_ME_ROOT.joinpath(f"timings/classify_{datetime.now():%Y-%m-%d_%H-%M}.tsv")

# !!  Need to initialize variables before importing modules  !!
from plot_me import parse_DB, classify, tools, bio


if __name__ == '__main__':
    print(f"PLoT-ME : Pre-classification of Long-reads for Memory Efficient Taxonomic assignment "
          f"Information in the readme file and on https://github.com/sylvain-ri/PLoT-ME. ")

