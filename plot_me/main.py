#!/usr/bin/env python3
"""
#############################################################################
Main script, with all packages

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-18
Reads Binning Project
#############################################################################
               FLOW
todo: SyntheticReads from RefSeq DB
      - go through SOME files
      - load genome, make CustomRead, split the sequence [start:stop],
      - write it to output fastQ

todo: parse RefSeq DB
      - go through all files
      - load fastQ (genome), make "CustomReads", split them in windows, count k-mers,
      - save frequency into similar file structure, only .fna => .4mer

todo: find bins for DB
      - go through all frequency counts
      - load all k-mer freq into one file
      - cluster all windows
      - write bin number
      - windows ???
todo: build kraken DB
      - use kraken2 builder

todo: microbiome community
      - load fastQ file, make CustomeRead for each read
      - for each Read: Reads.bin(), Reads.write(path_bin_x)
      - for each new FastQ binned file
          - fastQ.classify()
          - fastQ.understand_report()
      - load all reports
      - AUC()

#############################################################################
         Classes
CustomRead: additional functions for Read (split, count_kmer, change description)
FilesInDir: Go through all files
FastQ     : load fastQ (genome or community), read all reads/chromosomes,
            do for each: bin()
Report    : combine, compute AUC
"""

# #############################################################################
# Imports
import os
import os.path as osp

# Script






















