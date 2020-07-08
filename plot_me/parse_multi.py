#!/usr/bin/env python3
"""
#############################################################################
Run parse_DB.py with multiple parameters

#############################################################################
Sylvain @ GIS / Biopolis / Singapore
Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
Started on 2019-12-11
Reads Binning Project
#############################################################################
"""

import argparse
import datetime as dt
from multiprocessing.pool import Pool
import pandas as pd
from tqdm import tqdm

import parse_DB
from tools import ArgumentParserWithDefaults


def main(k, cores, skip_existing):
    """ record which param have been done
        date	cluster	k	w	clf_param	omit
    """
    history = "/home/ubuntu/Scripts/Reads_Binning/plot_me/tried_param.tsv"

    # with Pool(8) as pool:  # file copy don't need many cores (main.cores)
    #     results = list(tqdm(pool.imap(pll_copy_segments_to_bin, islice(df_per_fna, stop if stop > 0 else None)),
    #                         total=len(df_per_fna)))

    n_clusters = 20
    classifier_param = ["kraken2", "35", "31", "7"]  # ["kraken2", "25", "22", "4"]
    omit_folders = ("plant", "vertebrate")

    for k in (k, ):
        for w in (5000, 10000, 25000):
            omit = ",".join(omit_folders)
            f_clf_param = ",".join(classifier_param)
            s_param = f"{dt.datetime.now():%Y-%m-%d_%H-%M}\t{n_clusters}\t{k}\t{w}\t{f_clf_param}\t{omit}"
            print(s_param)

            # skip already done
            df = pd.read_csv(history, sep="\t")
            if df[(df.k == k) & (df.w == w) & (df.cluster == n_clusters) & (df.clf_param == f_clf_param)].shape[0] > 0:
                print("skipped")
                continue

            parse_DB.main(
                folder_database="/mnt/data/NCBI/20190704/refseq",
                folder_output="/mnt/data/Segmentation",
                n_clusters=n_clusters,
                k=k,
                window=w,
                cores=cores,
                skip_existing=skip_existing,
                force_recount=False,
                early_stop=6,
                omit_folders=omit_folders,
                path_taxonomy="/mnt/data/taxonomy",
                ml_model="minikm",
                full_DB=False,
                classifier_param=classifier_param,
                k2_clean=False)

            print("done")
            with open(history, 'a') as f:
                f.write(s_param + "\n")



if __name__ == '__main__':
    # Option to display default values, metavar='' to remove ugly capitalized option's names
    parser = ArgumentParserWithDefaults(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-k', '--kmer',      default=4, type=int, help='Size of the kmers', metavar='')
    parser.add_argument('-c', '--cores',     default=4, type=int, help='number of cores', metavar='')
    parser.add_argument('-s', '--skip_existing',     default="110000", type=str, help='skip steps', metavar='')
    args = parser.parse_args()

    main(args.kmer, args.cores, args.skip_existing)
