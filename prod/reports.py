#!/usr/bin/env python3
# #############################################################################
# Sylvain @ GIS / Biopolis / Singapore
# Sylvain Jun-Zhe RIONDET <Riondet_Sylvain_from.tp@gis.a-star.edu.sg>
# Started on 2019-12-17
# Reads Binning Project
#
# #############################################################################
#
#
# About reports from Kraken2
#
import os
import os.path as osp
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc

from prod.bio import ncbi, get_list_rank
from prod.tools import PATHS

# pd.set_option('precision', 2)


class Report:

    def __init__(self, folder, filename, bin_full, nb_reads=None):
        self.folder = folder
        self.filename = filename
        self.bin_full = bin_full  # either "full" or "bins"
        self.nb_reads = nb_reads
        self.report = None
        self.all_reports = None

        self.recall = None
        self.precision = None
        self.auc = None

        self.load()

    def load(self, ):
        """Load and filter to species only"""
        if self.bin_full == "full":
            self.report = pd.read_csv(osp.join(self.folder, self.filename), sep="\t", names=[
                "%_clade_covered", "reads_clade_covered", "reads_direct_taxo", "Rank", "NCBI_Taxo_id",
                "Scientific_name"])
            self.report = self.report[self.report["Rank"] == "S"].sort_values(by=["%_clade_covered"], ascending=False)

        elif self.bin_full == "bins":
            self.aggregate_bins()
        else:
            raise NotImplementedError(f"'{self.bin_full}' isn't implemented")

    def aggregate_bins(self):
        """ Reports should have an identifier <.full.> or <.bin-x.> to identify them """

        assert self.nb_reads != 0, "The total number of reads is necessary to load and compute for bins stats"
        reports_bins = []
        for rep in self.filename:
            tmp_df = pd.read_csv(osp.join(self.folder, rep), sep="\t", names=[
                "%_clade_covered", "reads_clade_covered", "reads_direct_taxo", "Rank", "NCBI_Taxo_id",
                "Scientific_name"])
            tmp_df["bin"] = int(rep.split(".bin-")[1].split(".")[0])
            reports_bins.append(tmp_df)
        report_bins = pd.concat(reports_bins, ignore_index=True)

        self.all_reports = report_bins[report_bins.Rank == "S"].sort_values("reads_clade_covered", ascending=False)
        aggregated_bins = self.all_reports.groupby("NCBI_Taxo_id").agg('sum').sort_values("reads_clade_covered",
                                                                                          ascending=False)
        aggregated_bins["Scientific_name"] = ncbi.translate_to_names(aggregated_bins.index)
        aggregated_bins["%_clade_covered"] = aggregated_bins.reads_clade_covered / self.nb_reads * 100
        aggregated_bins.drop(columns=["bin"], inplace=True)
        aggregated_bins.reset_index(inplace=True)
        self.report = aggregated_bins

    def assigned_reads(self):
        return self.report["%_clade_covered"].sum()

    def prec_recall(self, gt_species):
        thresholds = self.report["%_clade_covered"].round(2).unique()
        thresholds.sort()
        data = []

        for i, threshold in enumerate(thresholds):
            found = set(self.report[self.report["%_clade_covered"] >= threshold].NCBI_Taxo_id)
            tp = len(set.intersection(found, gt_species))
            fn = len(gt_species) - tp
            fp = len(found) - tp
            data.append((threshold, tp, fn, fp,))

        df_auc = pd.DataFrame(data, columns=["threshold", "tp", "fn", "fp"])
        df_auc["recall"] = df_auc.tp / (df_auc.tp + df_auc.fn)
        df_auc["precision"] = df_auc.tp / (df_auc.tp + df_auc.fp)
        df_auc[["recall", "precision"]] = df_auc[["recall", "precision"]].fillna(0)
        self.recall = df_auc["recall"].tolist()
        self.precision = df_auc["precision"].tolist()
        self.auc = auc(self.recall, self.precision)

    def plot_pr(self):
        plt.plot(self.recall, self.precision)
        plt.axis([0, 1, 0, 1])

    def __repr__(self):
        return f"Report from {self.bin_full} DB classification, {self.filename}"


class ReportsAnalysis:

    def __init__(self, folder, string_full, string_bins, path_ground_truth):
        """ string* are matching string to find the full and bin reports """
        self.folder = folder
        self.string_full = string_full
        self.string_bins = string_bins
        self.path_ground_truth = path_ground_truth
        self.path_report_full = ""
        self.path_report_bins = ""

        self.selected_r = None
        self.gt = None
        self.gt_stats = None
        self.reports = {}

        self.nb_reads = None
        self.gt_species = None
        self.auc = None
        self._recall = {}
        self._precision = {}

    @property
    def report(self):
        if self.selected_r is None:
            self.selected_r = 0
        if self.selected_r in self.reports.keys():
            return self.reports[self.selected_r]
        else:
            return None

    @property
    def recall(self):
        if self.selected_r not in self._recall.keys():
            return None
        else:
            return self._recall[self.selected_r]

    @property
    def precision(self):
        if self.selected_r not in self._precision.keys():
            return None
        else:
            return self._precision[self.selected_r]

    def load_gt(self):
        self.gt = pd.read_pickle(self.path_ground_truth)
        self.nb_reads = self.gt.shape[0]
        self.gt["count"] = 1
        gt_stats = self.gt.groupby("taxon").sum()
        gt_stats.reset_index(inplace=True)
        gt_stats["species"] = get_list_rank(gt_stats.taxon)
        gt_stats = gt_stats.groupby("species").sum()[["count"]].sort_values("count", ascending=False)
        gt_stats.reset_index(inplace=True)
        gt_stats["name"] = ncbi.translate_to_names(gt_stats.species)
        gt_stats["percentage"] = gt_stats["count"] * 100 / gt_stats["count"].sum()
        self.gt_stats = gt_stats

    def load_reports(self):
        found = [f for f in os.listdir(self.folder) if self.string_full in f and "report" in f]
        assert len(
            found) == 1, f"Multiple matches ({len(found)}) for full report file, with string ({self.string_full})"
        self.reports[0] = Report(self.folder, found[0], "full")

        found = [f for f in os.listdir(self.folder) if self.string_bins in f and "report" in f]
        assert len(
            found) > 1, f"Not enough matches ({len(found)}) for bins report files, with string ({self.string_bins})"
        self.reports[1] = Report(self.folder, found, "bins", self.nb_reads)

    def prec_recall(self, select=-1):
        self.gt_species = set(self.gt_stats.species.unique())
        if select < 0:
            for k in self.reports.keys():
                self.reports[k].prec_recall(self.gt_species)
        else:
            self.reports[select].prec_recall(self.gt_species)

    def plot_pr(self):
        plt.plot(self.report.recall, self.report.precision)
        plt.axis([0, 1, 0, 1])
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"{self.report.bin_full} database, auc={self.report.auc:.2f}")


# todo: add nice call
if __name__ == '__main__':
    folder_reports = "/home/ubuntu/Disks/HDD1000/Reports/2019-12-05_100000-WindowReads_20-BacGut"
    analysis = ReportsAnalysis(folder_reports, "full", ".bin-",
                               "~/Data/Segmentation/Test-Data/Synthetic_from_Genomes/" \
                               "2019-12-05_100000-WindowReads_20-BacGut/2019-12-05_100000-WindowReads_20-BacGut.GT.pd")
    analysis.load_gt()
    analysis.load_reports()
    analysis.prec_recall()

    analysis.plot_pr()

    analysis.selected_r = 1
    analysis.plot_pr()



