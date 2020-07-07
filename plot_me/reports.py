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

# import ete3.ncbi_taxonomy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import auc

from plot_me.bio import ncbi, get_list_rank
from plot_me.tools import PATHS


pd.set_option('precision', 5)
# ncbi = ete3.ncbi_taxonomy.NCBITaxa()


class Report:
    kraken_cols = ["per_clade_covered", "reads_clade_covered", "reads_direct_taxo", "rank", "taxon", "scientific_name"]
    line_style = ['-.', '--', ':']
    obj_counter = 0

    def __init__(self, title, folder, nb_reads=None):
        self.obj_id = Report.obj_counter
        Report.obj_counter += 1
        self.title       = title
        self.folder      = folder
        self.nb_reads    = nb_reads
        self.nb_assigned = None
        self.report      = None
        self.all_reports = None

        self.thresholds  = None
        self.recall      = None
        self.precision   = None
        self.auc         = None

    def load_full(self, filename):
        """Load and filter to species only"""
        self.report = pd.read_csv(osp.join(self.folder, filename), sep="\t", names=self.kraken_cols)
        self.nb_reads = self.report[self.report.taxon <= 1].reads_clade_covered.sum()
        self.report = self.report[self.report["rank"] == "S"][["reads_clade_covered", "taxon"]].groupby(["taxon"]).sum()
        #         self.report["cluster"] = "full"
        self.report.rename(columns={"reads_clade_covered": "full_DB"}, inplace=True)
        self.assigned_reads()

    def load_multi(self, list_files):
        """ Reports should have an identifier <.bin-x.> to identify them """

        reports_bins = []
        bins = []
        list_total_nb_reads = []
        for file_name in list_files:
            tmp_df = pd.read_csv(osp.join(self.folder, file_name), sep="\t", names=self.kraken_cols)
            list_total_nb_reads.append(tmp_df[tmp_df.taxon <= 1].reads_clade_covered.sum())
            tmp_df["cluster"] = file_name.split(".bin-")[1].split(".")[0]
            reports_bins.append(tmp_df)
        print(self.title, list_total_nb_reads)
        self.nb_reads = sum(list_total_nb_reads)
        report_bins = pd.concat(reports_bins, ignore_index=True)
        report_selected = report_bins[report_bins["rank"] == "S"][["reads_clade_covered", "taxon", "cluster"]]
        # todo: losing the cluster provenance by summing everything :/
        aggregated = report_selected.groupby(["taxon"]).sum()
        aggregated.rename(columns={"reads_clade_covered": self.title}, inplace=True)
        #         aggregated.sort_values("reads_clade_covered", ascending=False, inplace=True)
        self.report = aggregated
        self.assigned_reads()

    def load_gt(self, file_path):
        gt_tmp = pd.read_pickle(file_path)
        gt_counting = pd.DataFrame(gt_tmp.taxon.value_counts())
        gt_counting.rename(columns={"taxon": "ground_truth"}, inplace=True)
        gt_counting["taxon"] = get_list_rank(gt_counting.index)
        self.report = gt_counting.groupby(["taxon"]).sum()
        self.assigned_reads()

    #         self.report["cluster"] = "gt"

    def assigned_reads(self):
        self.nb_assigned = int(self.report.iloc[:, 0].sum())

    def normalize(self):
        self.report.iloc[:, 0] /= self.nb_assigned

    def prec_recall(self, gt_species):
        """ Get a set containing the species present. change to ratio instead of absolute number after working on the report """
        #         print(self.title)
        # Floor the numbers by multiplying by 'rounding', then floor with numpy, then dividing again.
        rounding = 10 ** 5
        thresholds = ((self.report.iloc[:, 0] * rounding).apply(np.floor) / rounding).unique()
        thresholds.sort()
        data = []

        for i, threshold in enumerate(thresholds):
            found = set(self.report[self.report.iloc[:, 0] >= threshold].index)
            tp = len(set.intersection(found, gt_species))
            fn = len(gt_species) - tp
            fp = len(found) - tp
            data.append((threshold, tp, fn, fp,))

        df_auc = pd.DataFrame(data, columns=["threshold", "tp", "fn", "fp"])
        df_auc["recall"] = df_auc.tp / (df_auc.tp + df_auc.fn)
        df_auc["precision"] = df_auc.tp / (df_auc.tp + df_auc.fp)
        df_auc[["recall", "precision"]] = df_auc[["recall", "precision"]].fillna(0)
        # Extend the last precision to 0 recall, as we don't have abundance threshold down to 0%
        df_auc.loc[df_auc.index.max() + 1] = df_auc.iloc[-1]
        df_auc.recall.iloc[-1] = 0

        self.df_auc = df_auc
        self.thresholds = thresholds
        self.recall = df_auc["recall"].tolist()
        self.precision = df_auc["precision"].tolist()
        self.auc = auc(self.recall, self.precision)

    def plot_pr(self, nb=5, total=10, string_gt=""):
        # todo: thicker line, dotted line
        ratio = (total - nb) / total
        label = f"auc={self.auc:.3f}, ({self.nb_assigned}/{self.nb_reads}) : {self.title}"
        plt.plot(self.recall, self.precision,  # alpha=0.7,
                 linewidth=2 + 3 * ratio, linestyle=self.line_style[self.obj_id % len(self.line_style)],
                 marker='+', markersize=10 + 4 * ratio, markeredgewidth=1 + 2 * ratio,
                 label=label)  # Change line style so see them despite overlaying each other
        self.legend = label

    def __repr__(self):
        return f"Report from {self.title} DB classification, {self.folder}"


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


def compare_setups(reports):
    comparison = pd.concat(reports, axis=1)
    comparison.fillna(0, inplace=True)                   # non found taxon are NaN
    comparison = comparison[(comparison.T != 0).any()]   # Remove full zeros rows
    comparison["species_name"] = ncbi.get_taxid_translator(comparison.index).values()
#     comparison.drop(columns=["cluster"], inplace=True)
    return comparison


def load_all(folder_reports, path_gt, settings=["3mer_s5", "3mer_s10", "4mer_s10", ]):
    gt = Report("GT", folder_reports)
    gt.load_gt(path_gt)
    gt.normalize()

    reports = {}
    full_matches = [f for f in os.listdir(folder_reports) if "full" in f and f.endswith(".report")]
    if len(full_matches) >= 1:
        full_k25 = [f for f in full_matches if "k25_l22_s4" in f]
        full_k35 = [f for f in full_matches if "kraken2_full" in f]
        if len(full_k25) > 0:
            full = Report("full_k25", folder_reports, )
            full.load_full(full_k25[0])
            full.normalize()
            reports["full_k25_l22_s4"] = full
        if len(full_k35) > 0:
            full = Report("full_k35", folder_reports, )
            full.load_full(full_k35[0])
            full.normalize()
            reports["full_k35_l31_s7"] = full
    else:
        print(f"wrong number of files matching for full report: {len(full_matches)}")

    #     params = ["3mer_s5", "3mer_s10", "4mer_s10", ] # "3mer_s10",
    for param in settings:
        reports[param] = Report(param, folder_reports, )
        reports[param].load_multi([f for f in os.listdir(folder_reports)
                                   if f"clustered_by_minikm_{param}000" in f and f.endswith(".report")])
        reports[param].normalize()

        k, s = param.split("mer_s")
        files_k25 = [f for f in os.listdir(folder_reports) if f"minikm_b10_k{k}_s{s}000" in f and f.endswith(".report")]
        if len(files_k25) > 0:
            key_p = param + "_k25"
            reports[key_p] = Report(key_p, folder_reports, )
            reports[key_p].load_multi(files_k25)
            reports[key_p].normalize()

    comparison = compare_setups([gt.report] + [reports[k].report for k in reports.keys()])
    return comparison, reports, gt


def plot_auc_comparison(gt, reports, fig_path, title=f"AUC for various parameters of k and w"):
    # Add the ground truth species
    r = gt.report
    r = r[r.ground_truth > 0].copy()
    r["percentage"] = round(r.ground_truth * 100, 2)
    r["species"] = ncbi.translate_to_names(r.index)
    r.drop(columns=["ground_truth"], inplace=True)
    r.sort_values(by=["percentage", "species"], ascending=[False, True], inplace=True)
    string_gt = "\n * GroundTruth * \n" + r[r.percentage > 0].to_string(index=False)

    legend = []
    gt_set = set(gt.report.index)
    fig = plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 12})
    for i, k in enumerate(reports.keys()):
        reports[k].prec_recall(gt_set)
        reports[k].plot_pr(i, len(reports))
        legend.append(reports[k].legend)

    plt.plot([], [], ' ', label=string_gt)
    legend.append(string_gt)
    #     plt.legend(legend, loc='lower left')
    plt.legend(legend, title="AUC, (classified @species/total nb of reads), parameter", loc='center left',
               bbox_to_anchor=(1, 0.5))
    plt.axis([0, 1.1, 0, 1.1])
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(title)
    plt.savefig(osp.join(fig_path, f"{title}"), bbox_inches='tight')


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



