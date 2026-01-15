from datetime import datetime
import json
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
import random # <- added on 250506
# marked library (<-) added on 250408 to interpolate distributions using cubic spline rather than polynomial fit
from scipy.interpolate import CubicSpline # <-
import scipy.stats as stats
from sklearn.metrics import mean_squared_error
import statsmodels.api as sm
import string
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")
import time

from mask_predictions import *
from shared_utils import *

seed_time = time.time()
random.seed(0)


class Analyze_predictions_utils(Mask_predictions):
    def __init__(self, params, mutation_stats):
        super().__init__(params, mutation_stats)

        self.params             = params

        self.base_dict          = None
        self.blocks             = None
        self.expression_stats   = None
        self.global_predictions = None
        self.masking_stats      = pd.DataFrame()
        self.means              = {"gene id": [], "mean": []}
        self.mutation_stats     = mutation_stats
        self.predictions        = None
        self.seqs               = None

        if "id_filter" in self.params and self.params["id_filter"] == "CPTAC3":
            self.params["selection_cols"] = {"gene symbol": "ID:gene symbol", "project": "ID:project", "case id": "ID:case id", "sample id": "ID:sample id"}
            self.projects                 = ["Adenomas and Adenocarcinomas", "Ductal and Lobular Neoplasms", "Epithelial Neoplasms, NOS", "Gliomas", "Lipomatous Neoplasms",
                                             "Nevi and Melanomas", "Squamous Cell Neoplasms", "Transitional Cell Papillomas and Carcinomas"]
        
        elif "id_filter" in self.params and self.params["id_filter"] == "MSK":
            self.params["selection_cols"] = {"gene symbol": "ID:gene symbol", "sample id": "ID:sample id", "cancer type": "ID:cancer type"}
            self.projects                 = ["Breast Cancer", "Colorectal Cancer", "Non-Small Cell Lung Cancer", "Pancreatic Cancer", "Prostate Cancer"]

        elif "id_filter" in self.params and self.params["id_filter"] == "TCGA":
            self.params["selection_cols"] = {"gene symbol": "ID:gene symbol", "project": "ID:project", "case id": "ID:case id", "sample id": "ID:sample id"}
            self.projects                 = ["TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
                                             "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
                                             "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
                                             "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"]
            
        self.stats              = None
        self.variant_filter     = None
        self.variants           = None

        # plot settings
        self.cmap = plt.get_cmap('viridis')

        rcParams['font.family']         = 'sans-serif'
        rcParams['font.sans-serif']     = ['Arial']
        rcParams["axes.labelsize"]      = 13
        rcParams["xtick.labelsize"]     = 11
        rcParams["ytick.labelsize"]     = 11

        if params["create_newdir"] == True:
            self.newdir           = self.create_newdir()
            self.print_params()


    def _analyze_adaptation(self, rules_analysis, exon_probs, exon_pos, exon_scores, index, total_cds, category, side, show=False):
        if index[1] <= index[0]: return rules_analysis
        rules_analysis[category]["expected_"+side] += np.nansum(exon_probs[index[0]:index[1]])
        rules_analysis[category]["observed_"+side] += np.sum([1 for i in exon_pos if i-total_cds >= index[0] and i-total_cds < index[1]])
        rules_analysis[category][side+"_size"]     += index[1]-index[0]
        rules_analysis[category][side+"_preds"].extend([exon_scores[k] for k in range(len(exon_scores)) if exon_pos[k] >= index[0] and exon_pos[k] < index[1]])
        if show == True: print(category, {key: rules_analysis[category][key] for key in rules_analysis[category] if "preds" not in key})
        return rules_analysis


    def analyze_adaptation(self):
        stop_codon_distribution = [0 for _ in range(100)]
        start_cutoffs           = [0, 50, 100, 150, 200, 250, 300, 500]
        ptc_distribution        = {"occurrence": [0 for _ in range(100)], "susceptibility": [[] for _ in range(100)]}
        rules_analysis          = {"canonical": {"expected_escape": 0, "expected_target": 0, "observed_escape": 0, "observed_target": 0, "escape_size": 0, "target_size": 0, "escape_preds": [], "target_preds": [],
                                                 "binomial statistic": None, "binomial pvalue": None, "fisher exact pvalue": None, "fisher exact statistic": None,
                                                 "median escape": None, "median target": None, "mw statistic": None, "mw pvalue": None},
                                   **{"lindeboom_"+str(start_cutoff): {"expected_escape": 0, "expected_target": 0, "observed_escape": 0, "observed_target": 0, "escape_size": 0, "target_size": 0, "escape_preds": [], "target_preds": [],
                                                                       "binomial statistic": None, "binomial pvalue": None, "fisher exact pvalue": None, "fisher exact statistic": None,
                                                                       "median escape": None, "median target": None, "mw statistic": None, "mw pvalue": None}
                                                                       for start_cutoff in start_cutoffs}}
        
        total_cds_test = 0
        bar = IncrementalBar(set_bar("analyzing adaptation"), max=self.predictions.shape[0])
        for i in range(self.predictions.shape[0]):
            for exon in self.predictions.iloc[i].loc["variant targets"]["LABEL:NMD score"]:
                for j in range(len(self.predictions.iloc[i].loc["variant targets"]["LABEL:NMD score"][exon])):
                    if pd.isna(self.predictions.iloc[i].loc["variant targets"]["LABEL:NMD score"][exon][j]) == False:
                        ptc_pos = int(self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"][exon][j])
                        pos     = int(100*(ptc_pos/ len(self.predictions.iloc[i].loc["probabilities"])))
                        ptc_distribution["occurrence"][pos] += 1
                        ptc_distribution["susceptibility"][pos].append(self.predictions.iloc[i].loc["variant targets"]["LABEL:NMD score"][exon][j])

            for j in range(len(self.predictions.iloc[i].loc["probabilities"])):
                if pd.isna(self.predictions.iloc[i].loc["probabilities"][j]) == False:
                    pos = int(100*(j / len(self.predictions.iloc[i].loc["probabilities"])))
                    stop_codon_distribution[pos] += self.predictions.iloc[i].loc["probabilities"][j]

            total_cds = 0
            for j, exon in enumerate(self.predictions.iloc[i].loc["probabilities_by_exon"]):
                exon_probs  = self.predictions.iloc[i].loc["probabilities_by_exon"][exon]
                exon_pos    = self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"][exon]
                exon_scores = self.predictions.iloc[i].loc["variant targets"]["LABEL:NMD score"][exon]

                # canonical rules
                if j < len(self.predictions.iloc[i].loc["probabilities_by_exon"])-2:
                    rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, len(exon_probs)), total_cds, "canonical", "target")
                
                # penultimate exon
                if j == len(self.predictions.iloc[i].loc["probabilities_by_exon"])-2:
                    rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, max(0, len(exon_probs)-50)), total_cds, "canonical", "target")
                    rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (max(0, len(exon_probs)-50), len(exon_probs)), total_cds, "canonical", "escape")

                # last exon
                if j == len(self.predictions.iloc[i].loc["probabilities_by_exon"])-1:
                    rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, len(exon_probs)), total_cds, "canonical", "escape")
                
                if total_cds_test+len(exon_probs) != rules_analysis["canonical"]["escape_size"]+rules_analysis["canonical"]["target_size"]:
                    print("< size error1 occurred @analyze_adaptation:", total_cds_test, "/", rules_analysis["canonical"]["escape_size"]+rules_analysis["canonical"]["target_size"])
                    exit()

                # lindeboom rules
                for start_cutoff in start_cutoffs:
                    if j < len(self.predictions.iloc[i].loc["probabilities_by_exon"])-2:
                        if total_cds < start_cutoff:
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, min(start_cutoff-total_cds, len(exon_probs))),
                                                                      total_cds, "lindeboom_"+str(start_cutoff), "escape")
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (min(start_cutoff-total_cds, len(exon_probs)),
                                                                      len(exon_probs)), total_cds, "lindeboom_"+str(start_cutoff), "target")
                            
                        if total_cds >= start_cutoff:
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, len(exon_probs)),
                                                                      total_cds, "lindeboom_"+str(start_cutoff), "target")

                    # penultimate exon
                    if j == len(self.predictions.iloc[i].loc["probabilities_by_exon"])-2:
                        if total_cds < start_cutoff:
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, min(start_cutoff-total_cds, len(exon_probs))),
                                                                      total_cds, "lindeboom_"+str(start_cutoff), "escape")
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (min(start_cutoff-total_cds, len(exon_probs)),
                                                                      max(0, len(exon_probs)-50)), total_cds, "lindeboom_"+str(start_cutoff), "target")
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (max(min(start_cutoff-total_cds, len(exon_probs)),
                                                                      max(0, len(exon_probs)-50)), len(exon_probs)), total_cds, "lindeboom_"+str(start_cutoff), "escape")

                            if len(exon_probs)-max(min(start_cutoff-total_cds, len(exon_probs)), max(0, len(exon_probs)-50)) > 50:
                                print("< size error2 occurred with", max(min(start_cutoff-total_cds, len(exon_probs))), "/", len(exon_probs))
                                exit()

                        if total_cds >= start_cutoff:
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, max(0, len(exon_probs)-50)),
                                                                      total_cds, "lindeboom_"+str(start_cutoff), "target")
                            rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (max(0, len(exon_probs)-50), len(exon_probs)),
                                                                      total_cds, "lindeboom_"+str(start_cutoff), "escape")

                    # last exon
                    if j == len(self.predictions.iloc[i].loc["probabilities_by_exon"])-1:
                        rules_analysis = self._analyze_adaptation(rules_analysis, exon_probs, exon_pos, exon_scores, (0, len(exon_probs)), total_cds,
                                                                  "lindeboom_"+str(start_cutoff), "escape")
                    
                    if (total_cds_test+len(exon_probs) != rules_analysis["lindeboom_"+str(start_cutoff)]["escape_size"]
                        +rules_analysis["lindeboom_"+str(start_cutoff)]["target_size"]):
                        print("< size error3 occurred @analyze_adaptation:", total_cds_test, "/", rules_analysis["lindeboom_"+str(start_cutoff)]["escape_size"]
                              +rules_analysis["lindeboom_"+str(start_cutoff)]["target_size"])
                        exit()
                        
                total_cds      += len(exon_probs)
                total_cds_test += len(exon_probs)

            #bar.next()
        bar.finish()

        # analyze trajectory
        ptc_distribution["susceptibility"] = [np.mean(susceptibilities) if len(susceptibilities) > 0 else 0 for susceptibilities in ptc_distribution["susceptibility"]]
        corr_ptc_occurrence                = [value / np.sum(ptc_distribution["occurrence"]) for value in ptc_distribution["occurrence"]]
        corr_stop_codon_distribution       = [value / np.sum(stop_codon_distribution) for value in stop_codon_distribution]
        corr_susceptibilities              = [value / np.sum(ptc_distribution["susceptibility"]) for value in ptc_distribution["susceptibility"]]

        adaptation_analysis = pd.DataFrame({"susceptibility": ptc_distribution["susceptibility"], "rel. susceptibility": corr_susceptibilities,
                                            "stop codon occurrence": stop_codon_distribution, "rel. stop codon occurrence": corr_stop_codon_distribution,
                                            "ptc occurrence": ptc_distribution["occurrence"], "rel. ptc occurrence": corr_ptc_occurrence})
        
        adaptation_analysis.to_csv(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_adaptation_analysis.txt", index=False, sep=",")

        plt.plot(corr_stop_codon_distribution, label="stop codon candidates")
        plt.plot(corr_ptc_occurrence, label="observed PTCs")
        plt.plot(corr_susceptibilities, label="mean susceptibility")
        plt.legend()
        plt.show()

        # analyze categories
        for key in rules_analysis:
            total  = rules_analysis[key]["observed_escape"]+rules_analysis[key]["observed_target"]
            fisher = stats.fisher_exact([[rules_analysis[key]["expected_escape"]*total, rules_analysis[key]["expected_target"]*total],
                                         [rules_analysis[key]["observed_escape"], rules_analysis[key]["observed_target"]]])
            rules_analysis[key]["fisher exact statistic"] = fisher.statistic
            rules_analysis[key]["fisher exact pvalue"]    = fisher.pvalue

            # binomial test
            expected_probability = rules_analysis[key]["expected_escape"] / (rules_analysis[key]["expected_escape"]+rules_analysis[key]["expected_target"])
            observed_probability = rules_analysis[key]["observed_escape"] / (rules_analysis[key]["observed_escape"]+rules_analysis[key]["observed_target"])       
            rules_analysis[key]["binomial statistic"] = observed_probability-expected_probability

            if observed_probability >= expected_probability:
                rules_analysis[key]["binomial pvalue"] = 1-stats.binom.cdf(rules_analysis[key]["observed_escape"]-1,
                                                                          (rules_analysis[key]["observed_escape"]+rules_analysis[key]["observed_target"]), expected_probability)

            if observed_probability < expected_probability:
                rules_analysis[key]["binomial pvalue"] = stats.binom.cdf(rules_analysis[key]["observed_escape"],
                                                                        (rules_analysis[key]["observed_escape"]+rules_analysis[key]["observed_target"]), expected_probability)
                
            # Mann-Whitney U-test
            rules_analysis[key]["median escape"] = np.median(rules_analysis[key]["escape_preds"])
            rules_analysis[key]["median target"] = np.median(rules_analysis[key]["target_preds"])
            mw                                   = stats.mannwhitneyu(rules_analysis[key]["escape_preds"], rules_analysis[key]["target_preds"])
            rules_analysis[key]["mw statistic"]  = mw.statistic
            rules_analysis[key]["mw pvalue"]     = mw.pvalue

        for key1 in rules_analysis:
            for key2 in rules_analysis[key1]:
                if type(rules_analysis[key1][key2]) != list:
                    print(key1, key2, rules_analysis[key1][key2])

        with open(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_rules_analysis.txt", "w") as f:
            f.write(json.dumps(rules_analysis, indent=4))


    def analyze_blocks(self, status):
        block_stats   = {**{block_target: [] for block_target in self.params["block_targets"]}, **{"mean": [], "std": [], "count": [], "mean_count": []}}
        self.variants = self.variants.sort_values(by=self.params["block_targets"])
        blocks        = create_blocks(self.variants, self.params["block_targets"])

        bar = IncrementalBar(set_bar("analyzing blocks"), max=len(blocks))
        for i in range(len(blocks)):
            selected_block = self.variants.iloc[blocks[i]["index"]]
            
            for block_target in self.params["block_targets"]:
                block_stats[block_target].append(blocks[i]["block id"][block_target])

            block_stats["mean"].append(np.mean(selected_block["FEATURE:prediction"]))
            block_stats["std"].append(np.std(selected_block["FEATURE:prediction"]))
            block_stats["count"].append(selected_block["FEATURE:prediction"].shape[0])
            block_stats["mean_count"].append(selected_block["FEATURE:prediction"].shape[0])

            bar.next()
        bar.finish()

        block_stats = pd.DataFrame(block_stats)
        block_stats.to_csv(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_block_stats.txt", index=False, sep=",")

        # create block stats for first block target only (e.g. project)
        outer_block_stats   = {self.params["block_targets"][0]: [], "mean": [], "std": [], "count": [], "mean_count": [], "cases": []}
        outer_block_targets = self.variants.drop_duplicates(subset=self.params["block_targets"][0])[self.params["block_targets"][0]].tolist()
        
        bar = IncrementalBar(set_bar("analyzing outer blocks"), max=len(outer_block_targets))
        for outer_block_target in outer_block_targets:
            selected_block_stats = block_stats[block_stats[self.params["block_targets"][0]] == outer_block_target]

            outer_block_stats[self.params["block_targets"][0]].append(outer_block_target)
            outer_block_stats["mean"].append(np.mean(selected_block_stats["mean"]))
            outer_block_stats["std"].append(np.std(selected_block_stats["mean"]))
            outer_block_stats["count"].append(np.sum(selected_block_stats["count"]))
            outer_block_stats["mean_count"].append(selected_block_stats["mean"].shape[0])
            outer_block_stats["cases"].append(self.get_cases(status, outer_block_target))

            bar.next()
        bar.finish()

        outer_block_stats = pd.DataFrame(outer_block_stats)
        outer_block_stats.to_csv(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_outer_block_stats.txt", index=False, sep=",")
        return


    def _analyze_selection(self, genes, nmd_cutoffs, project_filter=None):
        project_tag = ""
        if project_filter != None:
            project_tag = "".join(project+"-" if i < len(project_filter)-1 else project for i, project in enumerate(project_filter))
            project_tag = "_" + project_tag

        new_cols = {**{col+project_tag: [] for col in ["escape", "target"]}, **{"relative " + col+project_tag: [] for col in ["escape", "target"]}}

        for i in range(genes.shape[0]):
            real_predictions = json.loads(genes.iloc[i].loc["real predictions"])
            # filter real_predictions accoring to project filter (if passed)
            if project_filter != None:
                # exception to exclude "total" section
                if "project" in genes and type(genes.iloc[i].loc["project"]) == str:
                    projects         = json.loads(genes.iloc[i].loc["project"].replace("'", "\""))
                    real_predictions = [real_prediction for j, real_prediction in enumerate(real_predictions) if projects[j] in project_filter]

            # testing against pre-defined nmd cutoffs
            if self.params["selection_mode"] == "absolute":
                escape = len([j for j in range(len(real_predictions)) if real_predictions[j] < nmd_cutoffs[0]])
                target = len([j for j in range(len(real_predictions)) if real_predictions[j] > nmd_cutoffs[1]])

            # testing against the individual average
            if self.params["selection_mode"] == "relative":
                if self.params["selection_method"] == "binomial":
                    escape = len([j for j in range(len(real_predictions)) if genes.iloc[i].loc["binomial-statistic "+self.params["variant_test_targets"][0]] > 0])
                    target = len([j for j in range(len(real_predictions)) if genes.iloc[i].loc["binomial-statistic "+self.params["variant_test_targets"][0]] < 0])

                else:
                    escape = len([j for j in range(len(real_predictions)) if real_predictions[j] < genes.iloc[i].loc["avg. model prediction"]])
                    target = len([j for j in range(len(real_predictions)) if real_predictions[j] > genes.iloc[i].loc["avg. model prediction"]]) # <- until here


            new_cols["escape"+project_tag].append(escape)
            new_cols["target"+project_tag].append(target)
            
            if len(real_predictions) > 0:
                new_cols["relative escape"+project_tag].append(float(escape)/float(len(real_predictions)))
                new_cols["relative target"+project_tag].append(float(target)/float(len(real_predictions)))
            
            else:
                new_cols["relative escape"+project_tag].append(0)
                new_cols["relative target"+project_tag].append(0)
            
        for col in new_cols:
            if col in genes: genes[col] = new_cols[col]
            else:            genes.insert(genes.shape[1], col, new_cols[col])

        return genes


    def analyze_selection(self, genes, block=None):
        # check whether total block is last entry of genes
        for col in genes:
            if "pvalue" in col:
                print(col, genes[col][genes[col].isna()].shape[0], genes[col][0:genes[col].shape[0]-1].min(),
                      genes[col][0:genes[col].shape[0]-1].max(), genes[col][0:genes[col].shape[0]-1].mean())
                
        # <- test shifted from below on 251115
        if genes.shape[0] > 0 and genes.iloc[-1].loc["block id"] != "total":
            print("< total value not found.")
            exit()

        # <- added on 251021 to allow exclusion of low-counting variants
        if self.params["selection_minsize"] != None:
            init_shape = genes.shape[0]
            genes      = genes[genes["avg. real size"] >= self.params["selection_minsize"]]
            print("< filtering using selection_minsize reduced dataset from", init_shape, "to", genes.shape[0])

        if genes.shape[0] > 0 and genes.iloc[-1].loc["block id"] == "total": # <- "genes.shape[0] > 0 and" added on 251104 from 
            # create false-discovery-rate-corrected pvalues
            for col in genes:
                if "pvalue" in col:
                    # currently scipy function not working (probably scipy error, 241026)
                    # genes[col] = scipy.stats.false_discovery_control(genes[col])
                    # import was difficult (241026), probably due to old Python version
                    # .api needed to be used (import statsmodel.api as sm instead of statsmodel as sms)
                    # sm.stats.fdrcorrection was used instead of sm.stats.multitest.fdrcorrection
                    # total value must be excluded from FDR correction
                    index      = [i for i in range(genes.shape[0]-1) if pd.isna(genes.iloc[i].loc[col]) == False] 
                    pvalues    = [genes.iloc[i].loc[col] for i in range(genes.shape[0]-1) if pd.isna(genes.iloc[i].loc[col]) == False]
                    pvalues    = sm.stats.fdrcorrection(pvalues, alpha=0.05)[1]
                    genes[col] = [*[None for _ in range(genes.shape[0]-1)], genes.iloc[-1].loc[col]]
                    
                    for i, mapped_i in enumerate(index):
                        genes.at[genes.index[mapped_i], col] = pvalues[i]

        for col in genes:
            if "pvalue" in col:
                print(col, genes[col][0:genes[col].shape[0]-1].min(), genes[col][0:genes[col].shape[0]-1].max(), genes[col][0:genes[col].shape[0]-1].mean())

        # filter for genes below significance threshold
        init_shape = genes.shape[0]

        if "escape" not in self.params["selection_method"] and "target" not in self.params["selection_method"]:
            genes = genes[genes[self.params["selection_method"]+"-pvalue "+self.params["variant_test_targets"][0]] <= self.params["variant_test_cutoff"]]

        else:
            escape_genes = genes[genes[self.params["selection_method"].replace("target", "escape")+"-pvalue "+self.params["variant_test_targets"][0]] <= self.params["variant_test_cutoff"]]
            target_genes = genes[genes[self.params["selection_method"].replace("escape", "target")+"-pvalue "+self.params["variant_test_targets"][0]] <= self.params["variant_test_cutoff"]]
            genes        = pd.concat([escape_genes, target_genes])
            genes        = genes.drop_duplicates(subset="gene symbol")
            
        print("< applying selection cutoff", self.params["variant_test_cutoff"], "reduced data from", init_shape, "to", genes.shape[0])
        
        if genes.shape[0] > 0:
            genes = self._analyze_selection(genes, [self.params["nmd_escape_cutoff"], self.params["nmd_target_cutoff"]])

            if block != None: print("<", block)

            if block == None:
                bar = IncrementalBar(set_bar("analyzing selection"), max=len(self.projects))
                for project in self.projects:
                    genes = self._analyze_selection(genes, [self.params["nmd_escape_cutoff"], self.params["nmd_target_cutoff"]], [project])
                    bar.next()
                bar.finish()

            self.show_selection(genes.sort_values(by="escape", ascending=False).iloc[0:10], "escape", project_mode=False)
            self.show_selection(genes.sort_values(by="target", ascending=False).iloc[0:10], "target", project_mode=False)

            if block == None:
                genes.to_csv(self.newdir+self.params["os_sep"]+"selection_stats.txt", index=False, sep=",")
                # generate escape/target specific outputs
                escape = genes[[True if genes.iloc[i].loc["escape"] > genes.iloc[i].loc["target"] else False for i in range(genes.shape[0])]]
                escape.to_csv(self.newdir+self.params["os_sep"]+"selection_stats_escape.txt", index=False, sep=",")
                target = genes[[True if genes.iloc[i].loc["escape"] < genes.iloc[i].loc["target"] else False for i in range(genes.shape[0])]]
                target.to_csv(self.newdir+self.params["os_sep"]+"selection_stats_target.txt", index=False, sep=",")

            if block != None:
                genes.to_csv(self.newdir+self.params["os_sep"]+"selection_stats_"+block+".txt", index=False, sep=",")

        else:
            print("< significance filtering removed all entries.")
    

    def _analyze_stop_codons(self):
        bases = ["A", "C", "G", "T"]

        tripletts = []
        for base1 in bases:
            for base2 in bases:
                for base3 in bases:
                    tripletts.append(base1+base2+base3)

        codon_efficiencies = {triplett: [] for triplett in tripletts}

        bar = IncrementalBar(set_bar("analyzing codons-resolved efficiencies"), max=self.variants.shape[0])
        for i in range(self.variants.shape[0]):
            # determine sequence position from CDS_position
            if "-" not in self.variants.iloc[i].loc["ID:CDS position"].split("/")[0]:
                seq_pos = int(self.variants.iloc[i].loc["ID:CDS position"].split("/")[0])-1

            else:
                seq_pos = self.variants.iloc[i].loc["ID:CDS position"].split("/")[0]
                seq_pos = int(seq_pos.split("-")[0])-1
                
            if len(get_numbers(self.variants.iloc[i].loc["ID:HGVSc"])) > 0 and get_numbers(self.variants.iloc[i].loc["ID:HGVSc"])[0]-1 == seq_pos:
                triplett = self.variants.iloc[i].loc["FEATURE:all cds"][seq_pos-1:seq_pos+2]
                codon_efficiencies[triplett].append(self.variants.iloc[i].loc["FEATURE:prediction"])

            bar.next()
        bar.finish()

        codon_efficiency_comparison = pd.DataFrame({"triplett": [triplett for triplett in codon_efficiencies],
                                                    "mean":     [None for _ in codon_efficiencies],
                                                    "median":   [None for _ in codon_efficiencies],
                                                    "size":     [None for _ in codon_efficiencies],
                                                    "ks":       [None for _ in codon_efficiencies],
                                                    "mw":       [None for _ in codon_efficiencies]},
                                                    index=[triplett for triplett in codon_efficiencies])
        
        for triplett1 in codon_efficiencies:
            preds1 = codon_efficiencies[triplett1]
            preds2 = []
            [preds2.extend(codon_efficiencies[triplett2]) for triplett2 in codon_efficiencies if triplett2 != triplett1]
            codon_efficiency_comparison.at[triplett1, "mean"]   = np.mean(preds1)
            codon_efficiency_comparison.at[triplett1, "median"] = np.median(preds1)
            codon_efficiency_comparison.at[triplett1, "size"]   = len(preds1)
            codon_efficiency_comparison.at[triplett1, "ks"] = scipy.stats.kstest(preds1, preds2).pvalue
            codon_efficiency_comparison.at[triplett1, "mw"] = scipy.stats.mannwhitneyu(preds1, preds2).pvalue

        # calculate adjusted pvalues
        codon_efficiency_comparison["ks"] = sm.stats.fdrcorrection(codon_efficiency_comparison["ks"], alpha=0.05)[1]
        codon_efficiency_comparison["mw"] = sm.stats.fdrcorrection(codon_efficiency_comparison["mw"], alpha=0.05)[1]
        return codon_efficiency_comparison.sort_index()
    

    # function (<-) added on 250528
    def analyze_stop_codons(self):
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)

        # conduct calculation for redundant variants
        codon_efficiency_comparison = self._analyze_stop_codons()
        print(codon_efficiency_comparison)
        codon_efficiency_comparison.to_csv(self.newdir+self.params["os_sep"]+"codon_efficiency_comparison.txt", sep=",", index=False)

        # repeat calculation with non-redundant variants
        self.variants               = self.variants.drop_duplicates(subset="ID:transcript id")
        codon_efficiency_comparison = self._analyze_stop_codons()
        print(codon_efficiency_comparison)
        codon_efficiency_comparison.to_csv(self.newdir+self.params["os_sep"]+"codon_efficiency_comparison_non_redundant.txt", sep=",", index=False)

        # calculate idealized distributions
        bases          = ["A", "C", "G", "T"]
        tripletts = []
        for base1 in bases:
            for base2 in bases:
                for base3 in bases:
                    tripletts.append(base1+base2+base3)

        codon_distribution      = {triplett:       [0 for _ in range(self.params["transformation_steps"])] for triplett in tripletts}
        stop_codon_distribution = {"nonsense":     [0 for _ in range(self.params["transformation_steps"])],
                                   "frameshift_1": [0 for _ in range(self.params["transformation_steps"])],
                                   "frameshift_2": [0 for _ in range(self.params["transformation_steps"])],
                                   "control":      [0 for _ in range(self.params["transformation_steps"])]}

        bar = IncrementalBar(set_bar("analyzing stop codons"), max=self.variants.shape[0])
        for i in range(self.variants.shape[0]):
            cds = self.variants.iloc[i].loc["FEATURE:all cds"]
            
            # track stop codon distribution
            for j in range(3, len(cds)-3, 3):
                # calculate relative position in shared virtual coordinate system
                pos = int((j-3)/(len(cds)-6)*self.params["transformation_steps"])
                stop_codon_distribution["control"][pos] += 1

                # detect nonsense candidates
                if (cds[j:j+2] == "TA" or cds[j:j+2] == "TG"
                    or cds[j+1:j+3] == "AA" or cds[j+1:j+3] == "AG" or cds[j+1:j+3] == "GA"
                    or (cds[j] == "T" and cds[j+2] == "A") or (cds[j] == "T" and cds[j+2] == "G")):
                    stop_codon_distribution["nonsense"][pos] += 1

                # detect candidates for frameshift +1
                if j+4 < len(cds) and cds[j+1:j+4] in ["TAA", "TAG", "TGA"]:
                    stop_codon_distribution["frameshift_1"][pos] += 1

                # detect candidates for frameshift +2
                if j+5 < len(cds) and cds[j+2:j+5] in ["TAA", "TAG", "TGA"]:
                    stop_codon_distribution["frameshift_2"][pos] += 1

            # track codon distribution
            for j in range(2, len(cds)-4, 1):
                pos = int((j-2)/(len(cds)-6)*self.params["transformation_steps"])

                if cds[j:j+3] in codon_distribution:
                    codon_distribution[cds[j:j+3]][pos] += 1

            bar.next()
        bar.finish()

        # plot and store stop codon distribution
        stop_codon_distribution = pd.DataFrame(stop_codon_distribution)
        plt.plot(stop_codon_distribution["nonsense"], color="red", label="nonsense")
        plt.plot(stop_codon_distribution["frameshift_1"], color="blue", label="frameshift_1")
        plt.plot(stop_codon_distribution["frameshift_2"], color="magenta", label="frameshift_2")
        plt.plot(stop_codon_distribution["control"], color="black", label="control")
        plt.legend()
        plt.show()
        stop_codon_distribution.to_csv(self.newdir+self.params["os_sep"]+"stop_codon_distribution.txt", sep=",", index=False)

        # plot and store codon distribution
        codons                    = sorted(list(codon_distribution.keys()))
        sorted_codon_distribution = {}
        for codon in codons:
            sorted_codon_distribution[codon] = codon_distribution[codon]
            plt.plot(np.arange(100), codon_distribution[codon], label=codon)
        
        plt.legend()
        plt.show()

        codon_distribution = pd.DataFrame({codon: sorted_codon_distribution[codon] for codon in sorted_codon_distribution})
        codon_distribution.to_csv(self.newdir+self.params["os_sep"]+"codon_distribution.txt", sep=",", index=False)


    def analyze_total_selection(self):
        model_predictions = []
        [model_predictions.extend(self.predictions.iloc[i].loc["predictions"]) for i in range(self.predictions.shape[0])]
        variant_predictions = []
        [variant_predictions.extend(self.predictions.iloc[i].loc["real predictions"]) for i in range(self.predictions.shape[0])]

        # marked (<-) added / removed on 250603
        # test_results = self.conduct_test(model_predictions, variant_predictions, self.params["variant_test_targets"][0]) # <- removed
        if self.params["apply_mutation_stats"] == False: # <- added
            # marked (<-) added on 250821 to allow weights for masked calculations
            if self.params["mutation_stats_ptc_weights"] == True or self.params["apply_mutation_weights"] == True:
                model_probabilities = []
                [model_probabilities.extend(self.predictions.iloc[i].loc["probabilities"]) for i in range(self.predictions.shape[0])]
                test_results = self.conduct_test_w_mutation_stats(model_predictions, model_probabilities, variant_predictions, self.params["variant_test_targets"][0]) # <- added

            else:
                test_results = self.conduct_test(model_predictions, variant_predictions, self.params["variant_test_targets"][0])

        if self.params["apply_mutation_stats"] == True:
            model_probabilities = []
            [model_probabilities.extend(self.predictions.iloc[i].loc["probabilities"]) for i in range(self.predictions.shape[0])]
            test_results = self.conduct_test_w_mutation_stats(model_predictions, model_probabilities, variant_predictions, self.params["variant_test_targets"][0]) # <- added

        total_predictions                                = pd.DataFrame({col: [None] for col in self.predictions.columns})
        total_predictions.at[0, "block id"]              = "total"
        # avg. real size must be present in order to avoid elimination before storage
        total_predictions.at[0, "avg. model prediction"] = np.nanmean(model_predictions)
        total_predictions.at[0, "avg. real prediction"]  = np.nanmean(variant_predictions)

        # marked (<-) added on 2502523
        total_predictions.at[0, "median model prediction"] = np.nanmedian(model_predictions)
        total_predictions.at[0, "median real prediction"]  = np.nanmedian(variant_predictions)
        total_predictions.at[0, "avg. real size"]          = len(variant_predictions)
        # "predictions" must be represented for later processing / empty to save storage
        total_predictions.at[0, "predictions"]           = []
        total_predictions.at[0, "real predictions"]      = variant_predictions

        for key in test_results:
            total_predictions.at[0, key] = test_results[key]
        
        self.predictions = pd.concat([self.predictions, total_predictions])


    def apply_variant_filter(self, data, target_col1, target_col2, mode="variant_filter"):
        if mode == "variant_filter": variant_filter = self.variant_filter
        if mode == "variants":       variant_filter = self.variants

        data[target_col1]           = [data.iloc[i].loc[target_col1].split(".")[0] for i in range(data.shape[0])]
        variant_filter[target_col2] = [variant_filter.iloc[i].loc[target_col2].split(".")[0] for i in range(variant_filter.shape[0])]
        init_shape                  = data.shape[0]

        if target_col1 in data and target_col2 in variant_filter:
            data = data[data[target_col1].isin(variant_filter[target_col2])]

        else:
            print("< variant filter column", self.params["variant_filter_col"], "was not found.")
            exit()

        print("< applying variant filter reduced data from", init_shape, "to", data.shape[0])
        return data
    

    # <- function added on 250605
    def calculate_distribution(self, values, probabilities):
        if len(values) != len(probabilities):
            print("< dimension error occurred @calculate_distribution ("+str(len(values)) + "/" + str(len(probabilities)) + ")")
            exit()

        updated_values = []

        if len(probabilities) > 0:
            max_value      = np.max(probabilities)
            random.seed(int(time.time()-seed_time))

            if self.params["bootstrapping_steps"] == None:
                self.params["bootstrapping_steps"] = len(values)

            last_index = 0
            #cutoffs = []; cutoffs2 = []
            bar = IncrementalBar(set_bar("calculating distribution"), max=self.params["bootstrapping_steps"]/10000)
            while len(updated_values) < self.params["bootstrapping_steps"]:
                rand_index      = int(random.random()*len(values))
                decision_cutoff = random.random()

                if probabilities[rand_index]/max_value >= decision_cutoff:
                    updated_values.append(values[rand_index])

                if len(updated_values) % 10000 == 0 and len(updated_values) != last_index:
                    last_index = len(updated_values)
                    bar.next()

            bar.finish()

        return updated_values


    def conduct_test(self, model_predictions, variant_predictions, target_col):
        # marked (<-) added on 250528
        if len(self.params["masks"]) > 0: # <- added on 250528
            model_predictions = [pred for pred in model_predictions if pd.isna(pred) == False] # <- added on 250528

        test_results = {}
        if len(model_predictions) > 0: # added on 250529
            pvalue, _, _ = bootstrapping_test(model_predictions, np.mean(variant_predictions), len(variant_predictions), simulation_steps=100) # added on 250529
            test_results["bootstrapping-pvalue "+target_col] = pvalue # added on 250529

        model_median         = np.median(model_predictions) # <- added from here
        model_count_below    = np.where(np.array(model_predictions) < model_median)[0].shape[0]
        model_count_above    = np.where(np.array(model_predictions) > model_median)[0].shape[0]
        real_count_below     = np.where(np.array(variant_predictions) < model_median)[0].shape[0]
        real_count_above     = np.where(np.array(variant_predictions) > model_median)[0].shape[0]

        expected_probability = None; observed_probability = None
        if model_count_below+model_count_above > 0:
            expected_probability = model_count_below / (model_count_below+model_count_above)

        if real_count_below+real_count_above > 0:
            observed_probability = real_count_below / (real_count_below+real_count_above)            

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False and observed_probability >= expected_probability:
            # test_results["binomial-pvalue "+target_col] = 1-stats.binom.cdf(real_count_below, (real_count_below+real_count_above), expected_probability) # <- removed on 250701
            test_results["binomial-pvalue "+target_col] = 1-stats.binom.cdf(real_count_below-1, (real_count_below+real_count_above), expected_probability) # <- added on 250701

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False and observed_probability < expected_probability:
            test_results["binomial-pvalue "+target_col] = stats.binom.cdf(real_count_below, (real_count_below+real_count_above), expected_probability)

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False:
            test_results["binomial-statistic "+target_col] = observed_probability-expected_probability # <- until here
        
        if "binomial-pvalue "+target_col in test_results:
            print("<", model_median, expected_probability, "/", len(model_predictions), real_count_below, real_count_above,
                test_results["binomial-pvalue "+target_col], test_results["binomial-statistic "+target_col])

        # Fisher's exact test on section-wise deviations
        model_mean   = np.mean(model_predictions)
        variant_mean = np.mean(variant_predictions)

        model_count = np.where(np.array(model_predictions) < model_mean)[0].shape[0]
        real_count  = np.where(np.array(variant_predictions) < model_mean)[0].shape[0]

        if len(variant_predictions) > 0 and len(model_predictions) > 0:
            fishers_exact = stats.fisher_exact([[real_count, len(variant_predictions)-real_count], [model_count, len(model_predictions)-model_count]])
            test_results["fishers exact-statistic "+target_col] = fishers_exact.statistic
            test_results["fishers exact-pvalue "+target_col]    = fishers_exact.pvalue

        # Fisher's exact test on cutoff-wise deviations: escape
        model_count   = np.where(np.array(model_predictions) <= self.params["nmd_escape_cutoff"])[0].shape[0]
        real_count    = np.where(np.array(variant_predictions) <= self.params["nmd_escape_cutoff"])[0].shape[0]
        fishers_exact = stats.fisher_exact([[real_count, len(variant_predictions)-real_count], [model_count, len(model_predictions)-model_count]])

        if len(variant_predictions) > 0 and len(model_predictions) > 0:
            test_results["fishers exact escape-statistic "+target_col] = fishers_exact.statistic
            test_results["fishers exact escape-pvalue "+target_col]    = fishers_exact.pvalue

        # Fisher's exact test on cutoff-wise deviations: target
        model_count   = np.where(np.array(model_predictions) >= self.params["nmd_target_cutoff"])[0].shape[0]
        real_count    = np.where(np.array(variant_predictions) >= self.params["nmd_target_cutoff"])[0].shape[0]

        fishers_exact = stats.fisher_exact([[real_count, len(variant_predictions)-real_count], [model_count, len(model_predictions)-model_count]])

        if len(variant_predictions) > 0 and len(model_predictions) > 0:
            test_results["fishers exact target-statistic "+target_col] = fishers_exact.statistic
            test_results["fishers exact target-pvalue "+target_col]    = fishers_exact.pvalue

        return test_results


    # <- function added on 250603
    def conduct_test_w_mutation_stats(self, model_predictions, model_probabilities, variant_predictions, target_col):
        test_results        = {}

        model_predictions   = [pred for pred in model_predictions if pd.isna(pred) == False]
        model_probabilities = [model_prob for model_prob in model_probabilities if pd.isna(model_prob) == False]

        model_median        = np.median(model_predictions)
        model_prob_below    = np.sum([model_probabilities[i] for i in range(len(model_probabilities)) if model_predictions[i] < model_median])
        model_prob_above    = np.sum([model_probabilities[i] for i in range(len(model_probabilities)) if model_predictions[i] > model_median])
        real_count_below    = np.where(np.array(variant_predictions) < model_median)[0].shape[0]
        real_count_above    = np.where(np.array(variant_predictions) > model_median)[0].shape[0]

        expected_probability = None; observed_probability = None
        if model_prob_below+model_prob_above > 0:
            expected_probability = model_prob_below / (model_prob_below+model_prob_above)

        if real_count_below+real_count_above > 0:
            observed_probability = real_count_below / (real_count_below+real_count_above)            

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False and observed_probability >= expected_probability:
            test_results["binomial-pvalue "+target_col] = 1-stats.binom.cdf(real_count_below-1, (real_count_below+real_count_above), expected_probability)

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False and observed_probability < expected_probability:
            test_results["binomial-pvalue "+target_col] = stats.binom.cdf(real_count_below, (real_count_below+real_count_above), expected_probability)

        if pd.isna(expected_probability) == False and pd.isna(observed_probability) == False:
            test_results["binomial-statistic "+target_col] = observed_probability-expected_probability

        return test_results


    def correct_predictions_by_exon(self, predictions_by_exon):
        predictions_by_exon = predictions_by_exon.replace("{", "{\"").replace("], ", "], \"").replace(":", "\":")
        return predictions_by_exon
    

    def _create_appris_selection(self, blocks, stats, temp):
        block_predictions = self.predictions.iloc[temp]

        if block_predictions.shape[0] != block_predictions.drop_duplicates(subset="transcript id").shape[0]:
            print("< duplicates detected @", block_predictions["gene id"].tolist(), block_predictions["transcript id"].tolist())
        
        if len(self.params["appris_priorities"]) > 0:
            appris_priorities = {key: [] for key in self.params["appris_priorities"]}
            [appris_priorities[block_predictions.iloc[i].loc["appris annotation"]].append(i) for i in range(block_predictions.shape[0])
             if block_predictions.iloc[i].loc["appris annotation"] in appris_priorities]

            temp = []
            i = 0
            while i < len(list(appris_priorities.keys())) and len(temp) == 0:
                if len(appris_priorities[list(appris_priorities.keys())[i]]) > 0:
                    temp                                      = appris_priorities[list(appris_priorities.keys())[i]]
                    stats[list(appris_priorities.keys())[i]] += 1

                i += 1

            # if multiple entries are available, only first one is selected
            if len(temp) > 1: temp = [temp[0]]

            block_predictions = block_predictions.iloc[temp]
        
        blocks.append(block_predictions)
        return blocks, stats


    def create_appris_selection(self):
        self.predictions = self.predictions[self.predictions["appris annotation"] != "-"]
        self.predictions = self.predictions.sort_values(by="gene id")

        init_shape       = self.predictions.shape[0]
        blocks   = []
        stats    = {key: 0 for key in self.params["appris_priorities"]}
        last_ids = {"gene id": None}
        temp     = []
        bar      = IncrementalBar(set_bar("creating appris selection"), max=self.predictions.shape[0])

        for i in range(self.predictions.shape[0]):
            check = [block_target for block_target in last_ids if last_ids[block_target] == self.predictions.iloc[i].loc[block_target]]

            if len(temp) == 0 or len(check) == len(last_ids):
                temp.append(i)

            elif len(temp) > 0:
                blocks, stats = self._create_appris_selection(blocks, stats, temp)
                temp   = [i]

            last_ids = {block_target: self.predictions.iloc[i].loc[block_target] for block_target in last_ids}

            bar.next()
        bar.finish()

        # create last block
        if len(temp) > 0: blocks, stats = self._create_appris_selection(blocks, stats, temp)

        self.predictions = pd.concat(blocks)
        print("< appris stats")
        print(json.dumps(stats, indent=4))
        appris_sum      = np.sum(stats[key] for key in stats)
        if appris_sum != self.predictions.shape[0]: print("< inconsistent sizes @apply_filter", appris_sum, "/", self.predictions.shape[0])
        duplicates_test = self.predictions[self.predictions.duplicated(subset="gene id")]

        print("< appris selection reduced entries from", init_shape, "to", self.predictions.shape[0], ". duplicates: ", duplicates_test.shape[0]) # new version
    

    def convert_lindeboom_predictions(self, lines):
        lindeboom_dictionary = {}
        cds_exists           = {}

        exception = "numeric" # exception found in Lindeboom preditions for hg19
        bar   = IncrementalBar(set_bar("creating Lindeboom dictionary"), max=len(lines)/1000)
        last_line = ""; steps = 0

        for line in lines:
            if line.split("\t")[2] == "exon" or line.split("\t")[2] == "CDS":
                residual = line.split("\t")
                chr      = residual[0]
                cdsstart = residual[3]
                cdsend   = residual[4]
                strand   = residual[6]
                gene_id  = residual[8].split(" ")[1].strip(";")
                gene_id  = gene_id.split(".")[0] # remove version number

            if line.split("\t")[2] == "exon":
                if strand != "-" and strand != "+":
                    print("strand error occurred @line:")
                    print(line)

                else:
                    if gene_id not in lindeboom_dictionary.keys():  lindeboom_dictionary[gene_id] = [{"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                                      "strand": strand, "scores": []}]
                    else:                                           lindeboom_dictionary[gene_id].append({"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                                          "strand": strand, "scores": []})

            if line.split("\t")[2] == "CDS":
                # check if previous line contains exon (as expected)
                if last_line.split("\t")[2] != "exon":
                    print("< format error @convert_lindeboom_predictions")
                    print(last_line)
                    exit()

                scores   = residual[8].split(" ")[9].strip(";\n").split(",")
                scores   = [float(score) for score in scores if exception not in score]

                if strand != "-" and strand != "+":
                    print("strand error occurred @line:")
                    print(line)

                else:
                    cds_exists[gene_id] = True

                    if gene_id not in lindeboom_dictionary.keys():  lindeboom_dictionary[gene_id] = [{"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                                      "strand": strand, "scores": scores}]
                    else:                                           lindeboom_dictionary[gene_id].append({"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                                          "strand": strand, "scores": scores})

            last_line = line
            if steps % 1000 == 0: bar.next()
            steps += 1

        bar.finish()

        updated_lindeboom_dictionary = {}
        # re-iterate over gene dictionary to remove excessive keys
        bar = IncrementalBar(set_bar("removing excessive keys"), max=len(lindeboom_dictionary))
        for gene_id in lindeboom_dictionary:
            if gene_id in cds_exists:
                # remove exons directly before cds (identical exon)
                selected_index                        = [i for i in range(len(lindeboom_dictionary[gene_id]))
                                                         if len(lindeboom_dictionary[gene_id][i]["scores"]) > 0
                                                         or (i < len(lindeboom_dictionary[gene_id])-1
                                                             and (len(lindeboom_dictionary[gene_id][i]["scores"]) == 0 and len(lindeboom_dictionary[gene_id][i+1]["scores"]) == 0))
                                                         or i == len(lindeboom_dictionary[gene_id])-1] 

                updated_lindeboom_dictionary[gene_id] = [lindeboom_dictionary[gene_id][i] for i in selected_index]

            bar.next()
        bar.finish()

        lindeboom_dictionary = updated_lindeboom_dictionary

        # create output identical to prediction scores
        predictions = {"gene id": [], "transcript id": [], "uc id": [], "chr": [], "strand": [], "predictions": [], "predictions_by_exon": [], "fails": []}
        bar = IncrementalBar(set_bar("rearranging lindeboom predictions"), max=len(lindeboom_dictionary))

        for id in lindeboom_dictionary:
            all_scores  = []
            exon_scores = {}

            for i in range(len(lindeboom_dictionary[id])):
                if lindeboom_dictionary[id][0]["strand"] == "+":
                    exon_scores[lindeboom_dictionary[id][i]["cdsstart"]] = convert_labels(lindeboom_dictionary[id][i]["scores"])
                    all_scores.extend(convert_labels(lindeboom_dictionary[id][i]["scores"]))

                if lindeboom_dictionary[id][0]["strand"] == "-":
                    exon_scores[lindeboom_dictionary[id][i]["cdsstart"]] = convert_labels(invert(lindeboom_dictionary[id][i]["scores"]))
                    if len(exon_scores[lindeboom_dictionary[id][i]["cdsstart"]]) > 0 and exon_scores[lindeboom_dictionary[id][i]["cdsstart"]][0] != convert_labels(lindeboom_dictionary[id][i]["scores"])[-1]:
                        print("< conversion error1 @", id)
                        exit()

            # reverse order of scores for -strand
            if lindeboom_dictionary[id][0]["strand"] == "-":
                rev_exon_scores = {}
                keys            = list(exon_scores.keys())

                for i in range(len(keys)):
                    rev_exon_scores[keys[len(keys)-1-i]] = exon_scores[keys[len(keys)-1-i]]
                    all_scores.extend(exon_scores[keys[len(keys)-1-i]])
                
                if len(all_scores) > 0 and len(exon_scores[list(exon_scores.keys())[-1]]) > 0 and all_scores[0] != exon_scores[list(exon_scores.keys())[-1]][0]:
                    print("< conversion error2 @", id)
                    exit()

                if len(all_scores) > 0 and len(exon_scores[list(exon_scores.keys())[0]]) > 0 and all_scores[-1] != exon_scores[list(exon_scores.keys())[0]][-1]:
                    print("< conversion error3 @", id)
                    exit()

                exon_scores = rev_exon_scores

            predictions["gene id"].append(None)
            predictions["transcript id"].append(None)
            predictions["uc id"].append(id)
            predictions["chr"].append(lindeboom_dictionary[id][0]["chr"])
            predictions["strand"].append(lindeboom_dictionary[id][0]["strand"])
            predictions["predictions"].append(all_scores)
            predictions["predictions_by_exon"].append(exon_scores)
            predictions["fails"].append([])
            bar.next()
        bar.finish()


        # load knownGene dictionary
        knowngene                        = pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["knowngene_fname"], delimiter=",", index_col=False).sort_values(by=["chr", "cdsstart"])
        
        predictions                      = pd.DataFrame(predictions)
        predictions["uc id"]             = [predictions.iloc[i].loc["uc id"].split(".")[0] for i in range(predictions.shape[0])]
        knowngene["uc id"]               = [knowngene.iloc[i].loc["uc id"].split(".")[0] for i in range(knowngene.shape[0])]
        knowngene["index"]               = [i for i in range(knowngene.shape[0])]
        selected_index                   = append_df_with_mapping([predictions, knowngene], "uc id", "uc id", "index", text="mapping uc ids", reverse=True, show_progress=True)
        predictions["gene id"]           = [knowngene.iloc[int(i)].loc["gene id"] if i != "-" else "-" for i in selected_index]
        predictions["transcript id"]     = [knowngene.iloc[int(i)].loc["transcript id"].split(".")[0] if i != "-" else "-" for i in selected_index]
        predictions["appris annotation"] = [knowngene.iloc[int(i)].loc["appris annotation"] if i != "-" else "-" for i in selected_index]

        predictions.to_csv(path_or_buf=self.params["data_dir"]+self.params["os_sep"]+"h38_lindeboom_predictions.txt", sep=",", index=False)
        return
    

    def create_newdir(self):
        datestr = str(datetime.now())
        index   = datestr.index(".")
        dirname = datestr[:index].replace(" ", "_").replace(":", "-")
        newdir  = self.params["data_dir"] + self.params["os_sep"] + dirname + "_" + self.params["file_tag"]
        print("< new directory:", newdir)
        if not os.path.exists(newdir): os.mkdir(newdir)
        return newdir


    def filter_errors(self):
        bar            = IncrementalBar(set_bar("filtering errors"), max=self.predictions.shape[0])
        filtered_index = []

        for i in range(self.predictions.shape[0]):
            fails  = json.dumps(self.predictions.iloc[i].loc["fails"])
            errors = [fails[key] for key in fails if key in self.params["error_filter"]]
       
            if np.sum(errors) == 0: filtered_index.append(i)
            bar.next()
        bar.finish()

        init_shape       = self.predictions.shape[0]
        self.predictions = self.predictions.iloc[filtered_index]
        self.predictions = self.predictions.reset_index()
        print("< error filtering reduced data from", init_shape, "to", self.predictions.shape[0])


    # <- modified on 251114 to recognize both capital and non-capital letters
    def get_biallelic_only(self, data):
        init_shape = data.shape[0]
        data       = data[data["ID:chromosome"] != "chrY"]
        data       = data[[True if data.iloc[i].loc["ID:chromosome"] != "chrX" or (data.iloc[i].loc["ID:chromosome"] == "chrX"
                           and type(data.iloc[i].loc["ID:gender"]) == str and data.iloc[i].loc["ID:gender"].upper() == "FEMALE") else False for i in range(data.shape[0])]]
        print("< biallelic filtering reduced dataset from", init_shape, "to", data.shape[0])
        return data


    # count cases in status file that were included in PTC extraction (containing both WXS and RNA data)
    def get_cases(self, status, outer_block_target):
        cases = 0
        for case in status["file_ids"][outer_block_target]:
            if len(status["file_ids"][outer_block_target][case]["RNA"]) > 0 and len(status["file_ids"][outer_block_target][case]["WXS"]) > 0:
                cases += 1

        return cases


    def interpolate(self, transformed_predictions, predictions, it, positions=None, transformed_targets=None, exon_category=""):
        max_position = len(predictions)
        if len(self.params["masks"]) > 0:
            transformed_positions = [float(i)/max_position for i in range(len(predictions)) if pd.isna(predictions[i]) == False]
            predictions           = [pred for pred in predictions if pd.isna(pred) == False]       

        else:
            transformed_positions = [float(i)/max_position for i in range(len(predictions))]

        if len(transformed_positions) < self.params["projection_cutoff"]:
            print("< sample size too small @interpolate ("+str(len(transformed_positions))+"/"+str(self.params["projection_cutoff"])+")")
            exit()

        fit_successful = False
        try:
            fit_params     = np.polyfit(transformed_positions, predictions, deg=self.params["polynomial_degree"])

            # marked entries (<-) added/removed on 250408 to interpolate distributions using cubic spline rather than polynomial fit
            cs_fit = CubicSpline(transformed_positions, predictions) # <- added
            fit_successful = True # <- shifted from above

        except:
            print("< fit failed @", it, exon_category)
            
        if fit_successful == True:
            transformation_steps        = [float(i)/self.params["transformation_steps"] for i in range(self.params["transformation_steps"])]
            transformed_predictions_col = cs_fit(transformation_steps)
            if len(exon_category) == 0: [self.stats["projected_mean_hist"][i].append(transformed_predictions_col[i]) for i in range(len(transformed_predictions_col))]
            if len(exon_category) > 0:  [self.stats["projected_"+exon_category+"_exon_mean_hist"][i].append(transformed_predictions_col[i]) for i in range(len(transformed_predictions_col))]

            if transformed_targets != None:
                for variant_key in self.params["variant_targets"]:
                    if len(positions) == len(transformed_targets[variant_key]):
                        if len(exon_category) == 0:
                            [self.stats[variant_key+"_projected_mean_hist"][int(self.params["transformation_steps"]*positions[i]/max_position)].append(transformed_targets[variant_key][i])
                             for i in range(len(positions))]
                            
                        if len(exon_category) > 0:
                            if np.max(positions) >= max_position:
                                print("< dimension error occurred @interpolate. exon category", exon_category,
                                      "@", it, self.predictions.iloc[it].loc["gene id"], self.predictions.iloc[it].loc["transcript id"], "max. position", max_position)
                                print("  positions")
                                print(positions)

                            [self.stats[variant_key+"_projected_"+exon_category+"_exon_mean_hist"][int(self.params["transformation_steps"]*positions[i]/max_position)].append(transformed_targets[variant_key][i])
                             for i in range(len(positions))]

                    else:
                        print("< mismatching sizes @get_projection",  len(positions), "/", len(transformed_targets[variant_key]))
                        print(positions)
                        print(transformed_targets[variant_key])

            if len(exon_category) == 0: 
                if transformed_predictions.shape[0] > 0: transformed_predictions = np.column_stack((transformed_predictions, transformed_predictions_col))
                else:                                    transformed_predictions = np.array(transformed_predictions_col)

        return transformed_predictions, fit_successful
    

    def get_projection(self):
        projection_stats        = {"full": {"successful": 0, "not successful": 0}, "exons": {"successful": 0, "not successful": 0}}
        transformed_predictions = np.zeros((0))
        # prior to pca, dimensions must be unified
        bar = IncrementalBar(set_bar("calculating projections"), max=self.predictions.shape[0])
        for i in range(self.predictions.shape[0]):
            # get strand-corrected list of exons
            exons = list(self.predictions.iloc[i].loc["predictions_by_exon"].keys())
            
            # interpolate all predictions for cds
            if ((len(self.params["masks"]) == 0 and len(self.predictions.iloc[i].loc["predictions"]) >= self.params["projection_cutoff"])
                or (len(self.params["masks"]) > 0 and len([pred for pred in self.predictions.iloc[i].loc["predictions"] if pd.isna(pred) == False]) >= self.params["projection_cutoff"])): # <- added
                if pd.isna(self.predictions.iloc[i].loc["variant targets"]) == False:
                    positions = []
                    if "FEATURE:ptc cds position" in self.predictions.iloc[i].loc["variant targets"]:
                        [positions.extend(self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"][exon])
                         for exon in self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"]]
                
                    transformed_targets = {variant_key: [] for variant_key in self.params["variant_targets"]}
                    for variant_key in self.params["variant_targets"]:
                        [transformed_targets[variant_key].extend(self.predictions.iloc[i].loc["variant targets"][variant_key][exon])
                         for exon in self.predictions.iloc[i].loc["variant targets"][variant_key]]
                    
                    transformed_predictions, fit_successful = self.interpolate(transformed_predictions, self.predictions.iloc[i].loc["predictions"],
                                                                               i, positions=positions, transformed_targets=transformed_targets)

                else:
                    transformed_predictions, fit_successful = self.interpolate(transformed_predictions, self.predictions.iloc[i].loc["predictions"], i)

                if fit_successful == True:  projection_stats["full"]["successful"]     += 1
                if fit_successful == False: projection_stats["full"]["not successful"] += 1

                # interpolate predictions per exon
                preceeding_length = 0
                if self.params["exon_resolution"] == True:
                    for j in range(len(exons)):
                        if j == 0 and len(exons) > 1:                     exon_category = "first"
                        if j > 0 and j < len(exons)-1 and len(exons) > 1: exon_category = "middle"
                        if j == len(exons)-1:                             exon_category = "last"

                        # marked (<-) added / removed on 250528
                        if ((len(self.params["masks"]) == 0 and len(self.predictions.iloc[i].loc["predictions_by_exon"][exons[j]]) >= self.params["projection_cutoff"])
                            or (len(self.params["masks"]) > 0 and len([pred for pred in self.predictions.iloc[i].loc["predictions_by_exon"][exons[j]] if pd.isna(pred) == False]) >= self.params["projection_cutoff"])): # <- added

                            if (pd.isna(self.predictions.iloc[i].loc["variant targets"]) == False
                                and "FEATURE:ptc cds position" in self.predictions.iloc[i].loc["variant targets"] and len(self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"][exons[j]]) > 0):
                                transformed_targets = {variant_key: [] for variant_key in self.params["variant_targets"]}

                                for variant_key in self.params["variant_targets"]:
                                    transformed_targets[variant_key] = self.predictions.iloc[i].loc["variant targets"][variant_key][exons[j]]
                                
                                positions = [value-preceeding_length for value in self.predictions.iloc[i].loc["variant targets"]["FEATURE:ptc cds position"][exons[j]]]

                                if np.min(positions) < 0 or np.max(positions) >= len(self.predictions.iloc[i].loc["predictions"]):
                                    print("< dimension error occurred @interpolate. min:",  np.min(positions), "max:", np.max(positions),
                                          "full size:", len(self.predictions.iloc[i].loc["predictions"]))

                                else:
                                    transformed_predictions, fit_successful = self.interpolate(transformed_predictions, self.predictions.iloc[i].loc["predictions_by_exon"][exons[j]], i,
                                                                                               positions=positions, transformed_targets=transformed_targets, exon_category=exon_category)
                                
                            else:
                                transformed_predictions, fit_successful = self.interpolate(transformed_predictions, self.predictions.iloc[i].loc["predictions_by_exon"][exons[j]], i,
                                                                                           exon_category=exon_category)

                            if fit_successful == True:  projection_stats["exons"]["successful"]   += 1
                            if fit_successful == False: projection_stats["exons"]["not successful"] += 1

                        preceeding_length += len(self.predictions.iloc[i].loc["predictions_by_exon"][exons[j]])
            
            bar.next()
        bar.finish()
        print(json.dumps(projection_stats, indent=4))
    

    def get_polynomial_(self, x, fit_params):
        return np.sum([fit_params[i]*math.pow(x, len(fit_params)-1-i) for i in range(len(fit_params))])


    def get_polynomial(self, x, fit_params):
        return [self.get_polynomial_(x[i], fit_params) for i in range(len(x))]


    def get_prediction_stats(self, it):
        if len(self.predictions.iloc[it].loc["predictions"]) > 0:
            exons = list(self.predictions.iloc[it].loc["predictions_by_exon"].keys())

            if self.params["variants_only"] == False:
                # calculate means (probably self.means is unused)
                self.means["gene id"].append(self.predictions.iloc[it].loc["gene id"])
                self.means["mean"].append(np.mean(self.predictions.iloc[it].loc["predictions"]))

                self.stats["mean"].append(np.mean(self.predictions.iloc[it].loc["predictions"]))
                self.stats["size"].append(len(self.predictions.iloc[it].loc["predictions"]))

                if self.params["exon_resolution"] == True:
                    if len(exons) > 1:
                        self.stats["first_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["predictions_by_exon"][exons[0]]))
                        self.stats["first_exon_size"].append(len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[0]]))

                        [self.stats["middle_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]])) for i in range(1, len(exons)-1)
                        if len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]]) > 0]
                        [self.stats["middle_exon_size"].append(len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]])) for i in range(1, len(exons)-1)
                        if len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]]) > 0]

                    self.stats["last_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["predictions_by_exon"][exons[-1]]))
                    self.stats["last_exon_size"].append(len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[-1]]))

                if self.params["stats_mode"] == "values":
                    self.stats["probabilities"].extend(self.predictions.iloc[it].loc["probabilities"]) # <- added on 250605
                    self.stats["values"].extend(self.predictions.iloc[it].loc["predictions"])

                    if self.params["exon_resolution"] == True:
                        if len(exons) > 1:
                            self.stats["first_exon_values"].extend(self.predictions.iloc[it].loc["predictions_by_exon"][exons[0]])
                            [self.stats["middle_exon_values"].extend(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]]) for i in range(1, len(exons)-1)
                            if len(self.predictions.iloc[it].loc["predictions_by_exon"][exons[i]]) > 0]

                            if type(self.predictions.iloc[it].loc["probabilities_by_exon"]) == dict: # <- added on 250607
                                self.stats["first_exon_probabilities"].extend(self.predictions.iloc[it].loc["probabilities_by_exon"][exons[0]]) # <- added on 250607
                                [self.stats["middle_exon_probabilities"].extend(self.predictions.iloc[it].loc["probabilities_by_exon"][exons[i]]) for i in range(1, len(exons)-1)
                                 if len(self.predictions.iloc[it].loc["probabilities_by_exon"][exons[i]]) > 0] # <- added on 250607

                        self.stats["last_exon_values"].extend(self.predictions.iloc[it].loc["predictions_by_exon"][exons[-1]])

                        if type(self.predictions.iloc[it].loc["probabilities_by_exon"]) == dict: # <- added on 250607
                            self.stats["last_exon_probabilities"].extend(self.predictions.iloc[it].loc["probabilities_by_exon"][exons[-1]]) # <- added on 250607

            # unclear whether None needs to be replaced by np.isnan
            if pd.isna(self.predictions.iloc[it].loc["variant targets"]) == False:
                for variant_key in self.params["variant_targets"]:
                    values = []
                    [[values.append(i) for i in self.predictions.iloc[it].loc["variant targets"][variant_key][exon]] for exon in exons]
                    
                    if len(values) > 0:
                        self.stats[variant_key + "_mean"].append(np.mean(values))
                        self.stats[variant_key + "_size"].append(len(values))
                        self.stats[variant_key + "_values"].extend(values)

                    if self.params["exon_resolution"] == True:
                        if len(exons) > 1 and len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[0]]) > 0:
                            self.stats[variant_key + "_first_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[0]]))
                            self.stats[variant_key + "_first_exon_size"].append(len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[0]]))
                            self.stats[variant_key + "_first_exon_values"].extend(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[0]])
                        
                        if len(exons) > 1:
                            [self.stats[variant_key + "_middle_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]])) for i in range(1, len(exons)-1)
                            if len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]]) > 0]
                            [self.stats[variant_key + "_middle_exon_size"].append(len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]])) for i in range(1, len(exons)-1)
                            if len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]]) > 0]
                            [self.stats[variant_key + "_middle_exon_values"].extend(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]]) for i in range(1, len(exons)-1)
                            if len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[i]]) > 0]
                                
                        if len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[-1]]) > 0:
                            self.stats[variant_key + "_last_exon_mean"].append(np.mean(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[-1]]))
                            self.stats[variant_key + "_last_exon_size"].append(len(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[-1]]))
                            self.stats[variant_key + "_last_exon_values"].extend(self.predictions.iloc[it].loc["variant targets"][variant_key][exons[-1]])
        
        
    def init_stats(self):
        if self.params["exon_resolution"] == False:
            self.base_dict = {
                              "mean"                             : [],
                              "probabilities"                    : [],
                              "size"                             : [],
                              "values"                           : [],
                              "projected_mean_hist"              : [[] for _ in range(self.params["transformation_steps"])]
                             }
            
        if self.params["exon_resolution"] == True:
            self.base_dict = {
                              "mean"                             : [],
                              "first_exon_mean"                  : [],
                              "last_exon_mean"                   : [],
                              "middle_exon_mean"                 : [],
                              "size"                             : [],
                              "first_exon_size"                  : [],
                              "last_exon_size"                   : [],
                              "middle_exon_size"                 : [],
                              "probabilities"                    : [],
                              "values"                           : [],
                              "first_exon_probabilities"         : [],
                              "last_exon_probabilities"          : [],
                              "middle_exon_probabilities"        : [],
                              "first_exon_values"                : [],
                              "last_exon_values"                 : [],
                              "middle_exon_values"               : [],
                              "projected_mean_hist"              : [[] for _ in range(self.params["transformation_steps"])],
                              "projected_first_exon_mean_hist"   : [[] for _ in range(self.params["transformation_steps"])],
                              "projected_last_exon_mean_hist"    : [[] for _ in range(self.params["transformation_steps"])],
                              "projected_middle_exon_mean_hist"  : [[] for _ in range(self.params["transformation_steps"])]
                             }
        

        self.stats = dict(self.base_dict)

        for variant_key in self.params["variant_targets"]:
            for base_key in self.base_dict:
                if "projected" not in base_key:
                    self.stats[variant_key + "_" + base_key]                           = []
                    self.stats[variant_key + "_" + base_key.replace("mean", "values")] = []

                if "projected" in base_key:
                    self.stats[variant_key + "_" + base_key] = [[] for _ in range(self.params["transformation_steps"])]


    def load(self, fnames, delimiter=",", mode="predictions", loading_key=None):
        data = []

        if mode == "global_predictions":
            with open(fnames[0], "r") as file:
                global_predictions = json.load(file)["values"]
            
            self.global_predictions = [global_predictions[i] for i in np.random.randint(low=0, high=len(global_predictions), size=10000)]


        else:
            if loading_key == None:
                for fname in fnames:
                    if os.path.isfile(self.params["data_dir"]+self.params["os_sep"]+fname) == True:
                        data.append(pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+fname, delimiter=delimiter, index_col=False))

                    else:
                        print("<", self.params["data_dir"]+self.params["os_sep"]+fname, "not found.")

            if loading_key != None:
                for fname in fnames:
                    data.append(load_by_key(self.params["data_dir"]+self.params["os_sep"]+fname, loading_key))

            if mode == "predictions":
                if len(data) == 1: self.predictions = data[0]
                else:              self.predictions = pd.concat(data, ignore_index=True)
            
                self.predictions = self.predictions[self.predictions["transcript id"] != "-"]
                if len(self.params["error_filter"]) > 0: self.filter_errors()

                # exclude Y-chromosomes
                init_shape       = self.predictions.shape[0]
                self.predictions = self.predictions[self.predictions["chr"] != "chrY"]
                print("< filtering of Y-chromosomes reduced dataset from", init_shape, "to", self.predictions.shape[0])

                # convert prediction and predition by exon to numerical list
                for i in range(self.predictions.shape[0]):
                    self.predictions.at[self.predictions.index[i], "predictions"] = json.loads(self.predictions.iloc[i].loc["predictions"])

                    if self.predictions.iloc[i].loc["predictions_by_exon"] != "{}":
                        corrected_predictions_by_exon                                         = self.correct_predictions_by_exon(self.predictions.iloc[i].loc["predictions_by_exon"])
                        corrected_predictions_by_exon                                         = json.loads(corrected_predictions_by_exon)
                        corrected_predictions_by_exon                                         = {key.replace("\'", ""): corrected_predictions_by_exon[key] for key in corrected_predictions_by_exon}
                        self.predictions.at[self.predictions.index[i], "predictions_by_exon"] = corrected_predictions_by_exon

                if self.params["appris_selection"] == True: self.create_appris_selection()
                self.predictions.index            = [self.predictions.iloc[i].loc["transcript id"].split(".")[0] for i in range(self.predictions.shape[0])]
                self.predictions["transcript id"] = [self.predictions.iloc[i].loc["transcript id"].split(".")[0] for i in range(self.predictions.shape[0])]

                if self.params["use_variant_filter"] == True:
                    self.predictions = self.apply_variant_filter(self.predictions, self.params["variant_filter_col"].replace("ID:", ""),
                                                                 self.params["variant_filter_col"], mode="variant_filter")
                    
                # <- added on 251111
                if ("apply_mutation_stats" in self.params and
                   (self.params["apply_mutation_stats"] == True or self.params["apply_mutation_weights"] == True or self.params["mutation_stats_ptc_weights"] == True)):
                    self.predictions = self.predictions[self.predictions["transcript id"].isin(self.mutation_stats["genes"]["total"])]
                    self.predictions = self.predictions[self.predictions["transcript id"].isin(self.seqs.index)]

            if mode == "variant_filter":
                if len(data) == 1: self.variant_filter = data[0]
                else:              self.variant_filter = pd.concat(data, ignore_index=True)

                if self.params["tcga_filter"] == "biallelic_only":
                    # filter out Y-chromosomes and X-chromosomes for men
                    self.variant_filter = self.get_biallelic_only(self.variant_filter)

                # now in shared_utils (250313)
                self.variant_filter = apply_value_filter(self.variant_filter, self.params["variant_filter_value_filter"], "variant filter")


            if mode == "variants":
                if len(data) == 1: self.variants = data[0]
                else:              self.variants = pd.concat(data, ignore_index=True)

                if self.params["use_variant_filter"] == True:
                    self.variants = self.apply_variant_filter(self.variants, self.params["variant_filter_col"], self.params["variant_filter_col"], mode="variant_filter")

                if "Cuomo" in self.params["id_filter"]:
                    init_shape    = self.variants.shape[0]
                    self.variants = self.variants.drop_duplicates(subset=["ID:cell id", "ID:variant id"])
                    print("< variants reduced from", init_shape, "to", self.variants.shape[0])

                if "MSK" in self.params["id_filter"] and self.params["tcga_filter"] == "biallelic_only":
                    # filter out Y-chromosomes and X-chromosomes for men
                    self.variants = self.get_biallelic_only(self.variants)

                if "Teran" in self.params["id_filter"]:
                    self.variants = self.variants.drop_duplicates(subset=["ID:subject id", "ID:variant id"])

                # <- modified on 251104
                if ("CPTAC3" in self.params["id_filter"] or "TCGA" in self.params["id_filter"]) and self.params["tcga_filter"] == "biallelic_only":
                    self.variants = self.get_biallelic_only(self.variants)

                if self.params["shared_variants_cases"] == True:
                    # performance is slow as 'duplicated' keeps all duplicated entries, second reduction step should accelarate (250321)
                    shared_variants          = self.variants[self.variants.duplicated(subset="ID:variant id")].drop_duplicates(subset="ID:variant id")["ID:variant id"].tolist() # new version
                    filtered_shared_variants = [shared_variant for shared_variant in shared_variants
                                                if self.variants[self.variants["ID:variant id"] == shared_variant].shape[0] >= self.params["shared_variants_filter"]]
                    shared_variants          = self.variants[self.variants["ID:variant id"].isin(filtered_shared_variants)]
                    shared_case_ids          = shared_variants.drop_duplicates(subset="ID:case id")["ID:case id"].tolist()
                    self.variants            = self.variants[self.variants["ID:case id"].isin(shared_case_ids)]

                if self.params["shared_variants"] == True:
                    # performance is slow as 'duplicated' keeps all duplicated entries, second reduction step should accelarate (250321)
                    shared_variants          = self.variants[self.variants.duplicated(subset="ID:variant id")].drop_duplicates(subset="ID:variant id")["ID:variant id"].tolist() # new version
                    filtered_shared_variants = [shared_variant for shared_variant in shared_variants
                                                if self.variants[self.variants["ID:variant id"] == shared_variant].shape[0] >= self.params["shared_variants_filter"]]
                    self.variants            = self.variants[self.variants["ID:variant id"].isin(filtered_shared_variants)]


                # filter out expression values below minimum (if specified)
                # now in shared_utils (250313)
                self.variants = apply_value_filter(self.variants, self.params["variant_value_filter"], "variants")
                
                # <- added on 251030
                if self.params["random_sample_size"] != None:
                    self.variants = self.variants.sample(frac=1).reset_index(drop=True)
                    self.variants = self.variants.iloc[0:self.params["random_sample_size"]]


                # report and remove missing values and check for expired placeholders ("-" and -1)
                for col in np.unique([*self.params["variant_targets"], *self.params["variant_test_targets"]]):
                    init_shape = self.variants.shape[0]

                    if col in self.variants.columns:
                        self.variants = self.variants[~self.variants[col].isna()]
                        print(init_shape-self.variants.shape[0], "values removed based on", col)

                    else:
                        print("< error.", col, "not found in variants.")
                        exit()
    
                    if self.variants[self.variants[col] == -1].shape[0] > 0:
                        print("< expired placeholder '-1' detected @", col)
                        exit()

                    if self.variants[self.variants[col] == "-"].shape[0] > 0:
                        print("< expired placeholder '-' detected @", col)
                        exit()

                        
    # mapping variant data to predictions and calculate significance of selection
    def map_variants(self):
        stats                  = {"first": 0, "middle": 0, "last": 0, "unassigned": 0}
        self.variants["index"] = [i for i in range(self.variants.shape[0])]

        # changed 241114, mapping to transcript id means perfect consistency of variants and full predictions with the disadvantage of loss of variants
        # mapping to gene id means no perfect consistency of variants and full predictions with the advantage of not losing variants
        mapped_index = append_df_with_mapping2([self.predictions, self.variants], self.params["mapping_target"], "ID:"+self.params["mapping_target"],
                                               "index", text=set_bar("mapping variants I"))

        self.predictions["avg. model prediction"]   = [None for _ in range(self.predictions.shape[0])]
        self.predictions["avg. real prediction"]    = [None for _ in range(self.predictions.shape[0])]
        # marked (<-) added on 250523
        self.predictions["median model prediction"] = [None for _ in range(self.predictions.shape[0])]
        self.predictions["median real prediction"]  = [None for _ in range(self.predictions.shape[0])]
        self.predictions["avg. real size"]          = [None for _ in range(self.predictions.shape[0])]
        self.predictions["block id"]                = [None for _ in range(self.predictions.shape[0])]
        self.predictions["probabilities"]           = [[] for _ in range(self.predictions.shape[0])]
        self.predictions["probabilities_by_exon"]   = [[] for _ in range(self.predictions.shape[0])]
        self.predictions["real predictions"]        = [[] for _ in range(self.predictions.shape[0])]
        self.predictions["variant index"]           = mapped_index
        self.predictions["variant targets"]         = [None for _ in range(self.predictions.shape[0])]

        for key in self.params["selection_cols"]:
            self.predictions[key] = [[] for _ in range(self.predictions.shape[0])]

        for variant_test_target in self.params["variant_test_targets"]:
            # marked entry (<-) added on 250408 to conduct binomial test using median as cutoff
            self.predictions["binomial-pvalue "+variant_test_target]                 = [None for _ in range(self.predictions.shape[0])]
            self.predictions["binomial-statistic "+variant_test_target]              = [None for _ in range(self.predictions.shape[0])]
            self.predictions["bootstrapping-pvalue "+variant_test_target]            = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact-statistic "+variant_test_target]         = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact-pvalue "+variant_test_target]            = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact escape-statistic "+variant_test_target]  = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact escape-pvalue "+variant_test_target]     = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact target-statistic "+variant_test_target]  = [None for _ in range(self.predictions.shape[0])]
            self.predictions["fishers exact target-pvalue "+variant_test_target]     = [None for _ in range(self.predictions.shape[0])]

        bar = IncrementalBar(set_bar("mapping variants II"), max=len(mapped_index))
        for i in range(len(mapped_index)):
            if len(mapped_index[i]) > 0:
                variant_targets = {key: {exon: [] for exon in self.predictions.iloc[i].loc["predictions_by_exon"]}
                                   for key in self.params["variant_targets"]}
                
                exons           = list(self.predictions.iloc[i].loc["predictions_by_exon"].keys())
                # reverse order if -strand for variant assignment to corresponding exon
                if len(exons) > 1 and int(exons[1]) < int(exons[0]):
                    exons = [exons[len(exons)-1-i] for i in range(len(exons))]

                variant_distribution = {key: [] for key in self.params["variant_test_targets"]}
                
                for j in range(len(mapped_index[i])):
                    if self.predictions.iloc[i].loc[self.params["mapping_target"]] != self.variants.iloc[mapped_index[i][j]].loc["ID:"+self.params["mapping_target"]]:
                        print("< mapping error occurred @map_variants: ", self.predictions.iloc[i].loc[self.params["mapping_target"]],
                              "/", self.variants.iloc[mapped_index[i][j]].loc["ID:"+self.params["mapping_target"]])

                    unassigned = True
                    for key_index, key in enumerate(self.params["variant_targets"]):
                        for k in range(len(exons)):
                            if ((k < len(exons)-1 and self.variants.iloc[mapped_index[i][j]].loc[self.params["variant_targets"][key]] >= int(exons[k])
                                and self.variants.iloc[mapped_index[i][j]].loc[self.params["variant_targets"][key]] < int(exons[k+1]))
                                or (k == len(exons)-1 and self.variants.iloc[mapped_index[i][j]].loc[self.params["variant_targets"][key]] >= int(exons[k]))):
                                if key_index == 0 and unassigned == False:
                                    print("< error @map_variants. feature was assigned to item twice.", key_index, key, "exon", k)
                                    exit()

                                if key_index == 0 and k == 0 and len(exons) > 1 :  stats["first"]  += 1; unassigned = False
                                if key_index == 0 and k == len(exons)-1         :  stats["last"]   += 1; unassigned = False
                                if key_index == 0 and k >= 1 and k < len(exons)-1: stats["middle"] += 1; unassigned = False
                                variant_targets[key][exons[k]].append(float(self.variants.iloc[mapped_index[i][j]].loc[key]))

                        if key_index == 0 and unassigned == True: stats["unassigned"] += 1
                                
                    # collect data for bootstrapping test to test for significant deviation between real features and theoretical distribution
                    [variant_distribution[variant_test_target].append(float(self.variants.iloc[mapped_index[i][j]].loc[variant_test_target]))
                     for variant_test_target in self.params["variant_test_targets"]]

                self.predictions.at[self.predictions.index[i], "variant targets"] = variant_targets

                # conduct bootstrapping test to test for significant deviation between real features and theoretical distribution
                for variant_test_target in self.params["variant_test_targets"]:
                    if len(self.params["masks"]) > 0:
                        predictions, predictions_by_exon, probabilities, probabilities_by_exon, self.predictions.at[self.predictions.index[i], "variant targets"] = self.mask_predictions(self.predictions.loc[self.predictions.index[i]].loc["predictions"],
                                                                                                                       self.predictions.loc[self.predictions.index[i]].loc["predictions_by_exon"],
                                                                                                                       variant_targets, self.seqs.loc[self.predictions.index[i]].loc["cds"],
                                                                                                                       self.predictions.index[i], self.params["mutation_stats_gene_target"], self.params["mutation_stats_pair_target"])

                        self.predictions.at[self.predictions.index[i], "predictions"]           = predictions
                        self.predictions.at[self.predictions.index[i], "predictions_by_exon"]   = predictions_by_exon
                        self.predictions.at[self.predictions.index[i], "probabilities"]         = probabilities
                        self.predictions.at[self.predictions.index[i], "probabilities_by_exon"] = probabilities_by_exon # <- until here

                    if self.params["apply_global_predictions"] == False: predictions = self.predictions.iloc[i].loc["predictions"]
                    else:                                                predictions = self.global_predictions

                    if len(variant_distribution[variant_test_target]) > 0:
                        self.predictions.at[self.predictions.index[i], "avg. model prediction"] = np.nanmean(predictions)
                        self.predictions.at[self.predictions.index[i], "avg. real prediction"]  = np.nanmean(variant_distribution[variant_test_target])
                        self.predictions.at[self.predictions.index[i], "median model prediction"] = np.nanmedian(predictions)
                        self.predictions.at[self.predictions.index[i], "median real prediction"]  = np.nanmedian(variant_distribution[variant_test_target])

                        self.predictions.at[self.predictions.index[i], "avg. real size"]        = len(variant_distribution[variant_test_target])
                        self.predictions.at[self.predictions.index[i], "block id"]              = self.variants.iloc[mapped_index[i][j]].loc["ID:gene symbol"]
                        self.predictions.at[self.predictions.index[i], "real predictions"]      = variant_distribution[variant_test_target]

                        # this probably works only for a single variant_test_target
                        for key in self.params["selection_cols"]:
                            self.predictions.at[self.predictions.index[i], key] = self.variants.iloc[mapped_index[i]][self.params["selection_cols"][key]].tolist()

                        if self.params["apply_mutation_stats"] == False and self.params["apply_mutation_weights"] == False and self.params["mutation_stats_ptc_weights"] == False: # <- added, recheck
                            test_results = self.conduct_test(predictions, variant_distribution[variant_test_target], variant_test_target) # <- added

                        if self.params["apply_mutation_stats"] == True or self.params["apply_mutation_weights"] == True or self.params["mutation_stats_ptc_weights"] == True: # recheck
                            test_results = self.conduct_test_w_mutation_stats(predictions, self.predictions.loc[self.predictions.index[i]].loc["probabilities"],
                                                                              variant_distribution[variant_test_target], variant_test_target) # <- added

                        for key in test_results:
                            self.predictions.at[self.predictions.index[i], key] = test_results[key]
            
            # <- added on 251111 
            elif len(mapped_index[i]) == 0:
                predictions, predictions_by_exon, probabilities, probabilities_by_exon, self.predictions.at[self.predictions.index[i], "variant targets"] = self.mask_predictions(self.predictions.loc[self.predictions.index[i]].loc["predictions"],
                                                                                                            self.predictions.loc[self.predictions.index[i]].loc["predictions_by_exon"],
                                                                                                            {variant_target: {} for variant_target in self.params["variant_targets"]}, self.seqs.loc[self.predictions.index[i]].loc["cds"],
                                                                                                            self.predictions.index[i], self.params["mutation_stats_gene_target"], self.params["mutation_stats_pair_target"])

                self.predictions.at[self.predictions.index[i], "predictions"]           = predictions
                self.predictions.at[self.predictions.index[i], "predictions_by_exon"]   = predictions_by_exon
                self.predictions.at[self.predictions.index[i], "probabilities"]         = probabilities
                self.predictions.at[self.predictions.index[i], "probabilities_by_exon"] = probabilities_by_exon # <- until here

            bar.next()
        bar.finish()
        print("< stats", stats)


    def print_params(self):
        params = json.dumps(self.params, indent=4)
        with open(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_params.json", "w") as file:
            file.write(params)       
    
    
    def save_stats(self):
        # marked if-condition (<-) added on 250528
        if len(self.params["masks"]) > 0:
            for key in self.stats:
                values          = pd.DataFrame({key: self.stats[key]})
                self.stats[key] = values[~values[key].isna()][key].tolist()
                print("<", values[values[key].isna()].shape[0], "removed for", key)

        # marked if-condition (<-) added on 250605, modified on 250821
        if ((self.params["apply_mutation_stats"] == True or (self.params["apply_mutation_stats"] == False and self.params["mutation_stats_ptc_weights"] == True))
            and self.params["variants_only"] == False):
            self.stats["values"] = self.calculate_distribution(self.stats["values"], self.stats["probabilities"])

            if self.params["exon_resolution"] == True:
                self.stats["first_exon_values"]  = self.calculate_distribution(self.stats["first_exon_values"], self.stats["first_exon_probabilities"])
                self.stats["last_exon_values"]   = self.calculate_distribution(self.stats["last_exon_values"], self.stats["last_exon_probabilities"])
                self.stats["middle_exon_values"] = self.calculate_distribution(self.stats["middle_exon_values"], self.stats["middle_exon_probabilities"])

        stats = json.dumps(self.stats, indent=4)

        with open(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_stats.json", "w") as file:
            file.write(stats)

        self.means = pd.DataFrame(self.means)
        self.means.to_csv(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_means.txt", index=False, sep=",")
    

    def show_selection(self, selection, header, project_mode=False):
        pd.set_option('display.max_columns', None)
        pd.options.display.width = 0
        print("<", header)
        # total values do not contain gene symbols
        selection["gene symbol"] = [json.loads(selection.iloc[i].loc["gene symbol"].replace("'", "\""))[0]
                                    if type(selection.iloc[i].loc["gene symbol"]) == str else None for i in range(selection.shape[0])]

        if project_mode == False:
            if self.params["selection_method"] == "fishers exact escape" or self.params["selection_method"] == "fishers exact target":
                selected_cols = [col for col in ["gene symbol", "gene id",
                                                 self.params["selection_method"]+"-statistic FEATURE:prediction",
                                                 self.params["selection_method"]+"-pvalue FEATURE:prediction",
                                                 "escape", "target", "relative escape", "relative target"]]

            else:
                selected_cols = [col for col in ["gene symbol", "gene id", "avg. model prediction", "avg. real prediction",
                                                 "median model prediction", "median real prediction", "avg. real size",
                                                 self.params["selection_method"]+"-pvalue FEATURE:prediction",
                                                 "escape", "target", "relative escape", "relative target"]]
            
        if project_mode == True:
            selected_cols = [*["gene symbol", "gene id"], *[col for col in selection.columns if ("escape" in col or "target" in col) and "relative" not in col]]
            
        selection = selection[selected_cols]
        print(selection)
        
        
    def store_test(self):
        cols = [
                "gene id", "transcript id", "avg. model prediction", "avg. real prediction", "median model prediction", "median real prediction", # <- added
                "avg. real size", "block id", "real predictions", "variant index", # <- added
                *self.params["selection_cols"], 
                *["binomial-pvalue "+variant_test_target for variant_test_target in self.params["variant_test_targets"]], # <-
                *["binomial-statistic "+variant_test_target for variant_test_target in self.params["variant_test_targets"]], # <-
                *["bootstrapping-pvalue "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact-statistic "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact-pvalue "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact escape-statistic "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact escape-pvalue "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact target-statistic "+variant_test_target for variant_test_target in self.params["variant_test_targets"]],
                *["fishers exact target-pvalue "+variant_test_target for variant_test_target in self.params["variant_test_targets"]]
                ]

        filtered_predictions = self.predictions[self.predictions["avg. real size"] > 0]
        filtered_predictions = filtered_predictions[cols]
        filtered_predictions.to_csv(self.newdir+self.params["os_sep"]+self.params["file_tag"]+"_selection_stats.txt", index=False)