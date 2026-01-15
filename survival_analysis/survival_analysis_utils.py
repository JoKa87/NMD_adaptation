from copy import deepcopy
from datetime import datetime
import json
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter # <- added on 250521
from lifelines.statistics import logrank_test # <- added on 250521
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import platform
import statsmodels.api as sm
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

sys.path.insert(0, parent_dir+"\\plots")
from plot_utils import *
sys.path.insert(0, parent_dir+"\\shared")
from shared_utils import *
sys.path.insert(0, parent_dir+"\\NMD_activity")
from tcga_tools_utils import *


class Survival_analysis_utils(Plot_utils):
    def __init__(self, params, newdir=None):
        self.params         = params
        self.data           = pd.DataFrame()

        os_name = platform.system()
        self.os_sep = "//"
        if os_name == "Windows":
            self.os_sep = "\\"

        if newdir == None: self.create_newdir()
        else:              self.newdir = newdir

        # print parameters to file
        with open(self.newdir+self.os_sep+"params.json", "w") as f:
            f.write(json.dumps(self.params, indent=4))


    def create_forest_plot(self, cox_regression_results):
        center = 1
        if self.params["cox_log"] == True: center = 0

        features = {
                    "categories"    : pd.unique([col.split("_")[1] for col in cox_regression_results.columns]),
                    "center"        : center,
                    "error_type"    : "absolute",
                    "pvalue_col"    : "p-value",
                    "xcol"          : "HR",
                    "xerror"        : ["HR-lower 95%", "HR-upper 95%"],
                    "xlabel"        : "HR"
                    }

        xlabel = features["xlabel"]

        fig = plt.figure(constrained_layout=True)
        gs  = fig.add_gridspec(int(cox_regression_results.shape[0]/2)+1, 2)
        dim = int(cox_regression_results.shape[0]/2)
        if 2*dim < cox_regression_results.shape[0]: dim += 1

        subplots = []
        for i in range(dim):
            subplots.append(fig.add_subplot(gs[i, 0:1]))
            subplots.append(fig.add_subplot(gs[i, 1:2]))

        for i in range(cox_regression_results.shape[0]):
            data = pd.DataFrame({"category"             : [category for category in features["categories"]],
                                 features["pvalue_col"] : [cox_regression_results.iloc[i].loc[col] for col in cox_regression_results.columns
                                                           if col.split("_")[0] == features["pvalue_col"]],
                                 features["xcol"]       : [cox_regression_results.iloc[i].loc[col] for col in cox_regression_results.columns
                                                           if col.split("_")[0] == features["xcol"]],
                                 features["xerror"][0]  : [cox_regression_results.iloc[i].loc[col] for col in cox_regression_results.columns
                                                           if col.split("_")[0] == features["xerror"][0]],
                                 features["xerror"][1]  : [cox_regression_results.iloc[i].loc[col] for col in cox_regression_results.columns
                                                           if col.split("_")[0] == features["xerror"][1]]})
            features["xlabel"] = xlabel + " (" + cox_regression_results.index[i].replace("FEATURE:", "") + ")"
            subplots[i]        = self.forest_plot(subplots[i], data, features)
        
        plt.show()

        fig.savefig(self.newdir+self.params["os_sep"]+"cox_regression_results.png", dpi=300)


    def create_newdir(self):
        datestr = str(datetime.now())
        index   = datestr.index(".")
        dirname = datestr[:index].replace(" ", "_").replace(":", "-")
        newdir  = self.params["target_dir"] + self.os_sep + dirname + "_" + self.params["tag"]
        if not os.path.exists(newdir): os.mkdir(newdir)
        # marked (<-) added / removed on 250605
        # return newdir # <- removed
        self.newdir = newdir # <- added
    

    def cox_regression(self, data=pd.DataFrame(), verbose=True):
        if data.shape[0] == 0: data = self.data
        if data.shape[0] == 0: print("< empty dataframe @cox_regression")

        features = [col for col in data.columns if ("FEATURE" in col and len(self.params["selected_features"]) == 0) or col in self.params["selected_features"]] # <- new
        scores   = pd.DataFrame({col: [None for _ in features] for col in ["HR", "HR-lower 95%", "HR-upper 95%", "p-value", "p-adj", "concordance", "sample size"]}, index=features) # <- added / modified on 250823
        print("<", self.params["fname"])

        # univariate analysis
        if self.params["cox_mode"] == "univariate":
            for feature in features:
                selected_data = data[~data[feature].isna()]

                if selected_data.shape[0] > 1 and selected_data[feature].min() != selected_data[feature].max():
                    cph = CoxPHFitter()
                    cph.fit(selected_data[[feature, *self.params["label_ids"]]], duration_col=self.params["label_ids"][1], event_col=self.params["label_ids"][0])
                    
                    if self.params["cox_log"] == False:
                        scores.at[feature, "HR"]           = cph.summary.loc[feature].loc["exp(coef)"]
                        scores.at[feature, "HR-lower 95%"] = cph.summary.loc[feature].loc["exp(coef) lower 95%"]
                        scores.at[feature, "HR-upper 95%"] = cph.summary.loc[feature].loc["exp(coef) upper 95%"]

                    if self.params["cox_log"] == True:
                        scores.at[feature, "HR"]           = cph.summary.loc[feature].loc["coef"]
                        scores.at[feature, "HR-lower 95%"] = cph.summary.loc[feature].loc["coef lower 95%"]
                        scores.at[feature, "HR-upper 95%"] = cph.summary.loc[feature].loc["coef upper 95%"]

                    scores.at[feature, "p-value"]     = cph.summary.loc[feature].loc["p"]
                    scores.at[feature, "concordance"] = cph.concordance_index_
                    scores.at[feature, "sample size"] = selected_data.shape[0]

        selected_data = data.dropna(subset=features)
    
        # multivariate analysis
        if self.params["cox_mode"] == "multivariate":
            cph = CoxPHFitter()
            if selected_data.shape[0] > 1:
                cph.fit(selected_data[[*features, *self.params["label_ids"]]], duration_col=self.params["label_ids"][1],
                        event_col=self.params["label_ids"][0])
                
                for feature in features:
                    if self.params["cox_log"] == False:
                        scores.at[feature, "HR"]           = cph.summary.loc[feature].loc["exp(coef)"]
                        scores.at[feature, "HR-lower 95%"] = cph.summary.loc[feature].loc["exp(coef) lower 95%"]
                        scores.at[feature, "HR-upper 95%"] = cph.summary.loc[feature].loc["exp(coef) upper 95%"]

                    if self.params["cox_log"] == True:
                        scores.at[feature, "HR"]           = cph.summary.loc[feature].loc["coef"]
                        scores.at[feature, "HR-lower 95%"] = cph.summary.loc[feature].loc["coef lower 95%"]
                        scores.at[feature, "HR-upper 95%"] = cph.summary.loc[feature].loc["coef upper 95%"]

                    scores.at[feature, "p-value"]     = cph.summary.loc[feature].loc["p"]
                    scores.at[feature, "concordance"] = cph.concordance_index_
                    scores.at[feature, "sample size"] = selected_data.shape[0]
        
        if verbose == True: print(scores.sort_index())
        return scores.sort_index()


    # <- function re-integrated in altered fashion on 250521
    def _kaplan_meier_estimator(self, x, y, ax, tag):
        kmf = KaplanMeierFitter(label=tag)
        kmf.fit(x, y)
        return kmf.plot(), pd.concat([kmf.survival_function_, kmf.confidence_interval_survival_function_])


    # <- function re-integrated in altered fashion on 250521
    def kaplan_meier_estimator(self):
        km_results = {"class": [], "cutoff": [], "size 1": [], "size 2": [], "statistic": [], "pvalue": []} # <- added on 250825

        # if filters are defined, multiple runs are conducted to perform single output
        if "class" in self.params["cox_filter"]: classes = self.data.sort_values(by=self.params["cox_filter"]["class"]).drop_duplicates(subset=self.params["cox_filter"]["class"])[self.params["cox_filter"]["class"]].tolist()
        else:                                    classes = []

        for class_key in [*classes, "total"]:
            if class_key != "total": data = self.data[self.data[self.params["cox_filter"]["class"]] == class_key]
            else:                    data = self.data
            
            print("<", str(class_key)+" ("+str(data.shape[0])+")")
            data    = data[~data[list(self.params["km_filter"].keys())[0]].isna()]
            data    = data.sort_values(by=list(self.params["km_filter"].keys())[0])

            if data.shape[0] > 1:
                fig, ax = plt.subplots(1)

                if self.params["mode"] != "kaplan_meier_max_stats":
                    # <- added on 251101
                    if type(self.params["km_filter"][list(self.params["km_filter"].keys())[0]]) != list:
                        index   = split_index(data.shape[0], list(self.params["km_filter"].values())[0])
                        cutoffs = (data.iloc[index[0][1]].loc[list(self.params["km_filter"].keys())[0]], data.iloc[index[-1][0]].loc[list(self.params["km_filter"].keys())[0]])
                        times1  = data[data[list(self.params["km_filter"].keys())[0]] <= cutoffs[0]][self.params["label_ids"][1]]
                        status1 = data[data[list(self.params["km_filter"].keys())[0]] <= cutoffs[0]][self.params["label_ids"][0]]
                        times2  = data[data[list(self.params["km_filter"].keys())[0]] >= cutoffs[1]][self.params["label_ids"][1]]
                        status2 = data[data[list(self.params["km_filter"].keys())[0]] >= cutoffs[1]][self.params["label_ids"][0]]

                    else:
                        cutoffs = self.params["km_filter"][list(self.params["km_filter"].keys())[0]]
                        cutoffs = (cutoffs[0], cutoffs[1])
                        times1  = data[data[list(self.params["km_filter"].keys())[0]] <= cutoffs[0]][self.params["label_ids"][1]]
                        status1 = data[data[list(self.params["km_filter"].keys())[0]] <= cutoffs[0]][self.params["label_ids"][0]]
                        times2  = data[data[list(self.params["km_filter"].keys())[0]] > cutoffs[1]][self.params["label_ids"][1]]
                        status2 = data[data[list(self.params["km_filter"].keys())[0]] > cutoffs[1]][self.params["label_ids"][0]]

                    if times1.shape[0] > 0 and times2.shape[0]> 0:
                        ax      = self._kaplan_meier_estimator(times1, status1, ax, list(self.params["km_filter"].keys())[0]+" <= "+str(cutoffs[0]))
                        ax      = self._kaplan_meier_estimator(times2, status2, ax, list(self.params["km_filter"].keys())[0]+" >(=) "+str(cutoffs[1]))

                        results = logrank_test(times1, times2, event_observed_A=status1, event_observed_B=status2)
                        km_results["class"].append(class_key) # <- added on 250825
                        km_results["size 1"].append(times1.shape[0]) # <- added on 251101
                        km_results["size 2"].append(times2.shape[0]) # <- added on 251101
                        km_results["cutoff"].append(cutoffs) # <- added on 250825
                        km_results["statistic"].append(results.test_statistic) # <- added on 250825
                        km_results["pvalue"].append(results.p_value) # <- added on 250825

                        plt.title(str(class_key)+", pvalue:"+str(round(results.p_value, 4)))
                        if self.params["show_plot"] == True: plt.show()
                
                else:
                    observed_opt_stat, observed_opt_pvalue, observed_opt_cutoff = self.kaplan_meier_max_stats(data)

                    # print results
                    print("< statistic:", observed_opt_stat, "pvalue:", observed_opt_pvalue, "cutoff:", observed_opt_cutoff)
                    fig, ax = plt.subplots(1)
                    times1  = data[data[list(self.params["km_filter"].keys())[0]] <= observed_opt_cutoff][self.params["label_ids"][1]]
                    status1 = data[data[list(self.params["km_filter"].keys())[0]] <= observed_opt_cutoff][self.params["label_ids"][0]]
                    times2  = data[data[list(self.params["km_filter"].keys())[0]] > observed_opt_cutoff][self.params["label_ids"][1]]
                    status2 = data[data[list(self.params["km_filter"].keys())[0]] > observed_opt_cutoff][self.params["label_ids"][0]]
                    ax      = self._kaplan_meier_estimator(times1, status1, ax, list(self.params["km_filter"].keys())[0]+" <= "+str(observed_opt_cutoff))
                    ax      = self._kaplan_meier_estimator(times2, status2, ax, list(self.params["km_filter"].keys())[0]+" >= "+str(observed_opt_cutoff))

                    results = logrank_test(times1, times2, event_observed_A=status1, event_observed_B=status2)
                    km_results["class"].append(class_key)
                    km_results["size 1"].append(times1.shape[0])
                    km_results["size 2"].append(times2.shape[0])
                    km_results["cutoff"].append(observed_opt_cutoff)
                    km_results["statistic"].append(observed_opt_stat)
                    km_results["pvalue"].append(observed_opt_pvalue)

                    if results.test_statistic != observed_opt_stat:
                        print("< mismatching test statistics:", results.test_statistic, "/", observed_opt_stat)

                    plt.title("statistic: "+str(round(observed_opt_stat, 4))+"\n pvalue:"+str(round(observed_opt_pvalue, 4))+"\n opt. cutoff:"+str(round(observed_opt_cutoff, 4)))
                    if self.params["show_plot"] == True: plt.show()

        return pd.DataFrame(km_results)
        

    # <- added on 250825
    def _kaplan_meier_max_stats(self, data, cutoffs):
        max_stat = -np.inf; max_cutoff = None

        for cutoff in cutoffs:
            times1  = data[data[list(self.params["km_filter"].keys())[0]] <= cutoff][self.params["label_ids"][1]]
            status1 = data[data[list(self.params["km_filter"].keys())[0]] <= cutoff][self.params["label_ids"][0]]
            times2  = data[data[list(self.params["km_filter"].keys())[0]] > cutoff][self.params["label_ids"][1]]
            status2 = data[data[list(self.params["km_filter"].keys())[0]] > cutoff][self.params["label_ids"][0]]
            results = logrank_test(times1, times2, event_observed_A=status1, event_observed_B=status2)

            if results.test_statistic > max_stat:
                max_stat   = results.test_statistic
                max_cutoff = cutoff

        return max_stat, max_cutoff


    # <- added on 250825
    def kaplan_meier_max_stats(self, data):
        # define cutoffs
        cutoffs = np.percentile(data[list(self.params["km_filter"].keys())[0]], self.params["km_percentiles"])
        
        # find optimum cutoff
        observed_opt_stat, observed_opt_cutoff = self._kaplan_meier_max_stats(data, cutoffs)

        # determine max stats
        max_stats = []
        bar = IncrementalBar(set_bar("conducting max stat analysis"), max=self.params["max_stat_steps"])
        for _ in range(self.params["max_stat_steps"]):
            permutated_data                              = data.copy()
            permutated_data[self.params["label_ids"][0]] = np.random.permutation(permutated_data[self.params["label_ids"][0]].values)
            permutated_data[self.params["label_ids"][1]] = np.random.permutation(permutated_data[self.params["label_ids"][1]].values)
            max_stat, _                                  = self._kaplan_meier_max_stats(permutated_data, cutoffs)
            max_stats.append(max_stat)
            bar.next()
        bar.finish()

        # determine p-value
        observed_opt_pvalue = (np.sum(max_stats >= observed_opt_stat) + 1) / (self.params["max_stat_steps"] + 1)
        return observed_opt_stat, observed_opt_pvalue, observed_opt_cutoff


    def load(self):
        if os.path.isfile(self.params["data_dir"]+self.os_sep+self.params["fname"]) == True:
            self.data = pd.read_csv(self.params["data_dir"]+self.os_sep+self.params["fname"], delimiter=self.params["separator"], index_col=False)

        else:
            print("<", self.params["data_dir"]+self.os_sep+self.params["fname"], "not found.")
            exit()

        # check if input contains invalid values (included because placeholder system was changed to None/np.nan with certain exception remaining with "-")
        # '-1' is not a valid placeholder anymore!
        for col in self.data.columns:
            if ("FEATURE" in col or "LABEL" in col) and col not in self.params["placeholder_exceptions"]:
                if self.data[self.data[col] == -1].shape[0] > 0:
                    print("<", self.data[self.data[col] == -1].shape[0], "times '-1' was found in", col+".")
                    exit()

                if self.data[self.data[col] == "-"].shape[0] > 0:
                    print("<", self.data[self.data[col] == "-"].shape[0], "times '-' was found in", col+".")
                    exit()

        init_shape = self.data.shape[0]
        self.data  = self.data[~self.data[self.params["label_ids"][0]].isna()]
        self.data  = self.data[~self.data[self.params["label_ids"][1]].isna()]
        print("< removal of missing labels reduced data from", init_shape, "to", self.data.shape[0])

        # if features contain strings, conduct encoding
        str_features = [col for col in self.data.columns if "FEATURE" in col and self.data[col].dtype == object]

        if len(str_features) > 0:
            for str_feature in str_features:
                unique_features        = self.data.drop_duplicates(subset=str_feature)[str_feature].tolist()
                self.data[str_feature] = [unique_features.index(self.data.iloc[i].loc[str_feature]) for i in range(self.data.shape[0])]

            print("<", len(str_features), "feature(s) was/were found that contained strings and were encoded as integers.")


    def _run_cox_regression(self, cox_regression_results, scores, tag=""):
        if len(tag) > 0: tag = "_" + tag
        for index in scores.index:
            for col in scores.columns:
                if col+tag in cox_regression_results: cox_regression_results[col+tag].append(scores.loc[index].loc[col])
                else:                                 cox_regression_results[col+tag] = [scores.loc[index].loc[col]]

        return cox_regression_results 
    

    def run_cox_regression(self):
        cox_regression_results = {}

        # if filters are defined, multiple runs are conducted to perform single output
        if "class" in self.params["cox_filter"]: classes = self.data.sort_values(by=self.params["cox_filter"]["class"]).drop_duplicates(subset=self.params["cox_filter"]["class"])[self.params["cox_filter"]["class"]].tolist()
        else:                                    classes = []

        for class_key in [*classes, "total"]:
            if class_key != "total": data = self.data[self.data[self.params["cox_filter"]["class"]] == class_key]
            else:                    data = self.data
            
            print("<", str(class_key)+" ("+str(data.shape[0])+")")
            scores = self.cox_regression(data)
            cox_regression_results = self._run_cox_regression(cox_regression_results, scores, str(class_key))

            # create second local copy of data
            value_data = data

            if "value" in self.params["cox_filter"]:
                for filter_value in self.params["cox_filter"]["value"][list(self.params["cox_filter"]["value"].keys())[0]]:
                    value_data             = value_data[value_data[list(self.params["cox_filter"]["value"].keys())[0]] >= filter_value]
                    print("<", class_key, str(filter_value)+" ("+str(value_data.shape[0])+")") # <- added on 250704
                    scores                 = self.cox_regression(value_data)
                    cox_regression_results = self._run_cox_regression(cox_regression_results, scores, class_key+"_"+str(filter_value))

        # <- added on 250823 to conduct FDR correction, from here
        cox_regression_results = pd.DataFrame(cox_regression_results, index=scores.index).sort_index()

        cols = np.unique([col.split("_")[1] for col in cox_regression_results.columns])
        for feature in self.params["selected_features"]:
            pvalues = []; total_pvalues = []

            for col in cols:
                if "total" not in col: pvalues.append(cox_regression_results.loc[feature].loc["p-value_"+col]); print(col)
                else:                  total_pvalues.append(cox_regression_results.loc[feature].loc["p-value_"+col])

            pvalues = [*sm.stats.fdrcorrection(pvalues, alpha=0.05)[1], *total_pvalues]
            
            for i, col in enumerate(cols):
                cox_regression_results.at[feature, "p-adj_"+col] = pvalues[i]
                
        pd.set_option('display.max_columns', None)
        print(cox_regression_results[[col for col in cox_regression_results.columns if "p-adj" in col]])
        return cox_regression_results