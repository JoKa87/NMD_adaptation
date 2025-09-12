import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from tcga_tools_utils import *


class Evaluate_cancer_scores_utils():
    def __init__(self, params):
        self.params                     = params

        self.class_plot                 = {}
        self.class_samples              = {"class 1": [], "class 2": []}
        self.cancer_scores              = pd.DataFrame()
        self.class_extracted_samples    = {}
        self.class_genes_test           = pd.DataFrame()
        self.class_test                 = pd.DataFrame()
        self.selected_sample_names      = []
        self.selected_samples           = None
        self.stats_summary              = pd.DataFrame()


    def analyze_class_genes_(self, ptc_variants):
        # create list of ptc variants occurring in affected samples and generate subframes for each class
        class_variants = {}
        selected_genes = []
        for class_key in self.class_extracted_samples:
            class_variants[class_key] = ptc_variants[ptc_variants["ID:case id"].isin(self.class_extracted_samples[class_key])]
            selected_genes.extend(class_variants[class_key].drop_duplicates(subset=["ID:gene id"])["ID:gene id"].tolist())

        selected_genes = np.unique(selected_genes)

        # initialize class test container
        class_index           = ["total", *selected_genes]
        self.class_genes_test = pd.DataFrame({"item": class_index, "class 1": [0 for _ in class_index], "class 2": [0 for _ in class_index],
                                              "statistic": [0 for _ in class_index], "p": [0 for _ in class_index], "padj": [0 for _ in class_index]}, index=class_index)
        
        # collect statistics
        bar = IncrementalBar(set_bar("collecting statistics"), max=len(selected_genes)*len(self.class_extracted_samples))
        for selected_gene in selected_genes:
            for class_key in self.class_extracted_samples:
                selected_variants                                  = class_variants[class_key][class_variants[class_key]["ID:gene id"] == selected_gene]
                selected_variants                                  = selected_variants.drop_duplicates(subset=["ID:case id"])
                self.class_genes_test.at[selected_gene, class_key] = selected_variants.shape[0]
                self.class_genes_test.at["total", class_key]      += selected_variants.shape[0]

                bar.next()
        bar.finish()

        
        # calculate statistics
        pvalues = {}
        bar = IncrementalBar(set_bar("calculating statistics"), max=self.class_genes_test.index.shape[0])
        for i in self.class_genes_test.index:
            if i != "total" and type(self.class_genes_test.loc[i].loc["class 1"]) != list:
                #print(i, [class_test.at[i, "class 1"], class_test.at["total", "class 1"]-class_test.at[i, "class 1"]], [class_test.at[i, "class 2"], class_test.at["total", "class 2"]-class_test.at[i, "class 2"]])
                result = stats.fisher_exact([[self.class_genes_test.at[i, "class 1"], self.class_genes_test.at["total", "class 1"]-self.class_genes_test.at[i, "class 1"]],
                                             [self.class_genes_test.at[i, "class 2"], self.class_genes_test.at["total", "class 2"]-self.class_genes_test.at[i, "class 2"]]])

                self.class_genes_test.at[i, "statistic"] = result.statistic
                self.class_genes_test.at[i, "p"]         = result.pvalue
                pvalues[i]                               = self.class_genes_test.loc[i].loc["p"]
            bar.next()
        bar.finish()
        
        padj = sm.stats.fdrcorrection([pvalues[key] for key in pvalues], alpha=0.05)[1]
        for i, index in enumerate(pvalues):
            self.class_genes_test.at[index, "padj"] = padj[i]
        

    def analyze_class_genes(self):
        ptc_variants = pd.read_csv(self.params["variant_path"], delimiter=",")

        for class_key in self.class_extracted_samples:
            self.class_extracted_samples[class_key] = np.unique(self.class_extracted_samples[class_key])

        if self.params["class_test_inclusive"] == True:
            self.class_extracted_samples["class 2"] = [self.class_extracted_samples["class 2"][j] for j in range(len(self.class_extracted_samples["class 2"]))
                                                       if self.class_extracted_samples["class 2"][j] not in self.class_extracted_samples["class 1"]]
        
        print("< sample sizes:", len(self.class_extracted_samples["class 1"]), "/", len(self.class_extracted_samples["class 2"]))
        self.analyze_class_genes_(ptc_variants)
    

    def conduct_class_test(self):
        pvalues1 = {}; pvalues2 = {}

        # if class 2 contains class 1, values must be subtracted
        if self.params["class_test_inclusive"] == True:
            self.class_test.at["total", "class 2"] -= self.class_test.loc["total"].loc["class 1"]

        for i in self.class_test.index:
            '''
            contingency table  (permutation yields same result):

            	        project	    other projects
            class1      
            class2
            '''

            if i != "total" and type(self.class_test.loc[i].loc["class 1"]) != list:
                # if class 2 contains class 1, values must be subtracted
                if self.params["class_test_inclusive"] == True: self.class_test.at[i, "class 2"] -= self.class_test.loc[i].loc["class 1"]

                result                             = stats.fisher_exact([[self.class_test.at[i, "class 1"], self.class_test.at["total", "class 1"]-self.class_test.at[i, "class 1"]],
                                                                         [self.class_test.at[i, "class 2"], self.class_test.at["total", "class 2"]-self.class_test.at[i, "class 2"]]])

                self.class_test.at[i, "statistic"] = result.statistic
                self.class_test.at[i, "p"]         = result.pvalue
                pvalues1[i]                        = self.class_test.loc[i].loc["p"]

            if type(self.class_test.loc[i].loc["class 1"]) == list and len(self.class_test.loc[i].loc["class 1"]) > 0 and len(self.class_test.loc[i].loc["class 2"]) > 0:
                # if class 2 contains class 1, values must be subtracted
                if self.params["class_test_inclusive"] == True: self.class_test.at[i, "class 2"] = [self.class_test.loc[i].loc["class 2"][j] for j in range(len(self.class_test.loc[i].loc["class 2"]))
                                                                                                    if self.class_test.loc[i].loc["class 2"][j] not in self.class_test.loc[i].loc["class 1"]]

                result                                = stats.kstest(self.class_test.loc[i].loc["class 1"], self.class_test.loc[i].loc["class 2"])
                self.class_test.at[i, "ks-statistic"] = result.statistic
                self.class_test.at[i, "ks-p"]         = result.pvalue
                result                                = stats.mannwhitneyu(self.class_test.loc[i].loc["class 1"], self.class_test.loc[i].loc["class 2"])
                u1                                    = result.statistic
                u2                                    = len(self.class_test.loc[i].loc["class 1"])*len(self.class_test.loc[i].loc["class 2"])-u1
                self.class_test.at[i, "statistic"]    = u1/u2
                self.class_test.at[i, "p"]            = result.pvalue
                pvalues2[i]                           = self.class_test.loc[i].loc["p"]

        # conduct Benjamini-Hochberg correction (projects and items considered separately)
        padj = sm.stats.fdrcorrection([pvalues1[key] for key in pvalues1], alpha=0.05)[1]
        for i, index in enumerate(pvalues1):
            self.class_test.at[index, "padj"] = padj[i]
            
        padj = sm.stats.fdrcorrection([pvalues2[key] for key in pvalues2], alpha=0.05)[1]
        for i, index in enumerate(pvalues2):
            self.class_test.at[index, "padj"] = padj[i]

        return
    

    def create_class_selection(self):
        x                 = []
        y                 = []
        full_sample_names = []

        for i in range(self.cancer_scores.shape[0]):
            values1, values2, _, sample_names = get_values(self.cancer_scores, self.params["feature_filter"], self.params["class_selector"][0], self.params["class_selector"][1], i,
                                                           self.params["ptc_target"], extract_ptc=True, get_sample_names=True)

            [x.append(values1[j]) for j in range(len(values1)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False]
            [y.append(values2[j]) for j in range(len(values2)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False]
            [full_sample_names.append(sample_names[j]) for j in range(len(sample_names)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False]

        x1, x2, y1, y2 = assign_cluster(x, y, self.params["categorization"][self.params["class_selector"][0]+"_"+self.params["class_selector"][1]]["intercept"],
                                        self.params["categorization"][self.params["class_selector"][0]+"_"+self.params["class_selector"][1]]["slope"])
        
        self.class_plot["x"]  = [x]
        self.class_plot["y"]  = [y] 
        self.class_plot["x1"] = [x1]
        self.class_plot["x2"] = [x2]
        self.class_plot["y1"] = [y1]
        self.class_plot["y2"] = [y2]

        params1, y_fit1 = fit(x1, y1)
        params2, y_fit2 = fit(x2, y2)
        plt.plot(x, y, "bo")
        plt.plot(x1, y_fit1, "--", color="red")
        plt.plot(x2, y_fit2, "--", color="magenta")

        equation = self.params["categorization"][self.params["class_selector"][0]+"_"+self.params["class_selector"][1]]
        plt.plot([0, max(x)], [0, equation["intercept"]+equation["slope"]*max(x)], "--", color="black")

        class1_x = []; class1_y = []; class2_x = []; class2_y = []
        for i in range(len(x)):
            dist1 = abs(y[i]-(params1[0]+params1[1]*x[i]))
            dist2 = abs(y[i]-(params2[0]+params2[1]*x[i]))

            if dist1 / dist2 >= self.params["class_cutoff"]:
                self.class_samples["class 1"].append(full_sample_names[i])
                plt.plot(x[i], y[i], "ro")
                class1_x.append(x[i])
                class1_y.append(y[i])

            if dist2 / dist1 >= self.params["class_cutoff"]:
                self.class_samples["class 2"].append(full_sample_names[i])
                plt.plot(x[i], y[i], "go")
                class2_x.append(x[i])
                class2_y.append(y[i])

        self.class_plot["class1_x"] = [class1_x]
        self.class_plot["class1_y"] = [class1_y]
        self.class_plot["class2_x"] = [class2_x]
        self.class_plot["class2_y"] = [class2_y]
        self.class_plot             = pd.DataFrame(self.class_plot)

        if self.params["show_plots"] == True: plt.show()

        self.selected_samples = self.class_samples[self.params["selected_class"]]
        print("< low class:", len(self.class_samples["class 1"]), ", high class:", len(self.class_samples["class 2"]))
        return
    

    def get_raw_values(self, feature1, feature2):
        x = []
        y = []
        z = []

        for i in range(self.cancer_scores.shape[0]):
            values1, values2, values3, sample_names = get_values(self.cancer_scores, self.params["feature_filter"], feature1, feature2, i, self.params["ptc_target"],
                                                                 class_samples=self.selected_samples, extract_ptc=True, get_sample_names=True)

            [x.append(values1[j]) for j in range(len(values1)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False and pd.isna(values3[j]) == False]
            [y.append(values2[j]) for j in range(len(values2)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False and pd.isna(values3[j]) == False]
            [z.append(values3[j]) for j in range(len(values3)) if pd.isna(values1[j]) == False and pd.isna(values2[j]) == False and pd.isna(values3[j]) == False]
            self.selected_sample_names.extend([sample_name for sample_name in sample_names if sample_name not in self.selected_sample_names])

        return x, y, z
    

    def get_stats(self, x, y, feature1, feature2):
        if len(x) > 1:
            if len([x[i] for i in range(len(x)) if x[i] == -1]) > 0 and len([y[i] for i in range(len(y)) if y[i] == -1]) > 0:
                print("< error with feature pair", feature1, "/", feature2, ": '-1' placeholder detected.")
            
            pearsonr  = stats.pearsonr(x, y)
            spearmanr = stats.spearmanr(x, y)

            if feature1+"_"+feature2 in self.stats_summary.index:
                self.stats_summary.at[feature1+"_"+feature2, "pearson-r"]  = pearsonr.statistic
                self.stats_summary.at[feature1+"_"+feature2, "pearson-p"]  = pearsonr.pvalue
                self.stats_summary.at[feature1+"_"+feature2, "spearman-r"] = spearmanr.statistic
                self.stats_summary.at[feature1+"_"+feature2, "spearman-p"] = spearmanr.pvalue
                self.stats_summary.at[feature1+"_"+feature2, "x"]          = x
                self.stats_summary.at[feature1+"_"+feature2, "y"]          = y

            elif feature2+"_"+feature1 in self.stats_summary.index:
                if (round(self.stats_summary.loc[feature2+"_"+feature1].loc["pearson-r"], 5) != round(pearsonr.statistic, 5)
                    or round(self.stats_summary.loc[feature2+"_"+feature1].loc["pearson-p"], 5) != round(pearsonr.pvalue, 5)
                    or round(self.stats_summary.loc[feature2+"_"+feature1].loc["spearman-r"], 5) != round(spearmanr.statistic, 5)
                    or round(self.stats_summary.loc[feature2+"_"+feature1].loc["spearman-p"], 5) != round(spearmanr.pvalue, 5)):
                    print("< inconsistent correlation detected @get_stats for", feature1+"_"+feature2)
                    print("  pearson-r", self.stats_summary.loc[feature2+"_"+feature1].loc["pearson-r"], "/", pearsonr.statistic,
                        "pearson-pvalue", self.stats_summary.loc[feature2+"_"+feature1].loc["pearson-p"], "/", pearsonr.pvalue)
                    print("  spearman-r", self.stats_summary.loc[feature2+"_"+feature1].loc["spearman-r"], "/", spearmanr.statistic,
                        "spearman-pvalue", self.stats_summary.loc[feature2+"_"+feature1].loc["spearman-p"], "/", spearmanr.pvalue)
                    

        else:
            print("< insufficient data size for feature pair", feature1, "/", feature2)
            

    def get_inversions(self):
        for target in self.params["inversion_targets"]:
            updated_values = []
            for i in range(self.cancer_scores.shape[0]):
                values        = self.cancer_scores.iloc[i].loc[target+"_"+self.params["data_type"]]
                ptc_mutations = self.cancer_scores.iloc[i].loc[self.params["ptc_target"]+"_"+self.params["data_type"]]
                values        = [values[j]*(-1) if pd.isna(values[j]) == False and pd.isna(ptc_mutations[j]) == False else None
                                 for j in range(len(values))]
                
                updated_values.append(values)

                if len(values) != len(ptc_mutations):
                    print("< size error occurred @get_inversions")
                    exit()

            self.cancer_scores[target+"_"+self.params["data_type"]] = updated_values
                

    def get_sums(self):
        for target in self.params["sum_targets"]:
            updated_values = []
            for i in range(self.cancer_scores.shape[0]):
                values        = self.cancer_scores.iloc[i].loc[target+"_"+self.params["data_type"]]
                ptc_mutations = self.cancer_scores.iloc[i].loc[self.params["ptc_target"]+"_"+self.params["data_type"]]
                values        = [values[j]*ptc_mutations[j] if pd.isna(values[j]) == False and pd.isna(ptc_mutations[j]) == False else values[j]
                                 for j in range(len(values))]
                
                updated_values.append(values)

                if len(values) != len(ptc_mutations):
                    print("< size error occurred @get_sums")
                    exit()

            self.cancer_scores[target+"_"+self.params["data_type"]] = updated_values
                

    # value conversion for mean or median data
    def get_values(self, feature1, feature2):
        if self.params["size_filter"] != None:
            selected_projects      = [i for i in range(self.cancer_scores.shape[0]) if self.cancer_scores.iloc[i].loc[feature1+"_count"] >= self.params["size_filter"]]
            selected_cancer_scores = self.cancer_scores.iloc[selected_projects]

            x = selected_cancer_scores[feature1+"_"+self.params["data_type"]]
            y = selected_cancer_scores[feature2+"_"+self.params["data_type"]]
            z = selected_cancer_scores[self.params["ptc_target"]+"_"+self.params["data_type"]]

        else:
            x = self.cancer_scores[feature1+"_"+self.params["data_type"]]
            y = self.cancer_scores[feature2+"_"+self.params["data_type"]]
            z = self.cancer_scores[self.params["ptc_target"]+"_"+self.params["data_type"]]

        if x[x.isna()].shape[0] > 0 or y[y.isna()].shape[0] > 0 or z[z.isna()].shape[0] > 0:
            print("< missing value detected @get_values")
            exit()

        if x[x == -1].shape[0] > 0 or y[y == -1].shape[0] > 0 or z[z == -1].shape[0] > 0:
            print("< '-1' found @get_values")
            exit()

        return x, y, z
    

    def init_stats(self):
        stats_list  = []
        stats_index = []
        for i in range(len(self.params["features"])):
            for j in range(i+1, len(self.params["features"])):
                if i != j:
                    stats_index.append(self.params["features"][i]+"_"+self.params["features"][j])
                    stats_list.append(None)

        self.stats_summary           = pd.DataFrame({"pair": stats_index, "pearson-r": list(stats_list), "pearson-p": list(stats_list), "pearson-padj": list(stats_list),
                                                     "spearman-r": list(stats_list), "spearman-p": list(stats_list), "spearman-padj": list(stats_list),
                                                     "x": list(stats_list), "y": list(stats_list)},
                                                     index=stats_index)

        self.class_extracted_samples = {"class 1": [], "class 2": []}
        class_index                  = [*[project for project in self.cancer_scores["project"]], *["total"], *[feature for feature in self.params["features"]]]
        self.class_test              = pd.DataFrame({"item"             : class_index,
                                                     "class 1"          : [0 if item not in self.params["features"] else [] for item in class_index],
                                                     "class 2"          : [0 if item not in self.params["features"] else [] for item in class_index],
                                                     "ks-statistic"     : [None for _ in class_index], "ks-p": [None for _ in class_index],
                                                     "statistic"        : [0 for _ in class_index], "p": [0 for _ in class_index], "padj": [0 for _ in class_index]}, index=class_index)
        return


    def load(self, mode):
        if mode == "cancer_scores":
            self.cancer_scores = pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["fname"], delimiter=",")
            self.cancer_scores = self.cancer_scores.sort_values("project", ascending=True)

            # load all lists stored as str to lists
            self.cancer_scores = convert_raw_data(self.cancer_scores)
            #self.cancer_scores = self.cancer_scores[self.cancer_scores["project"] == "TCGA-SKCM"]

            # temporary solution!
            #for i in range(self.cancer_scores.shape[0]):
                #print(i, self.cancer_scores.iloc[i].loc["age_at_initial_pathologic_diagnosis_raw"])
            #    self.cancer_scores.at[self.cancer_scores.index[i], "age_at_initial_pathologic_diagnosis_raw"] = [float(i) if pd.isna(i) == False else None for i in self.cancer_scores.iloc[i].loc["age_at_initial_pathologic_diagnosis_raw"]]
    

    def plot_contour(self, fig, ax, x, y, z, feature2, create_bar=False):
        max_x = max(x)
        max_y = max(y)
        max_z = max(z)

        contour_x, contour_y = np.meshgrid(np.arange(self.params["contour_bins"]), np.arange(self.params["contour_bins"]))

        data_3d = [[0 for _ in range(self.params["contour_bins"])] for _ in range(self.params["contour_bins"])]
        if self.params["contour_scale"] == "equidistant": index = np.arange(self.params["contour_bins"]+1)
        if self.params["contour_scale"] == "log":         index = [math.pow(2, i) for i in range(self.params["contour_bins"]+1)]; print(index)

        xticks = []

        for i in range(self.params["contour_bins"]):
            if self.params["contour_scale"] != "equisize":
                x_temp = [x[j] for j in range(len(z)) if z[j] <= (index[i+1])*max_z/self.params["contour_bins"] and z[j] > index[i]*max_z/self.params["contour_bins"]]
                y_temp = [y[j] for j in range(len(z)) if z[j] <= (index[i+1])*max_z/self.params["contour_bins"] and z[j] > index[i]*max_z/self.params["contour_bins"]]
                temp   = [[] for _ in range(self.params["contour_bins"])]
                [temp[min(len(temp)-1, int(self.params["contour_bins"]*y_temp[j]/max_y))].append(x_temp[j]) for j in range(len(x_temp))]

            if self.params["contour_scale"] == "equisize":
                sorted_z = np.argsort(z)
                x_temp   = [x[sorted_z[j]] for j in range(int(i*len(sorted_z)/self.params["contour_bins"]), int((i+1)*len(sorted_z)/self.params["contour_bins"]))]
                y_temp   = [y[sorted_z[j]] for j in range(int(i*len(sorted_z)/self.params["contour_bins"]), int((i+1)*len(sorted_z)/self.params["contour_bins"]))]
                xticks.append(np.mean([z[sorted_z[j]] for j in range(int(i*len(sorted_z)/self.params["contour_bins"]), int((i+1)*len(sorted_z)/self.params["contour_bins"]))]))
            
                sorted_y = np.argsort(y_temp)
                temp   = [[] for _ in range(self.params["contour_bins"])]
                [[temp[min(len(temp)-1, j)].append(x_temp[sorted_y[k]]) for k in range(int(j*len(sorted_y)/self.params["contour_bins"]), min(len(sorted_y), int((j+1)*len(sorted_y)/self.params["contour_bins"])))]
                 for j in range(self.params["contour_bins"])]
                yticks = []
                [yticks.append(np.mean([y_temp[sorted_y[k]] for k in range(int(j*len(sorted_y)/self.params["contour_bins"]), min(len(sorted_y), int((j+1)*len(sorted_y)/self.params["contour_bins"])))]))
                 for j in range(self.params["contour_bins"])]

            data_3d[i] = [np.mean(temp[j]) if len(temp[j]) > 0 else np.nan for j in range(self.params["contour_bins"])]


        data_3d = np.transpose(np.array(data_3d))
        ax.pcolor(contour_x, contour_y, data_3d, cmap=plt.get_cmap('viridis_r'), shading='nearest')
        ax.set_xlabel("ptc mutations")
        ax.set_ylabel(feature2)

        if self.params["contour_scale"] != "equisize":
            ax.set_xticks(np.arange(self.params["contour_bins"]), [round(i*max_z/self.params["contour_bins"], 3) for i in range(self.params["contour_bins"])], rotation="vertical")
            ax.set_yticks(np.arange(self.params["contour_bins"]), [round(i*max_y/self.params["contour_bins"], 3) for i in range(self.params["contour_bins"])])

        if self.params["contour_scale"] == "equisize":
            ax.set_xticks(np.arange(self.params["contour_bins"]), [round(xtick, 3) for xtick in xticks], rotation="vertical")
            ax.set_yticks(np.arange(self.params["contour_bins"]), [round(ytick, 3) for ytick in yticks])
        
        if create_bar == True:
            normalizer = mpl.colors.Normalize(vmin=0, vmax=max_x)
            cbar       = fig.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=plt.get_cmap('viridis_r')), location="right", shrink=0.5, format="%.1f", alpha=0.7)
            cbar.ax.locator_params(nbins=10)

        return fig, ax
    

    def plot_correlation(self, ax, x, y, feature2):           
        pearsonr  = stats.pearsonr(x, y)
        spearmanr = stats.spearmanr(x, y)
        stats_res = ("\n" + "  pearson:  " + str(round(pearsonr.statistic, 3)) + " " + str(round(pearsonr.pvalue, 3)) + "\n"
                            + "  spearman: " + str(round(spearmanr.statistic, 3)) + " " + str(round(spearmanr.pvalue, 3)))

        ax.plot(x, y, "bo", label=feature2+stats_res)
        ax.legend()
        return ax
    

    def plot_distribution(self):
        rows = int(self.cancer_scores.shape[0]/4)
        if rows < self.cancer_scores.shape[0]/4: rows += 1

        for feature in self.params["features"]:
            fig, ax = plt.subplots(rows, 4)
            fig.suptitle(feature, fontsize=14)

            col = 0; row = 0
            total = []
            for i in range(self.cancer_scores.shape[0]):
                x            = [value for value in self.cancer_scores.iloc[i].loc[feature+"_raw"] if pd.isna(value) == False]
                total.extend(x)
                label        = (self.cancer_scores.iloc[i].loc["project"] + "\n" + "  avg " + str(round(self.cancer_scores.iloc[i].loc[feature+"_median"], 4))
                                + "\n" + "  std " + str(round(self.cancer_scores.iloc[i].loc[feature+"_std"], 4)))

                ax[row][col].hist(x, bins=40, histtype="step", label=label)
                ax[row][col].legend(fontsize=8)

                col += 1
                if col == 4: col = 0; row += 1

            label        = ("total" + "\n" + "  avg " + str(round(np.mean(total), 4)) + "\n" + "  std " + str(round(np.std(total), 4)))
            
            ax[row][col].hist(total, bins=100, histtype="step", label=label) # [val for val in total if val <= 200]
            ax[row][col].legend(fontsize=8)
            plt.show()

        return
    

    def prepare_class_test(self, it):
        for i in range(self.cancer_scores.shape[0]):
            for class_key in self.params["class_filter"]:
                if len(self.class_samples[class_key]) > 0:
                    values, _, _, sample_names = get_values(self.cancer_scores, self.params["class_filter"][class_key], self.params["features"][it], self.params["features"][it],
                                                            i, self.params["ptc_target"], class_samples=self.class_samples[class_key], extract_ptc=True, get_sample_names=True)

                else:
                    values, _, _, sample_names = get_values(self.cancer_scores, self.params["class_filter"][class_key], self.params["features"][it], self.params["features"][it],
                                                            i, self.params["ptc_target"], extract_ptc=True, get_sample_names=True)

                self.class_extracted_samples[class_key].extend(sample_names)
                self.class_test.at[self.params["features"][it], class_key].extend([value for value in values if pd.isna(value) == False])

                if len(sample_names) != len(values):
                    print("< dimension error occurred @prepare_class_test:", len(sample_names), "/", len(values))

                if it == 0:
                    self.class_test.at[self.cancer_scores.iloc[i].loc["project"], class_key] = len([value for value in values if pd.isna(value) == False])
                    self.class_test.at["total", class_key]                                  += len([value for value in values if pd.isna(value) == False])
        
        return


    def run_correlation(self):
        # initialize containers
        self.init_stats()

        for i in range(len(self.params["features"])):
            rows = int(len(self.params["features"])/2)
            if rows == 1 or rows < len(self.params["features"])/2: rows += 1

            # initialize figure
            plt.close("all")
            fig, ax = plt.subplots(rows, 2, layout="constrained")
            fig.suptitle(self.params["features"][i], fontsize=14)

            # class test
            if self.params["conduct_class_test"] == True:
                self.prepare_class_test(i)

            # conduct correlation analysis
            col = 0; row = 0
            for j in range(len(self.params["features"])):                       
                if i != j:
                    if self.params["data_type"] == "median" or self.params["data_type"] == "mean":
                        x, y, z = self.get_values(self.params["features"][i], self.params["features"][j])

                    if self.params["data_type"] == "raw":
                        x, y, z = self.get_raw_values(self.params["features"][i], self.params["features"][j])

                    if self.params["smoothing"] == True:
                        _, x = smooth_data(z, x, intervals=self.params["smoothing_intervals"]) 
                        z, y = smooth_data(z, y, intervals=self.params["smoothing_intervals"])

                    # get statistics
                    self.get_stats(x, y, self.params["features"][i], self.params["features"][j])

                    # append plot
                    if self.params["plot_type"] == "correlation_plot":
                        ax[row][col] = self.plot_correlation(ax[row][col], x, y, self.params["features"][j])

                    if self.params["plot_type"] == "contour_plot":
                        if i < len(self.params["features"])-1 and j < len(self.params["features"])-1:
                            fig, ax[row][col] = self.plot_contour(fig, ax[row][col], x, y, z, self.params["features"][j])

                        elif i == len(self.params["features"])-1 and j < len(self.params["features"])-2:
                            fig, ax[row][col] = self.plot_contour(fig, ax[row][col], x, y, z, self.params["features"][j])

                        else:
                            fig, ax[row][col] = self.plot_contour(fig, ax[row][col], x, y, z, self.params["features"][j], create_bar=True)

                    col += 1
                    if col == 2: col = 0; row += 1

            if i == 0: print("< sample size", len(x), self.params["features"][i], self.params["features"][j])
            if self.params["show_plots"] == True: plt.show()

        # conduct class test (if selected), Bejamini-Hochberg correction included
        if self.params["conduct_class_test"] == True and self.params["data_type"] == "raw":
            self.conduct_class_test()

        elif self.params["conduct_class_test"] == True and self.params["data_type"] != "raw":
            print("< class test cannot be conducted for this datatype.")
    
        # conduct Bejamini-Hochberg correction
        self.stats_summary["pearson-padj"]  = sm.stats.fdrcorrection(self.stats_summary["pearson-p"], alpha=0.05)[1]
        self.stats_summary["spearman-padj"] = sm.stats.fdrcorrection(self.stats_summary["spearman-p"], alpha=0.05)[1]
       

    def store_results(self):
        new_dir = self.params["data_dir"]+self.params["os_sep"]+self.params["tag"]
        if os.path.isdir(new_dir) == False: os.mkdir(new_dir)

        with open(new_dir+"\\"+self.params["tag"]+"_params.json", "w") as f:
            f.write(json.dumps(self.params, indent=4))

        if self.params["block"] == None: tag = ""
        else:                            tag = "_"+self.params["block"]

        self.stats_summary.to_csv(new_dir+self.params["os_sep"]+"stats_summary"+tag+".txt", index=False, sep=",")
        self.stats_summary[[col for col in self.stats_summary.columns if col != "x" and col != "y"]].to_csv(new_dir+"\\stats_summary_overview"+tag+".txt", index=False, sep=",")
        
        if self.params["analyze_class_genes"] == True:
            self.class_genes_test.to_csv(new_dir+self.params["os_sep"]+"class_genes_test"+tag+".txt", sep=",", index=False)

        if self.params["class_selector"] != None:
            self.class_plot.to_csv(new_dir+self.params["os_sep"]+"class_plot"+tag+".txt", sep=",", index=False)

        if self.params["conduct_class_test"] == True:
            self.class_test.to_csv(new_dir+self.params["os_sep"]+"class_test"+tag+".txt", sep=",", index=False)
            self.class_test[[col for col in self.class_test.columns if col != "class 1" and col != "class 2"]].to_csv(new_dir+self.params["os_sep"]+"class_test_overview"+tag+".txt", sep=",", index=False)

        # store cancer scores with selected sample names
        if self.params["store_selection"] == True:
            feature_targets    = [
                                 #*["fpkm_unstranded", "ID:cnv total"],
                                 *self.params["features"],
                                 #*[col.replace("_mean", "") for col in self.cancer_scores.columns
                                 #  if "FEATURE" in col and "mean" in col and col.split("_")[len(col.split("_"))-2] != "n"]
                                 *[col.replace("_mean", "") for col in self.cancer_scores.columns
                                   if "FEATURE" in col and col in self.params["features"] and "mean" in col and col.split("_")[len(col.split("_"))-2] != "n"]
                                 ]

            clinical_features  = [#"age_at_initial_pathologic_diagnosis", "ajcc_pathologic_tumor_stage", "gender", "race",
                                  "OS", "OS.time", "DFI", "DFI.time", "DSS", "DSS.time", "PFI", "PFI.time"]
            immune_features    = ["Immune_score"]
            feature_targets    = [feature_target for feature_target in feature_targets
                                  if feature_target.replace("FEATURE:", "") not in clinical_features and feature_target not in immune_features]

            self.cancer_scores = convert_cancer_scores(self.cancer_scores, feature_targets, id_targets=["sample_names"], label_targets=[], datatype=self.params["data_type"],
                                                       randomize=False, selected_samples=self.class_extracted_samples[self.params["selected_class"]])

            immune_data        = pd.read_csv(self.params["immune_data_path"], delimiter=",")
            self.cancer_scores = map_data([self.cancer_scores, immune_data], immune_features, ["ID:sample_names", "case id"], "appending immune data", tag="FEATURE")

            clinical_data      = pd.read_csv(self.params["clinical_data_path"], delimiter=",")
            self.cancer_scores = map_data([self.cancer_scores, clinical_data], clinical_features, ["ID:sample_names", "case id"], "appending clinical data", tag="LABEL")

            # create random cohorts
            if self.params["cohorts"] != None:
                rand_index                       = np.arange(self.cancer_scores.shape[0])
                np.random.shuffle(rand_index)
                self.cancer_scores["rand_index"] = rand_index

                self.cancer_scores               = get_random_cohorts(self.cancer_scores, self.params["cohorts"], targets=["rand_index"])
                self.cancer_scores               = self.cancer_scores.drop(columns="rand_index")

            #self.cancer_scores.rename(columns={"LABEL:age_at_initial_pathologic_diagnosis": "FEATURE:age_at_initial_pathologic_diagnosis",
            #                                   "LABEL:ajcc_pathologic_tumor_stage": "FEATURE:ajcc_pathologic_tumor_stage",
            #                                   "LABEL:gender": "FEATURE:gender", 
            #                                   "LABEL:race": "FEATURE:race"}, inplace=True)
            self.cancer_scores.to_csv(new_dir+self.params["os_sep"]+"extracted_cancer_scores"+tag+".txt", sep=",", index=False)