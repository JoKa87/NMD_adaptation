import os
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from analyze_predictions_utils import *
from shared_utils import *


def _append_data(assembled_selection_stats, map, appended_data, message):
    appended_data["index"] = [i for i in range(appended_data.shape[0])]
    index                  = append_df_with_mapping([assembled_selection_stats, appended_data], "patient id", list(map.keys())[0],
                                                     "index", set_bar("mapping "+message+" data"), non_redundant=True, reverse=True)
    
    if len(map[list(map.keys())[0]]) == 0:
        map[list(map.keys())[0]] = [col for col in appended_data.columns if "LABEL" in col or "FEATURE" in col]

    for col in map[list(map.keys())[0]]:
        if col in appended_data.columns:
            assembled_selection_stats.insert(assembled_selection_stats.shape[1], col,
                                             [appended_data.iloc[int(i)].loc[col] if i != "-" else None for i in index])

    return assembled_selection_stats, map


def _evaluate_stats(assembled_selection_stats, map):
    stats = pd.DataFrame({col: [None for _ in map[list(map.keys())[0]]] for col in ["pearson", "pearson-p", "spearman", "spearman-p", "n"]},
                         index=[row for row in map[list(map.keys())[0]]])
    
    assembled_selection_stats = assembled_selection_stats.replace({"#NV": None})
    
    for row in map[list(map.keys())[0]]:
        if row in assembled_selection_stats.columns:
            assembled_selection_stats   = assembled_selection_stats.astype({row: 'float32'})
            selected_stats              = assembled_selection_stats[~assembled_selection_stats[row].isna()]

            if selected_stats.shape[0] > 0:
                pearson                     = scipy.stats.pearsonr(selected_stats["binomial statistic"], selected_stats[row])
                stats.at[row, "pearson"]    = pearson.statistic
                stats.at[row, "pearson-p"]  = pearson.pvalue

                spearman                    = scipy.stats.spearmanr(selected_stats["binomial statistic"], selected_stats[row])
                stats.at[row, "spearman"]   = spearman.statistic
                stats.at[row, "spearman-p"] = spearman.pvalue
                stats.at[row, "n"]          = selected_stats.shape[0]

    return stats


# parameters (default TCGA)
params = {
         "clinical_data_path"           : parent_dir+r"\data\tcga_survival_data.txt",
         "data_dir"                     : parent_dir+r"\data",
         "data_type"                    : "MSK", # "CPTAC3" "MSK" "TCGA"
         "map"                          : {"ID:sample_names": []}, # empty list defaults to selecting "FEATURE" and "LABEL" tagged cols
         "mode"                         : "analyze_selection",
         "os_sep"                       : "\\",
         "project_col"                  : "project",
         "significance_threshold"       : 1
         }

delimiter = "," # delimiter for clinical file

if params["data_type"] == "CPTAC3":
    params["clinical_path"] = parent_dir+r"\data\cptac3_survival_data.tsv"
    params["map"]           = {"cases.case_id": []}
    delimiter               = "\t"

if params["data_type"] == "MSK":
    params["clinical_path"] = parent_dir+r"\data\msk_chord_2024\msk_chord_survival_analysis.txt"
    params["map"]           = {"ID:PATIENT_ID": ["LABEL:OS_MONTHS", "LABEL:OS_STATUS", "FEATURE:expression", "FEATURE:frameshift",
                                                 "FEATURE:frameshift_mutations", "FEATURE:escape", "FEATURE:prediction", "FEATURE:target",
                                                 "FEATURE:ptc_mutations"]}
    params["project_col"]   = "cancer type"


if params["mode"] == "analyze_selection":
    fnames = os.listdir(params["data_dir"])
    fnames = [fname for fname in fnames if "selection_stats" in fname and "assembled_selection_stats" not in fname]

    # assemble patient selection stats
    assembled_selection_stats = {"patient id": [], "project": [], "binomial statistic": [], "binomial pvalue": []}

    bar = IncrementalBar(set_bar("assembling selection stats"), max=len(fnames))
    for fname in fnames:
        selection_stats = pd.read_csv(params["data_dir"]+params["os_sep"]+fname, sep=",")
        projects        = json.loads(selection_stats.iloc[0].loc[params["project_col"]].replace("'", "\""))

        if len(projects) == 1:
            assembled_selection_stats["patient id"].append(fname.split("_")[2])
            assembled_selection_stats["project"].append(projects[0]) # last entry does not provide project id
            assembled_selection_stats["binomial statistic"].append(selection_stats.iloc[-1].loc["binomial-statistic FEATURE:prediction"]*(-1))
            assembled_selection_stats["binomial pvalue"].append(selection_stats.iloc[-1].loc["binomial-pvalue FEATURE:prediction"])

        else:
            print("< error occurred. less or more than one project detected")
            print(" ", fname, "/", projects)
            exit()

        bar.next()
    bar.finish()

    # apply FDR-correction
    assembled_selection_stats["binomial padj"] = sm.stats.fdrcorrection(assembled_selection_stats["binomial pvalue"], alpha=params["significance_threshold"])[1]
    assembled_selection_stats                  = pd.DataFrame(assembled_selection_stats)
    assembled_selection_stats                  = assembled_selection_stats[assembled_selection_stats["binomial padj"] < params["significance_threshold"]]

    # plot distribution
    plt.hist(assembled_selection_stats["binomial statistic"], bins=40)
    plt.show()
    
    # clinical data
    clinical_data                            = pd.read_csv(params["clinical_path"], sep=delimiter)
    assembled_selection_stats, params["map"] = _append_data(assembled_selection_stats, params["map"], clinical_data, "clinical")

    # evaluate stats
    assembled_selection_stats.to_csv(params["data_dir"]+params["os_sep"]+"assembled_selection_stats.txt", sep=",")

    patient_stats = _evaluate_stats(assembled_selection_stats, params["map"])
    patient_stats.to_csv(params["data_dir"]+params["os_sep"]+"stats.txt", sep=",")
