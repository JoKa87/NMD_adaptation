import os
from survival_analysis_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def _run(sau):
    if sau.params["mode"] == "cox_regression":
        sau.load()
        cox_regression_results = sau.run_cox_regression()
        if sau.params["show_plot"] == True: sau.create_forest_plot(cox_regression_results)
        cox_regression_results.to_csv(sau.newdir+sau.params["os_sep"]+sau.params["fname"].split(".")[0]+"_cox_regression.txt")

    # marked (<-) added on 250521, modified on 250825
    if sau.params["mode"] == "kaplan_meier_estimator" or sau.params["mode"] == "kaplan_meier_max_stats":
        sau.load()
        km_results = sau.kaplan_meier_estimator()
        km_results.to_csv(sau.newdir+sau.params["os_sep"]+sau.params["fname"].split(".")[0]+"_kaplan_meier.txt", index=False)

    # return cox_regression_results # <- removed on 250826


def main():
    # parameters
    params = {
              "cox_filter"          : {"class": "ID:CANCER_TYPE"}, #"value": {"FEATURE:ptc_mutations": [1, 2, 5, 10, 20]}},
              "cox_log"             : True,
              "cox_mode"            : "multivariate", # "univariate" "multivariate"
              "data_dir"            : parent_dir+r"data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_projectwise", # put here the folder of target analysis (result from evaluate_cancer_scores)
              "fname"               : "selection_stats.txt",
              # key specifies splitting feature, value specifies splits (minimum value is 2, for >= 2 only extreme splits are considered (e.g. 0 and 2 for 3 splits))
              "km_filter"           : {"FEATURE:prediction": 2}, # <- added on 250521 
              "km_percentiles"      : [i for i in range(10, 90, 1)], # <- added on 250825 
              "label_ids"           : ["LABEL:OS", "LABEL:OS.time"],
              "max_stat_steps"      : 2, # <- added on 250825 
              "mode"                : "kaplan_meier_max_stats", # "cox_regression" "kaplan_meier_estimator" "kaplan_meier_max_stats"
              "os_sep"              : "\\",
              "placeholder_exceptions": ["FEATURE:fpkm_unstranded", "FEATURE:MSI_SCORE"], # defines categories to be excluded from placeholder analysis because negative values can exist
              "projectwise"         : False,
              # marked with <- edited on 250408 (added) to conduct survival analysis with selected features only
              "selected_features"   : ["FEATURE:ptc_mutations", "FEATURE:prediction"], # 'FEATURE:expression', 'FEATURE:ptc_mutations2', 'FEATURE:escape', 'FEATURE:target', 'FEATURE:prediction', 'FEATURE:fpkm_unstranded' 'FEATURE:expression', 'FEATURE:ptc_mutations2', 'FEATURE:exon', 'FEATURE:age_at_initial_pathologic_diagnosis', 'FEATURE:fpkm_unstranded', 'FEATURE:cnv total'], # ["FEATURE:expression", "FEATURE:target", "FEATURE:cnv total", "FEATURE:fpkm_unstranded", "FEATURE:ptc_mutations2"], # [] to switch off <-
              "separator"           : ",",
              "show_plot"           : False, # <- added on 250523
              "tag"                 : "msk",
              "target_dir"          : parent_dir+r"\data",
              "type"                : "msk_chord", # "mskcc" "msk_chord" "tcga" "tcga_cancer_scores"
             }

    
    # default values for tcga
    if params["type"] == "msk_chord":
        params["data_dir"]  = parent_dir+r"\data"
        params["fname"]     = "msk_chord_survival_analysis.txt"
        params["label_ids"] = ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"]


    if params["projectwise"] == False:
        if params["type"] != "tcga_cancer_scores":
            sau = Survival_analysis_utils(params)
            _run(sau)

        if params["type"] == "tcga_cancer_scores":
            subdirs         = os.listdir(params["data_dir"])
            subdirs         = [subdir for subdir in subdirs
                               if os.path.isdir(params["data_dir"]+params["os_sep"]+subdir) == True and "smoothing" not in subdir and "projectwise" not in subdir] # <- "smoothing" and "projectwise" conditions added on 250624
            data_dir        = params["data_dir"]
            params["fname"] = "extracted_cancer_scores.txt"
            params["tag"]   = params["fname"].split(".")[0]

            for subdir in subdirs:
                print("<", data_dir+params["os_sep"]+subdir)
                params["data_dir"] = data_dir+params["os_sep"]+subdir
                sau = Survival_analysis_utils(params, newdir=params["data_dir"])
                _run(sau)

    else:
        fnames = os.listdir(params["data_dir"])

        if params["type"] == "tcga_cancer_scores":
            fnames = [fname for fname in fnames if "extracted_cancer_scores" in fname]

            full_cox_regression_results = [] # <- added on 250605
            for fname in fnames:
                params["fname"] = fname
                params["tag"]   = fname.split(".")[0].split("_")[-1]

                sau = Survival_analysis_utils(params, newdir=params["data_dir"])
                try:
                    cox_regression_results       = _run(sau)
                    cox_regression_results.index = [params["tag"]+"_"+cox_regression_results.index[i] for i in range(cox_regression_results.shape[0])] # <- added on 250605
                    full_cox_regression_results.append(cox_regression_results) # <- added on 250605

                except:
                    print("<", fname, "is skipped.")

            full_cox_regression_results = pd.concat(full_cox_regression_results).sort_index() # <- added on 250605
            sau.params["tag"]      = sau.params["data_dir"].split(sau.params["os_sep"])[-1] # <- added on 250605
            sau.create_newdir() # <- added on 250605
            pd.set_option('display.max_rows', None)
            print(full_cox_regression_results[[True if "FEATURE:ptc_mutations2" in i else False for i in full_cox_regression_results.index]])
            full_cox_regression_results.to_csv(sau.newdir+sau.params["os_sep"]+sau.params["tag"]+"_cox_regression.txt", sep=",") # <- added on 250605

if __name__ == '__main__':
    main()