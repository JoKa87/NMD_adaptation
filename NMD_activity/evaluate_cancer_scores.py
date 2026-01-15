import os

from evaluate_cancer_scores_utils import *
from tcga_tools_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


# parameters
params = {
         "analyze_class_genes"          : False,
         "block"                        : None,
         "categorization"               : {"FEATURE:frameshift_mutations_FEATURE:ptc_mutations": {"intercept": 0, "slope": 5},
                                           "FEATURE:ptc_mutations_FEATURE:frameshift_mutations": {"intercept": 0, "slope": 0.33}},
         "class_cutoff"                 : 4,
         "class_filter"                 : {"class 1": {"FEATURE:ptc_mutations": (1000000, 0), "total_nmd": (1000000, 0)},
                                           "class 2": {"FEATURE:ptc_mutations": (1000000, 0), "total_nmd": (1000000, 0)}},
         "class_selector"               : None, # ("FEATURE:ptc_mutations", "FEATURE:frameshift_mutations"), None
         "class_test_inclusive"         : False, # True if class 2 is not-exclusive (contains class 1)
         "clinical_data_path"           : parent_dir+r"\data\tcga_survival_data.txt",
         "cohorts"                      : 5, # if None, no cohorts are created for extracted cancer scores, else n cohorts are created
         "conduct_class_test"           : False,
         "contour_bins"                 : 5,
         "contour_scale"                : "equisize", # "equidistant" "equisize" "log"
         "data_dir"                     : parent_dir+r"data\2025-08-11_23-06-51_TCGA_NMD_targets_analysis_FPKM_exp_ccorr", # insert directory here
         "data_type"                    : "raw",
         "feature_filter"               : {"FEATURE:ptc_mutations": (1000000, 0), "total_nmd": (1000000, 0)},
         "features"                     : [
                                          "FEATURE:ptc_mutations2",
                                          "FEATURE:prediction",
                                          "fpkm_unstranded",
                                          "ID:cnv total",
                                          ],
         "fname"                        : "cancer_scores_TCGA_NMD_targets_analysis_FPKM_exp_ccorr_frameshift",
         "immune_data_path"             : parent_dir+r"\data\tcga_immune_scores.txt",
         "inversion_targets"            : ["fpkm_unstranded"],
         "os_sep"                       : "\\",
         "plot_mode"                    : "correlation", # "correlation" "distribution"
         "plot_type"                    : "correlation_plot", # "contour_plot" "correlation_plot"
         "projectwise"                  : True,
         "ptc_target"                   : "ptc_mutations",
         "selected_class"               : "class 1",
         "size_filter"                  : 0, # applied to mean/median values (cancer scores)
         "smoothing"                    : False,
         "smoothing_intervals"          : 30,
         "show_plots"                   : False,
         "store_results"                : True,
         "store_selection"              : True,
         "sum_targets"                  : ["FEATURE:expression"],
         "tag"                          : "ptc0",
         "use_inversion_targets"        : True,
         "use_sum_targets"              : False,
         "variant_path"                 : parent_dir+r"\data\tcga_variants.txt"
         }

if params["class_selector"] != None: params["class_test_inclusive"] = False
if params["data_type"] != "raw" and params["store_selection"] == True:
    print("< for data type other than raw, calculation of cancer scores is currently not supported.")
    params["store_selection"] = False

def _run_correlation(params):
    # calculate inversion from selected averages
    if params["use_inversion_targets"] == True:
        ecsu.get_inversions()

    # calculate sums from selected averages
    if params["use_sum_targets"] == True:
        ecsu.get_sums()

    # class selection
    if params["class_selector"] != None and params["data_type"] == "raw":
        ecsu.create_class_selection()

    elif params["class_selector"] != None and params["data_type"] == "raw":
        print("< class selection currently only working for raw data.")
        exit()

    ecsu.run_correlation()

    # show results
    pd.set_option('display.max_rows', None)
    pd.options.display.width = 0
    ecsu.stats_summary       = ecsu.stats_summary.sort_values(by="spearman-padj")
    print(ecsu.stats_summary[[col for col in ecsu.stats_summary if col not in ["pair", "x", "y"]]])

    if params["conduct_class_test"] == True and params["data_type"] == "raw":
        ecsu.class_test = ecsu.class_test.sort_values(by="padj")
        print(ecsu.class_test[[col for col in ecsu.class_test.columns if "class" not in col]])

    if params["analyze_class_genes"] == True:
        ecsu.analyze_class_genes()
        ecsu.class_genes_test = ecsu.class_genes_test.sort_values(by="padj")
        print(ecsu.class_genes_test[ecsu.class_genes_test["padj"] <= 0.05])

    # store results to file
    if params["store_results"] == True:
        ecsu.store_results()


ecsu = Evaluate_cancer_scores_utils(params)

# loading data
ecsu.load("cancer_scores")
cancer_scores = ecsu.cancer_scores

# run analysis
if params["plot_mode"] == "correlation":
    if params["projectwise"] == False:
        _run_correlation(params)

    else:
        cancer_scores = ecsu.cancer_scores
        projects      = ecsu.cancer_scores.drop_duplicates(subset=["project"])["project"].tolist()

        for project in projects:
            print("project", project)
            #if project not in ["TCGA-CHOL", "TCGA-DLBC", "TCGA-MESO", "TCGA-SARC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM"]:
            #if project in ["TCGA-STAD"]:
            ecsu.cancer_scores = cancer_scores[cancer_scores["project"] == project]
            params["block"]    = project

            try:
                _run_correlation(params)

            except:
                print("<", project, "is skipped")


if params["plot_mode"] == "distribution":
    ecsu.plot_distribution()