import os
import platform
from multiprocessing import Process, Lock
import sys
import time

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from analyze_targets_utils import *
from extract_mutations_utils import *
from shared_utils import *


# parameters
params = {
            "appended_col"               : "NMD_targets",
            # marked (<-) added on 250429 to allow skipping of append_features
            "append_features"            : True, # <- added
            "append_sample_scores"       : False,
            "apply_filters"              : True, # for full test run @tcga_tools, filters need to be switched off (typically switched on)
            "clinical_data_path"         : parent_dir+r"\data\tcga_survival_data.txt",
            "cluster_key"                : None, # "Firehouse", None
            "create_calculations"        : True,
            "create_dataset"             : False,
            "data_dir"                   : parent_dir+r"\data", # here, TCGA data should be located
            "extended_features_path"     : parent_dir+r"\data\tcga_variants.txt",
            # marked (<-) added on 250429 to facilitate off-switching of feature integration
            "extend_features"            : True, # <- added
            "file_tag"                   : "_TCGA_NMD_targets_analysis_FPKM_exp_ccorr", # "_TCGA_NMD_targets_analysis_FPKM_exp_ccorr" "_TCGA_NMD_targets_analysis_FPKM_exp",
            "filter_mutations"           : False,
            "filter_mutation_margin"     : 0.25,
            "filter_mutations_type"      : "adjusted", # "adjusted" "threshold"
            "filter_samples"             : False,
            "hla_mutations_path"         : parent_dir+r"\data\tcga_hla_mutations.txt",
            "immune_data_path"           : parent_dir+r"\data\tcga_immune_scores.txt",
            "immune_editing_path"        : parent_dir+r"\data\tcga_immune_editing.txt",
            "info"                       : {
                                            "appended_features"     : ["Immune_score", "age_at_initial_pathologic_diagnosis", "ajcc_pathologic_tumor_stage",
                                                                       "HLAI", "HLAII", "CYT", "HED",
                                                                       "OS", "OS.time", "DFI", "DFI.time", "DSS", "DSS.time", "PFI", "PFI.time",
                                                                       "total_nmd", "frameshift_nmd", "missense_nmd", "nonsense_nmd",
                                                                       "total_hla", "frameshift_hla", "missense_hla", "nonsense_hla"],
                                            "extended_features"     : ["FEATURE:prediction", "ID:ptc reads", # "FEATURE:expression", 
                                                                       "ID:cnv total", "ID:cnv minor", "ID:cnv avg",
                                                                       "FEATURE:frameshift", "FEATURE:frameshift_mutations", "FEATURE:target", "FEATURE:escape",
                                                                       "FEATURE:ptc_mutations", "FEATURE:ptc_mutations2"],
                                            "ids"                   : ["project", "cluster_key", "variant_ids"],
                                            "combined_mutations"    : {"exon": ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                                                                                "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Silent"]},
                                            "mutations"             : ["total", "ptc_mutations", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "exon",
                                                                       "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "In_Frame_Del", "Nonstop_Mutation"],
                                            "rna"                   : ["fpkm_unstranded"], # from rnaseq data "fpkm_unstranded"
                                            "wxs"                   : ["Start_Position", "End_Position", "Tumor_Sample_Barcode", # from wxs data
                                                                       "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                                                                       "Variant_Classification", "HGVSc", "HGVSp", "HGVSp_Short"],
                                            "target_mutations"      : ["In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", # mutations in targets
                                                                       "Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "Nonstop_Mutation"]
                                            },
            "max_procs"                  : 4,
            "newdir"                     : None,
            "nmd_mutations_path"         : parent_dir+r"\data\tcga_nmd_mutations.txt",
            "non_targets_explicit"       : True, # if False non targets are defined as all entries that are not targets
            "os_sep"                     : "//",
            "processed_targets_fname"    : "TCGA_NMD_targets_analysis_FPKM_exp_ccorr.txt",
            "projects"                   : "all", # "all",
            "rna_calculation"            : ["avg", "avg", "avg"],
            "rna_threshold"              : 0,
            "run_dir"                    : parent_dir+r"\data",
            "sample_filter"              : {"FEATURE:frameshift": {"range": (0, 13), "inclusive": True, "include_misses": False}, "ptc_mutations": {"range": (50, 100000), "inclusive": True, "include_misses": False}},
            "score_threshold"            : {"FEATURE:escape": (0, 0.57), "FEATURE:target": (0.64, 1)},
            "shared_gene_fname"          : None,
            "stats_exceptions"           : {"FEATURE:escape": "relative_count_by_threshold", "FEATURE:target": "relative_count_by_threshold",
                                            "FEATURE:frameshift_mutations": "sum_wo_exception", "FEATURE:ptc_mutations2": "sum_wo_exception"},
            "status_fname"               : "tcga_data_info.json",
            "target_fname"               : "wt_nmd_targets.txt",
            "target_identifier"          : {"rna": "gene_name", "wxs": "Hugo_Symbol"},
            "target_normalization"       : True,
            "target_filter"              : 0,
            "target_size"                : None, # determined @ runtime, used for testing
            "target_threshold"           : 0, # expression of 0 is excluded
            "track_mutations"            : True,
            "transform_type"             : "RNA_ccorr_ptc", # "shared_RNA" used for SCLC-containing data "RNA_c_ptc" "RNA_ccorr_ptc" used for TCGA-only analysis
            "variant_value_filter"       : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0]}, # minimum expression for variants to be included in the analysis
            "wxs_identifier"             : "Transcript_ID" # wxs column that is used to create variant id and conduct ptc checks
        }

if platform.system() == "Windows": params["os_sep"] = "\\"
if params["track_mutations"] == False:
    params["info"]["mutations"] = []
    params["info"]["wxs"]       = []

if params["apply_filters"] == False:
    print("< warning. apply_filters is set to 'False'")

# marked (<-) added to facilitate skipping of integration of clinical and PTC data
if params["append_features"] == False: # <- starts here
    params["info"]["appended_features"] = []
    print("< appended_features set empty.")

if params["extend_features"] == False: 
    params["info"]["extended_features"] = []
    print("< extended_features set empty.") # <- ends here


def _extract(cluster, targets, non_targets):
    clusters = []
    for project_key in cluster:
        if "TCGA" in project_key or "SCLC" in project_key:
            for cluster_key in cluster[project_key]:
                if (project_key in params["projects"] or params["projects"] == "all"):
                    clusters.append({"project_key": project_key, "cluster_key": cluster_key})

    print("<", len(clusters), "selected.")

    rand_index = np.arange(len(clusters))
    np.random.shuffle(rand_index)
    clusters   = [{"project_key": clusters[rand_i]["project_key"], "cluster_key": clusters[rand_i]["cluster_key"]} for rand_i in rand_index]
    proc_index = split_index(len(clusters), params["max_procs"])
    lock  = Lock()
    procs = []

    for i in range(len(proc_index)):
        proc = Process(target=extract_targets, args=(cluster, clusters, targets, non_targets, proc_index[i], params, lock, i))
        procs.append(proc)
        proc.start()

    for i, proc in enumerate(procs):
        proc.join()

    return clusters


def main():
    start_time = time.time()

    if params["create_dataset"] == True:
        print("< calculation of", params["file_tag"].replace("_", " "))

        # load status data
        status  = load_status(params, fname=params["status_fname"])

        # initialize cluster dictionary
        cluster = init_cluster(status, params, cluster_key=params["cluster_key"])

        # filter cluster dictionary for case ids for which both RNA seq and WXS seq data exist
        if params["track_mutations"] == True:  cluster = filter_cluster(cluster, filter_mode="rna_and_wxs_present")
        if params["track_mutations"] == False: cluster = filter_cluster(cluster, filter_mode="rna_present")

        # load target file
        targets     = load_df(params["run_dir"]+params["os_sep"]+params["target_fname"], delimiter=",")
        targets     = targets.drop_duplicates(subset=[params["target_identifier"]["rna"]])
        non_targets = pd.DataFrame()

        if params["non_targets_explicit"] == True:
            non_targets = targets[targets["NMD_or_not"] == "N"]
            targets     = targets[targets["NMD_or_not"] == "Y"]

        if params["shared_gene_fname"] != None:
            shared_genes = pd.read_csv(params["run_dir"]+"\\"+params["shared_gene_fname"], delimiter=",")
            init_shape   = targets.shape[0]
            targets      = targets[targets["gene_name"].isin(shared_genes["gene"])]
            print("< removing unshared genes reduced dataset from", init_shape, "to", targets.shape[0])

        if params["shared_gene_fname"] == None:
            rna          = load_rna(params["run_dir"]+params["os_sep"]+"tcga_rna_example", params)
            init_shape   = targets.shape[0]
            targets      = targets[targets["gene_name"].isin(rna["gene_name"])]
            print("< removing unshared genes reduced dataset from", init_shape, "to", targets.shape[0])

        params["target_size"] = targets.shape[0]

        # extract expression data
        clusters = _extract(cluster, targets, non_targets)
        print("< dataset calculation time:", round(time.time()-start_time, 4))

        # loading from storage for aggregation
        cluster = aggregate(cluster, clusters, params, delimiter=",", folder="rnaseq"+params["file_tag"], path=params["run_dir"]+params["os_sep"]+params["processed_targets_fname"])

        # save parameters
        params_json = json.dumps(params, indent=4)
        with open(params["run_dir"]+params["os_sep"]+params["processed_targets_fname"].split(".")[0]+"_params.json", "w") as file:
            file.write(params_json)


    if params["create_calculations"] == True:
        print("< analysis of", params["file_tag"].replace("_", " "))

        # initialize folder for analysis run
        create_newdir(params)
        params_json = json.dumps(params, indent=4)
        with open(params["newdir"]+params["os_sep"]+"params.json", "w") as file:
            file.write(params_json)

        # load dataframe
        targets = load_df(params["run_dir"]+params["os_sep"]+params["processed_targets_fname"], delimiter=",")
        extended_features = pd.DataFrame()

        # marked (<-) added/removed to facilitate skipping of integration of clinical and PTC data
        #if len(params["info"]["extended_features"]) > 0: # <- removed
        if params["extend_features"] == True and len(params["info"]["extended_features"]) > 0: # <- added
            extended_features                           = load_df(params["extended_features_path"], delimiter=",")
            extended_features["FEATURE:ptc_mutations2"] = [1 for _ in range(extended_features.shape[0])]
            extended_features["FEATURE:escape"]         = extended_features["FEATURE:prediction"]
            extended_features["FEATURE:target"]         = extended_features["FEATURE:prediction"]

            if params["apply_filters"] == True:
                # filter out Y-chromosomes and X-chromosomes for men
                init_shape        = extended_features.shape[0]
                extended_features = extended_features[extended_features["ID:chromosome"] != "chrY"]
                extended_features = extended_features[[True if extended_features.iloc[j].loc["ID:chromosome"] != "chrX" or (extended_features.iloc[j].loc["ID:chromosome"] == "chrX"
                                                       and extended_features.iloc[j].loc["ID:gender"] == "FEMALE") else False for j in range(extended_features.shape[0])]]
                print("< filtering reduced dataset from", init_shape, "to", extended_features.shape[0])

                # apply further filters defined by variant_value_filter
                if len(params["variant_value_filter"]) > 0:
                    # changed 250313 (now same as in analyze_predictions)
                    #extended_features = filter_variants(extended_features, params)
                    extended_features = apply_value_filter(extended_features, params["variant_value_filter"], "variant filter")

                # check for expired placeholders
                check_placeholders(extended_features, params["info"]["extended_features"])
            
            
            counts = []
            for extended_feature in params["info"]["extended_features"]:
                counts.append(extended_feature + "_n")

            params["info"]["extended_features"].extend(counts)

        targets = remove_misses(targets, params)
        # marked (<-) added/removed on 250429 to allow skipping of append_features
        # targets = append_features(targets, params) # <- removed
        if params["append_features"] == True:
            targets = append_features(targets, params) # <- added

        # conduct calculations
        targets = calculate_target_expression(targets, extended_features, params)


if __name__ == '__main__':
    main()