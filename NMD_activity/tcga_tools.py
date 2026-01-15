import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
from scipy.cluster.hierarchy import linkage, dendrogram
import shutil
import sys

sys.path.insert(0, r"C:\Programming\Translational_genomics\NMD_analysis\shared")
from extract_mutations_utils import *
from tcga_tools_utils import *


clinical_data_path     = r"C:\Programming\Translational_genomics\NMD_analysis\data\TCGA-CDR-SupplementalTableS1_extracted_appended.txt"
data_dir               = r"C:\Programming\Translational_genomics\NMD_analysis\data"
fname                  = "filtered_status_cpatc3.json" # "cancer_scores_TCGA_SCLC_NMD_targets_analysis_FPKM_full_mutations_exp"
fnames                 = ["cptac3.txt", "cptac3_clinical.tsv"] #["tcga_primary.txt", "TCGA-CDR-SupplementalTableS1_extracted_appended.txt"] #["tcga.txt", "TCGA-CDR-SupplementalTableS1_extracted_appended.txt"]
immune_data_path       = r"C:\Programming\Translational_genomics\NMD_analysis\data\TCGA_ESTIMATE_appended.txt"
# "append_clinical_data" "assemble_test_results"
# "check_cnvs" "check_files" "check_sample_types" "check_sample_ids" "compare_status_files"
# "filter_sample_types" "filter_tcga_genes"
# "map_types"
# "split_data" "test_expressions" "test_expressions2" "test_ptc_mutations"
mode                   = "append_clinical_data"
status_path            = r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity\status_w_md5sum.json"
target_dir             = r"C:\Programming\Translational_genomics\TCGA"


#  used to create TCGA-CDR-SupplementalTableS1_extracted_appended.txt (appending case ids based on submitter id)
if mode == "append_case_ids":
    id_col = "bcr_patient_barcode" # clinical data (TCGA-CDR-SupplementalTableS1_extracted.csv)

    df = pd.read_csv(data_dir+"\\"+fnames[0], delimiter=";")

    with open(data_dir+"\\"+fnames[1], "r") as file:
        status = json.load(file)

    case_map = create_case_map(status)

    case_map["submitter_id_"] = [case_map.iloc[i].loc["submitter_id"].replace("-", "") for i in range(case_map.shape[0])]
    df["submitter_id_"]       = [df.iloc[i].loc[id_col].replace("-", "")[0:10] for i in range(df.shape[0])]

    case_ids = append_df_with_mapping([df, case_map], "submitter_id_", "submitter_id_", "case_id", "mapping submitter ids")
    df.insert(0, "case id", case_ids)
    df = df.drop(columns=["submitter_id_"])
    df.to_csv(data_dir+"\\"+fname.split(".")[0]+"_appended.txt", sep=",")


# located here as test tool for mode above:
if mode == "check_case_id_mapping":
    df = pd.read_csv(data_dir+"\\"+fname, delimiter=",")

    with open(status_path, "r") as file:
        status = json.load(file)

    bar = IncrementalBar("checking case id mapping", max=df.shape[0])
    for i in range(df.shape[0]):
        project = "TCGA-"+df.iloc[i].loc["type"]
        barcode = df.iloc[i].loc["bcr_patient_barcode"]
        case_id = df.iloc[i].loc["case id"]

        if len(status["case_ids"][project]["case_id"]) != len(status["case_ids"][project]["submitter_id"]):
            print("< case and submitter id sizes inconsistent for", project, len(status["case_ids"][project]["case_id"]), "/", len(status["case_ids"][project]["submitter_id"]))

        elif barcode in status["case_ids"][project]["submitter_id"]:
            pos = status["case_ids"][project]["submitter_id"].index(barcode)

            if pos >= 0 and pos < len(status["case_ids"][project]["submitter_id"]):
                if status["case_ids"][project]["case_id"][pos] != case_id:
                    print("< mismatch @", project, i, "case id", case_id)

            else:
                print("< exceeded position @", project, i, "case id", case_id)

        else:
            print("< barcode", barcode, "not found @", project, i)


# clinical data file (TCGA): TCGA-CDR-SupplementalTableS1_extracted_appended.txt
# clinical data file (CPTAC3): cptac3_clinical.tsv
# gender data needed for exclusion of monoallelic variants (male X)
if mode == "append_clinical_data":
    datatype  = "CPTAC3" # CPTAC3, TCGA

    if datatype == "CPTAC3":
        col_names  = {"LABEL:demographic.gender": "LABEL:gender"} # maps CPTAC3 names to TCGA names
        delimiters = [",", "\t"]
        features   = ["demographic.gender", "cases.disease_type"]
        target_col = "cases.case_id"

    if datatype == "TCGA":
        delimiters = [",", ","]
        features   = ["gender"] # ["age_at_initial_pathologic_diagnosis", "ajcc_pathologic_tumor_stage", "gender", "OS", "OS.time", "DFI", "DFI.time", "DSS", "DSS.time", "PFI", "PFI.time"]
        target_col = "case id"

    data = []
    for i, fname in enumerate(fnames):
        data.append(pd.read_csv(data_dir+"\\"+fname, delimiter=delimiters[i]))

        
    data[0] = map_data(data, features, ["ID:case id", target_col], "appending clinical data", tag="LABEL")
    if datatype == "CPTAC3":
        data[0] = data[0].rename(columns=col_names)
        # for CPTAC3 only, project (CPTAC3) is replaced with cancer types
        cancer_types          = np.unique(data[0]["LABEL:cases.disease_type"])
        data[0]["ID:project"] = data[0]["LABEL:cases.disease_type"]
        data[0]               = data[0].drop(columns=["LABEL:cases.disease_type"])
        print(json.dumps({cancer_type: data[0][data[0]["ID:project"] == cancer_type].shape[0]
                          for cancer_type in cancer_types}))
        
    data[0].to_csv(data_dir+"\\"+fnames[0].split(".")[0]+"_appended.txt", sep=",")


# appending TCGA_ESTIMATE_appended.txt
if mode == "append_features":
    datatype = "ESTIMATE"

    if datatype == "TCGA":
        features    = ["Variant_Classification", "RNASEQ_ptc_cnv", "RNASEQ_noptc_cnv"]
        delimiters  = ["," ,"\t"]
        mapping_col1 = "temp_id"
        mapping_col2 = "temp_id"

    if datatype == "ESTIMATE":
        features    = ["Stromal_score", "Immune_score", "ESTIMATE_score"]
        delimiters  = ["," ,","]
        mapping_col1 = "ID:case id"
        mapping_col2 = "case id"
    
    data = []

    for i in range(len(fnames)):
        data.append(pd.read_csv(data_dir+"\\"+fnames[i], delimiter=delimiters[i]))


    data[1]["index"] = [i for i in range(data[1].shape[0])]

    if datatype == "TCGA":
        data[0]["temp_id"] = [data[0].iloc[i].loc["ID:variant id"] + data[0].iloc[i].loc["ID:case id"] for i in range(data[0].shape[0])]
        data[1]["temp_id"] = [data[1].iloc[i].loc["variant_id"] + data[1].iloc[i].loc["case_id"] for i in range(data[1].shape[0])]

    test = map_data(data, features, [mapping_col1, mapping_col2], "appending immune data", tag="LABEL")
    data[0].to_csv(data_dir+"\\"+fnames[0].split(".")[0]+"_appended.txt", sep=",")


if mode == "assemble_test_results":
    correlation_targets = ["FEATURE:ptc_mutations2_FEATURE:cnv total", "FEATURE:ptc_mutations2_fpkm_unstranded", "FEATURE:ptc_mutations2_FEATURE:escape",
                           "FEATURE:ptc_mutations2_FEATURE:target", "FEATURE:ptc_mutations2_FEATURE:prediction", "FEATURE:ptc_mutations2_total_hla"]
    
    class_targets = ["ptc_mutations2", *[correlation_target.split("_")[-1] for correlation_target in correlation_targets]]
    
    fnames = os.listdir(data_dir)
    correlations = []
    for fname in fnames:
        if "stats_summary_overview" in fname:
            data = pd.read_csv(data_dir+"\\"+fname, delimiter=",")
            data = data[data["pair"].isin(correlation_targets)]
            data.insert(0, "project", [fname.split("_")[-1] for _ in range(data.shape[0])])
            correlations.append(data)

    correlations = pd.concat(correlations)
    pd.set_option('display.max_rows', None)
    print(correlations[correlations["spearman-p"] <= 0.05])


# comparison of cnv correction outputs from old and new version of hg38 build
if mode == "check_cnvs":
    projects = ["TCGA-ACC"]

    deviating_gene_ids = []

    for project in projects:
        fnames = os.listdir(target_dir+"\\"+project+"\\RNA_ccorr_ptc")

        for fname in fnames:
            data1 = load_rna("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\RNA_ccorr_ptc\\"+fname, {"transform_type": None})
            data2 = load_rna("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\RNA_ccorr_ptc_test\\"+fname, {"transform_type": None})
            data1.index = data1["gene_id"]
            data2.index = data2["gene_id"]

            deviations = data1[data1["fpkm_unstranded"] != data2["fpkm_unstranded"]]
            deviating_gene_ids.extend(deviations["gene_id"].tolist())

            if deviations.shape[0] > 0:
                print("<", project, fname)
                for i in range(deviations.shape[0]):
                    print(i, deviations.iloc[i].loc["gene_id"], deviations.iloc[i].loc["fpkm_unstranded"], "/", data2.loc[deviations.iloc[i].loc["gene_id"]].loc["fpkm_unstranded"],
                        deviations.iloc[i].loc["cnv_total"], "/", data2.loc[deviations.iloc[i].loc["gene_id"]].loc["cnv_total"])

    print("< deviating gene ids:")
    deviating_gene_ids = np.unique(deviating_gene_ids)
    for i, deviating_gene_id in enumerate(deviating_gene_ids):
        print(i, deviating_gene_id)


if mode == "check_files":
    subdirs = os.listdir(data_dir)
    for subdir in subdirs:
        if os.path.isdir(data_dir+"\\"+subdir) == True:
            subsubdirs = os.listdir(data_dir+"\\"+subdir)
            
            for subsubdir in subsubdirs:
                if os.path.isdir(data_dir+"\\"+subdir+"\\"+subsubdir) == True:
                    fnames = os.listdir(data_dir+"\\"+subdir+"\\"+subsubdir)
                    #if subsubdir == "RNA":            print("<", data_dir+"C"+subdir+"\\RNA", len(fnames))
                    #if  subsubdir == "RNA_c":         print("<", data_dir+"\\"+subdir+"\\RNA_c", len(fnames))
                    #if  subsubdir == "RNA_ccorr_ptc": print("<", data_dir+"\\"+subdir+"\\RNA_ccorr_ptc", len(fnames))
                    if  subsubdir == "RNA_c_ptc":     print("<", data_dir+"\\"+subdir+"\\RNA_c_ptc", len(fnames))
                    #if  subsubdir == "shared_RNA_c":  print("<", data_dir+"\\"+subdir+"\\shared_RNA_c", len(fnames))
                    #if  subsubdir == "shared_RNA_cq": print("<", data_dir+"\\"+subdir+"\\shared_RNA_cq", len(fnames))


# check for sample ids (downloaded on 250807)
if mode == "check_sample_ids":
    with open(data_dir+"\\"+fname, "r") as file:
        status = json.load(file)

    target_cols = ["CNV_ranges", "RNA", "WXS"]
    sample_stats = {"CNV_ranges_multiple": 0, "RNA_multiple": 0, "WXS_multiple": 0, "CNV_ranges_total": 0, "RNA_total": 0, "WXS_total": 0}

    total_ids = {}
    bar = IncrementalBar(set_bar("checking sample ids"), max=np.sum([len(status["file_ids"][project]) for project in status["file_ids"]]))
    for project in status["file_ids"]:
        for case_id in status["file_ids"][project]:
            current_case = status["file_ids"][project][case_id]
            #print(project, case_id)
            #print(json.dumps(current_case, indent=4))
            for target_col in target_cols:
                current_ids = {}
                for i in range(len(current_case[target_col+"_sample_type"])):
                    for j in range(len(current_case[target_col+"_sample_type"][i])):
                        #print("i", i, "j", j, target_col, current_case[target_col+"_sample_type"][i], current_case[target_col+"_sample_type"][i][j], len(current_case[target_col+"_sample_type"][i]), len(current_case[target_col+"_submitter_id"][i]))
                        if current_case[target_col+"_sample_type"][i][j] in current_ids: current_ids[current_case[target_col+"_sample_type"][i][j]].append(current_case[target_col+"_submitter_id"][i][j].split("-")[-1])
                        else:                                                            current_ids[current_case[target_col+"_sample_type"][i][j]] = [current_case[target_col+"_submitter_id"][i][j].split("-")[-1]]
                        if current_case[target_col+"_sample_type"][i][j] in total_ids:   total_ids[current_case[target_col+"_sample_type"][i][j]].append(current_case[target_col+"_submitter_id"][i][j].split("-")[-1])
                        else:                                                            total_ids[current_case[target_col+"_sample_type"][i][j]] = [current_case[target_col+"_submitter_id"][i][j].split("-")[-1]]
            
                for key in current_ids:
                    if len(np.unique(current_ids[key])) > 1:
                        print(project, case_id, key, current_ids[key])
                        sample_stats[target_col+"_multiple"] += 1
                    
                    sample_stats[target_col+"_total"] += 1

                #print(target_col)
                #print(json.dumps(current_ids, indent=4))
                
            bar.next()
    bar.finish()
    print(json.dumps(sample_stats, indent=4))


# comparison of unfiltered and filtered status
if mode == "check_sample_types":
    target_folders = ["RNA_c_ptc"] #, "shared_RNA", "shared_RNA_c", "shared_RNA_cq", "shared_RNA_q"]

    with open(data_dir+"\\"+fnames[0], "r") as file:
        unfiltered_status = json.load(file)

    with open(data_dir+"\\"+fnames[1], "r") as file:
        filtered_status   = json.load(file)

    for project in unfiltered_status["file_ids"]:
        if "TCGA" in project:
            paths           = os.listdir(target_dir+"\\"+project+"\\RNA")
            unfiltered_size = len(paths)
            
            unfiltered_status_count = 0
            for case_id in unfiltered_status["file_ids"][project]:
                if len(unfiltered_status["file_ids"][project][case_id]["WXS"]) > 0:
                    unfiltered_status_count += len(unfiltered_status["file_ids"][project][case_id]["RNA"])

            filtered_status_count = 0
            filtered_status_paths = []
            for case_id in filtered_status["file_ids"][project]:
                if len(filtered_status["file_ids"][project][case_id]["WXS"]) > 0:
                    filtered_status_count += len(filtered_status["file_ids"][project][case_id]["RNA"])
                    filtered_status_paths.extend(filtered_status["file_ids"][project][case_id]["RNA"])

            for target_folder in target_folders:
                filtered_paths = os.listdir(target_dir+"\\"+project+"\\"+target_folder)
                filtered_size  = len(filtered_paths)
                # test for different counts of existing files and files deposited in unfiltered status file
                gap1           = unfiltered_status_count-unfiltered_size
                # test for difference of count differences of files deposited in unfiltered and filtered status file and of existing unfiltered or filtered (modified) files
                gap2           = (unfiltered_status_count-filtered_status_count)-(unfiltered_size-filtered_size)

                #print(project, target_folder, "status size", filtered_status_count, "/", unfiltered_status_count, "file count", filtered_size, "/", unfiltered_size, "gap", gap1, "/", gap2)
                if gap2 != 0:
                    for filtered_status_path in filtered_status_paths:
                        if filtered_status_path not in filtered_paths:
                            print("  ", filtered_status_path)


# compare status files with respect to all levels and report differences
if mode == "compare_status_files":
    with open(data_dir+"\\"+fnames[0], "r") as file:
        status1 = json.load(file)

    with open(data_dir+"\\"+fnames[1], "r") as file:
        status2   = json.load(file)

    mismatches = compare_status(status1, status2)
    selected_mismatches = [entry for entry in mismatches["status1"] if "case_ids" in entry and "submitter_id" not in entry]
    misses = {"cases": len(selected_mismatches), "cnv": 0, "rna": 0, "wxs": 0, "all": 0}
    complete_misses = []; missing_cases = []
    for i, entry in enumerate(selected_mismatches):
        missing_cases.append(entry.split("_")[-1])
        case = status1["file_ids"][entry.split("_")[2]][entry.split("_")[-1]]
        if len(case["CNV_ranges"]) > 0: misses["cnv"] +=1
        if len(case["RNA"]) > 0:        misses["rna"] +=1
        if len(case["WXS"]) > 0:        misses["wxs"] +=1

        if len(case["CNV_ranges"]) > 0 and len(case["RNA"]) > 0 and len(case["WXS"]) > 0:
            complete_misses.append(entry.split("_")[-1])
            misses["all"] += 1

    print("< complete missing cases (represented in clinical data file)")
    clinical_data = pd.read_csv(clinical_data_path, delimiter=",")
    complete_misses = [complete_miss for complete_miss in complete_misses if complete_miss in clinical_data["case id"]]
    for i, complete_miss in enumerate(complete_misses):
        print(i, complete_miss)
    
    selected_mismatches = [entry for entry in mismatches["status1"] if "file_ids" in entry and entry.split("_")[3] not in missing_cases and "CNV_ranges" not in entry]
    selected_cases      = [entry.split("_")[3] for entry in selected_mismatches]
    selected_projects   = [entry.split("_")[2] for entry in selected_mismatches]

    print("< cases with missing expression data with no expression data at all")
    for i, selected_case in enumerate(selected_cases):
        if len(status2["file_ids"][selected_projects[i]][selected_case]["RNA"]) == 0:
            print(i, selected_case)


if mode == "convert_cancer_scores":
    datatype    = "raw"
    ids         = ["sample_names"]

    cancer_scores = pd.read_csv(data_dir+"\\"+fname, delimiter=",")

    features = [
               *["fpkm_unstranded",
                 "ID:noptc reads",
                 "ptc_mutations"],
               *[col.replace("_mean", "") for col in cancer_scores.columns if "FEATURE" in col and "mean" in col and col.split("_")[len(col.split("_"))-2] != "n"]
               ] 
    
    cancer_scores = convert_cancer_scores(cancer_scores, feature_targets=features, id_targets=ids, label_targets=["FEATURE:OS.time"], data_type=datatype, randomize=True)
    cancer_scores.to_csv(data_dir+"\\"+fname, sep=",", index=False)


if mode == "correlate_cancer_scores":
    data_dir = r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity\2025-06-12_17-29-04_TCGA_SCLC_NMD_targets_analysis_FPKM_full_mutations_exp_rerun"
    cancer_scores = pd.read_csv(data_dir+ "\\transformed_cancer_scores", delimiter=",")
    mutations     = pd.read_csv(data_dir+ "\\transformed_exon_mutations", delimiter=",")
    cancer_scores.index = cancer_scores["project"]
    mutations.index     = mutations["project"]
    projects            = np.unique(cancer_scores["project"])
    print(cancer_scores)
    rearranged_scores = pd.DataFrame({**{project+"_cancer_score": [cancer_scores.loc[project].loc[col] if pd.isna(cancer_scores.loc[project].loc[col]) == False else "NA"
                                                                   for col in cancer_scores.columns if "sample score" in col] for project in projects},
                                      **{project+"_mutations": [mutations.loc[project].loc[col] if pd.isna(mutations.loc[project].loc[col]) == False else "NA"
                                                                for col in mutations.columns if "sample score" in col] for project in projects}})
    rearranged_scores.to_csv(data_dir+"\\cancer_scores_exon_mutations.txt", sep=",", index=False)

    print(rearranged_scores)
    correlation_stats = pd.DataFrame({"r":   ["NA" for _ in [*projects, "total"]],
                                      "r-p": ["NA" for _ in [*projects, "total"]],
                                      "s":   ["NA" for _ in [*projects, "total"]],
                                      "s-p": ["NA" for _ in [*projects, "total"]]}, index=[*projects, "total"])
    
    total_cancer_scores = []; total_mutations = []
    for project in projects:
        cancer_scores = [float(cancer_score) for cancer_score in rearranged_scores[project+"_cancer_score"] if cancer_score != "NA"]
        mutations     = [float(mutation) for mutation in rearranged_scores[project+"_mutations"] if mutation != "NA"]
        total_cancer_scores.extend(cancer_scores)
        total_mutations.extend(mutations)
        print(project, len(cancer_scores), len(mutations))
        correlation_stats.at[project, "r"]   = scipy.stats.pearsonr(cancer_scores, mutations).statistic
        correlation_stats.at[project, "r-p"] = scipy.stats.pearsonr(cancer_scores, mutations).pvalue
        correlation_stats.at[project, "s"]   = scipy.stats.spearmanr(cancer_scores, mutations).statistic
        correlation_stats.at[project, "s-p"] = scipy.stats.spearmanr(cancer_scores, mutations).pvalue

    correlation_stats.at["total", "r"]   = scipy.stats.pearsonr(total_cancer_scores, total_mutations).statistic
    correlation_stats.at["total", "r-p"] = scipy.stats.pearsonr(total_cancer_scores, total_mutations).pvalue
    correlation_stats.at["total", "s"]   = scipy.stats.spearmanr(total_cancer_scores, total_mutations).statistic
    correlation_stats.at["total", "s-p"] = scipy.stats.spearmanr(total_cancer_scores, total_mutations).pvalue
    print(correlation_stats)
    correlation_stats.to_csv(data_dir+"\\cancer_scores_exon_mutations_correlation_stats.txt", sep=",")


# determine intersection of TCGA and SCLC genes (["log_normalised_counts.csv", "TCGA_rna_example"])
if mode == "create_tcga_cuomo_intersection":
    cuomo_rna = pd.read_csv(data_dir+"\\"+fnames[0], delimiter=",")
    tcga_rna  = pd.read_csv(data_dir+"\\"+fnames[1], delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5])
    
    cuomo_rna["gene_name"] = [cuomo_rna.iloc[i].loc[cuomo_rna.columns.tolist()[0]].split("_")[1] for i in range(cuomo_rna.shape[0])]
    selected_cuomo_rna     = pd.DataFrame({"gene": cuomo_rna[cuomo_rna["gene_name"].isin(tcga_rna["gene_name"])]["gene_name"]})
    selected_tcga_rna      = pd.DataFrame({"gene": tcga_rna[tcga_rna["gene_name"].isin(cuomo_rna["gene_name"])]["gene_name"]})
    shared_rna             = pd.concat([selected_cuomo_rna, selected_tcga_rna])
    shared_rna_genes       = shared_rna.drop_duplicates(subset=["gene"])
    print("sh", cuomo_rna.shape, tcga_rna.shape, selected_cuomo_rna.shape, selected_tcga_rna.shape, shared_rna.shape, shared_rna_genes.shape)
    shared_rna_genes.to_csv(data_dir+"\\Cuomo_TCGA_shared_genes.txt", sep=",", index=False)


# determine intersection of TCGA and SCLC genes (["NIHMS782739-supplement-Tables_S10.txt", "TCGA_rna_example"])
if mode == "create_tcga_sclc_intersection":
    sclc_rna = pd.read_csv(data_dir+"\\"+fnames[0], delimiter=";")
    tcga_rna = pd.read_csv(data_dir+"\\"+fnames[1], delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5])
    
    selected_sclc_rna = pd.DataFrame({"gene": sclc_rna[sclc_rna["gene"].isin(tcga_rna["gene_name"])]["gene"]})
    selected_tcga_rna = pd.DataFrame({"gene": tcga_rna[tcga_rna["gene_name"].isin(sclc_rna["gene"])]["gene_name"]})
    shared_rna        = pd.concat([selected_sclc_rna, selected_tcga_rna])
    shared_rna_genes  = shared_rna.drop_duplicates(subset=["gene"])
    print("sh", sclc_rna.shape, tcga_rna.shape, selected_sclc_rna.shape, selected_tcga_rna.shape, shared_rna.shape, shared_rna_genes.shape)

    shared_rna_genes.to_csv(data_dir+"\\SCLC_TCGA_shared_genes.txt", sep=",", index=False)


if mode == "filter_sample_types":
    extraction_size    = 50
    extraction_targets = [] # ["Solid Tissue Normal"] # define categories to be extracted to novel project
    # define categories to be preserved in the updated status
    #filter_targets     = ["Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"] # used for NMD activity calculation
    filter_targets     = ["Additional - New Primary", "Metastatic", "Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"] # used for TCGA PTC extraction
    print_extraction   = False
    print_status       = True
    tag                = "test"

    with open(data_dir+"\\"+fname, "r") as file:
        status = json.load(file)

    file_ids = {project: {case_id: {"CNV_genes": [], "CNV_ranges": [], "RNA": [], "WXS": [],
                                    "CNV_genes_sample_type": [], "CNV_ranges_sample_type": [],
                                    "RNA_sample_type": [], "WXS_sample_type": []}
                for case_id in status["file_ids"][project]} for project in status["file_ids"]}

    sample_type_stats = {project: {"CNV_genes": {}, "CNV_ranges": {}, "RNA": {}, "WXS": {}} for project in status["file_ids"]}
    filtering_stats   = {project: {"CNV_genes": {}, "CNV_ranges": {}, "RNA": {}, "WXS": {}} for project in status["file_ids"]}
    project_stats     = {}

    extracted_status  = {"case_ids": {}, "file_ids": {}} # container to transiently store extracted files (later checked for overall count)
    for project in status["file_ids"]:
        for case_id in status["file_ids"][project]:
            for key in sample_type_stats[project]:
                for i in range(len(status["file_ids"][project][case_id][key])):                  
                    sample_type = ""
                    #print(status["file_ids"][project][case_id])
                    #print(i, status["file_ids"][project][case_id][key+"_sample_type"])
                    
                    for j in range(len(status["file_ids"][project][case_id][key+"_sample_type"][i])):
                        if j < len(status["file_ids"][project][case_id][key+"_sample_type"][i])-1:
                            sample_type += status["file_ids"][project][case_id][key+"_sample_type"][i][j] + ", "

                        else:
                            sample_type += status["file_ids"][project][case_id][key+"_sample_type"][i][j]

                        #print(project, case_id, key, i, j, sample_type, status["file_ids"][project][case_id][key+"_sample_type"][i][j])

                    if sample_type in list(sample_type_stats[project][key].keys()):
                        sample_type_stats[project][key][sample_type].append(status["file_ids"][project][case_id][key][i])

                    else:
                        sample_type_stats[project][key][sample_type] = [status["file_ids"][project][case_id][key][i]]

                    filter_passed = False
                    for filter_target in filter_targets:
                        if filter_target in sample_type: filter_passed = True

                    if filter_passed == True:
                        file_ids[project][case_id][key].append(status["file_ids"][project][case_id][key][i])
                        file_ids[project][case_id][key+"_sample_type"].append(status["file_ids"][project][case_id][key+"_sample_type"][i])

                        if sample_type in list(filtering_stats[project][key].keys()):
                            filtering_stats[project][key][sample_type].append(status["file_ids"][project][case_id][key][i])

                        else:
                            filtering_stats[project][key][sample_type] = [status["file_ids"][project][case_id][key][i]]

                    if len(extraction_targets) > 0 and key == "RNA":
                        extraction_passed = False
                        for extraction_target in extraction_targets:
                            if extraction_target in sample_type: extraction_passed = True; matched_extraction_target = extraction_target

                        if extraction_passed == True:
                            extracted_project = project + "_" + matched_extraction_target.replace(" ", "_")
                            if extracted_project not in extracted_status["case_ids"]:
                                extracted_status["case_ids"][extracted_project] = {"case_id": [], "submitter_id": [], "tumor_code": [], "tumor_code_id": []}
                                extracted_status["file_ids"][extracted_project] = {}

                            if case_id not in extracted_status["case_ids"][extracted_project]["case_id"]:
                                extracted_status["case_ids"][extracted_project]["case_id"].append(case_id)
                                extracted_status["case_ids"][extracted_project]["submitter_id"].append(case_id)
                                extracted_status["case_ids"][extracted_project]["tumor_code"].append(extracted_project)
                                extracted_status["case_ids"][extracted_project]["tumor_code_id"].append(extracted_project)

                            if case_id in extracted_status["file_ids"][extracted_project]:
                                extracted_status["file_ids"][extracted_project][case_id]["RNA"].append(status["file_ids"][project][case_id][key][i])
                                extracted_status["file_ids"][extracted_project][case_id]["RNA_sample_type"].append(status["file_ids"][project][case_id][key+"_sample_type"][i])

                            else:
                                extracted_status["file_ids"][extracted_project][case_id] = {"CNV_genes": [], "CNV_ranges": [], "RNA": [status["file_ids"][project][case_id][key][i]], "WXS": [],
                                                                                            "CNV_genes_sample_type": [], "CNV_ranges_sample_type": [],
                                                                                            "RNA_sample_type": [status["file_ids"][project][case_id][key+"_sample_type"][i]], "WXS_sample_type": []}


    # if extraction targets are present, apply file size filter and create folders
    if print_extraction == True and len(extracted_status["case_ids"]) > 0:
        for extracted_project in extracted_status["case_ids"]:
            if len(extracted_status["case_ids"][extracted_project]["case_id"]) >= extraction_size:
                if os.path.isdir(target_dir+"\\"+extracted_project) == False:
                    os.mkdir(target_dir+"\\"+extracted_project)

                if os.path.isdir(target_dir+"\\"+extracted_project+"\\RNA") == False:
                    os.mkdir(target_dir+"\\"+extracted_project+"\\RNA")

                for case_id in extracted_status["file_ids"][extracted_project]:
                    for i in range(len(extracted_status["file_ids"][extracted_project][case_id]["RNA"])):
                        #print("from", target_dir+"\\"+extracted_project.split("_")[0]+"\\RNA\\"+extracted_status["file_ids"][extracted_project][case_id]["RNA"][i],
                        #      "to", target_dir+"\\"+extracted_project+"\\RNA"+extracted_status["file_ids"][extracted_project][case_id]["RNA"][i])
                        if os.path.isfile(target_dir+"\\"+extracted_project.split("_")[0]+"\\RNA\\"+extracted_status["file_ids"][extracted_project][case_id]["RNA"][i]) == True:
                            shutil.copyfile(target_dir+"\\"+extracted_project.split("_")[0]+"\\RNA\\"+extracted_status["file_ids"][extracted_project][case_id]["RNA"][i],
                                            target_dir+"\\"+extracted_project+"\\RNA\\"+extracted_status["file_ids"][extracted_project][case_id]["RNA"][i])

    for project in sample_type_stats:
        for key in sample_type_stats[project]:
            for sample_type in sample_type_stats[project][key]:
                if key == "RNA":
                    if sample_type in filtering_stats[project][key]: print(project, key, sample_type, len(sample_type_stats[project][key][sample_type]), "/", len(filtering_stats[project][key][sample_type]))
                    else:                                            print(project, key, sample_type, len(sample_type_stats[project][key][sample_type]), "/ 0")
                
    if print_status == True and len(filter_targets) > 0:
        status["file_ids"] = file_ids
        
        if len(extracted_status["case_ids"]) > 0:
            for extracted_project in extracted_status["case_ids"]:
                if len(extracted_status["case_ids"][extracted_project]["case_id"]) >= extraction_size:
                    status["case_ids"][extracted_project] = extracted_status["case_ids"][extracted_project]
                    status["file_ids"][extracted_project] = extracted_status["file_ids"][extracted_project]
        
        status = json.dumps(status, indent=4)
        if len(filter_targets) > 0 and len(extraction_targets) == 0: fname_tag = "filtered_"
        if len(filter_targets) == 0 and len(extraction_targets) > 0: fname_tag = "extracted_"
        if len(filter_targets) > 0 and len(extraction_targets) > 0:  fname_tag = "filtered_extracted_"

        with open(data_dir+"\\"+tag+"_"+fname_tag+fname, "w") as file:
            file.write(status)


# files: ["status.json", "Cuomo_TCGA_shared_genes.txt"]
# files: ["status.json", "SCLC_TCGA_shared_genes.txt"]
if mode == "filter_tcga_genes":
    source_folder = "RNA_c"
    target_folder = "shared_RNA_c"
    print_files   = True
    # load shared genes file for filtering
    with open(data_dir+"\\"+fnames[1], "r") as file:
        shared_genes = pd.read_csv(data_dir+"\\"+fnames[1], delimiter=",")

    with open(target_dir+"\\"+fnames[0], "r") as file:
        status = json.load(file)

    params  = {"data_dir": data_dir, "os_sep": "\\"}
    cluster = init_cluster(status, params, cluster_key="Firehouse")

    init      = True
    showcount = 0
    duplicate_genes = [] # list of gene names / symbols that are ambiguous
    for project in cluster:
        if "TCGA" in project:# or "TCGA-H" in project:
            if os.path.isdir(target_dir+"\\"+project+"\\"+target_folder) == False:
                os.mkdir(target_dir+"\\"+project+"\\"+target_folder)

            print(project)
            for cluster_key in cluster[project]:
                for i in range(len(cluster[project][cluster_key]["input"])):
                    for case_id in cluster[project][cluster_key]["input"][i]:
                        for rna_fname in cluster[project][cluster_key]["input"][i][case_id]["RNA"]:
                            if os.path.isfile(target_dir+"\\"+project+"\\"+source_folder+"\\"+rna_fname) == True:
                                tcga_rna = pd.read_csv(target_dir+"\\"+project+"\\"+source_folder+"\\"+rna_fname, delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5])
                                init_sum = tcga_rna["fpkm_unstranded"].sum()
                                tcga_rna = tcga_rna[tcga_rna["gene_name"].isin(shared_genes["gene"])].sort_values(by="gene_name")
                                tcga_rna = pd.DataFrame({"gene_name": tcga_rna["gene_name"], "fpkm_unstranded": tcga_rna["fpkm_unstranded"]})
                                selected_tcga_rna       = tcga_rna.drop_duplicates(subset="gene_name")
                                selected_tcga_rna.index = selected_tcga_rna["gene_name"]

                                if init == True:
                                    duplicated_genes = tcga_rna[tcga_rna.duplicated(subset="gene_name")].drop_duplicates(subset="gene_name")["gene_name"].tolist()
                                    init = False

                                for duplicated_gene in duplicated_genes:
                                    duplicated_rna = tcga_rna[tcga_rna["gene_name"] == duplicated_gene]
                                    duplicated_rna = duplicated_rna[duplicated_rna["fpkm_unstranded"] > 0]
                                    if duplicated_rna.shape[0] > 0: selected_tcga_rna.at[duplicated_gene, "fpkm_unstranded"] = duplicated_rna["fpkm_unstranded"].mean()
                                    else:                           selected_tcga_rna.at[duplicated_gene, "fpkm_unstranded"] = 0

                                    if project == "TCGA-CHOL" and showcount < 50:
                                        print(duplicated_genes)
                                        print(project, rna_fname, duplicated_gene, tcga_rna[tcga_rna["gene_name"] == duplicated_gene]["fpkm_unstranded"].shape, duplicated_rna["fpkm_unstranded"].shape, duplicated_rna["fpkm_unstranded"].mean())
                                        showcount += 1

                                if print_files == True:
                                    with open(target_dir+"\\"+project+"\\"+target_folder+"\\"+rna_fname, "w") as f:
                                        f.write("#" + "\n")
                                        f.write("gene_name" + "\t" + "fpkm_unstranded" + "\n")
                                        f.write("#" + "\n")
                                        f.write("#" + "\n")
                                        f.write("#" + "\n")
                                        f.write("#" + "\n")
                                        selected_tcga_rna.to_csv(path_or_buf=f, sep="\t", header=False, index=False)


# in the current form, it has to be post-processed to remove undesired signs
if mode == "map_types":
    with open(data_dir+"\\"+fname, "r") as file:
        status = json.load(file)

    params   = {"data_dir": data_dir, "os_sep": "\\"}
    cluster  = init_cluster(status, params, cluster_key="Firehouse")

    type_map = {"project": [], "cluster": [], "case_id": [], "cnv": [], "rna": [], "wxs": []}

    for project_key in cluster:
        print(project_key)
        for cluster_key in cluster[project_key]:
            for i in range(len(cluster[project_key][cluster_key]["input"])):
                for case_id in cluster[project_key][cluster_key]["input"][i]:
                    if len(cluster[project_key][cluster_key]["input"][i][case_id]["RNA"]) > 0 and len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"]) > 0:
                        type_map["project"].append(project_key)
                        type_map["cluster"].append(cluster_key)
                        type_map["case_id"].append(case_id)
                        type_map["cnv"].append(json.dumps(cluster[project_key][cluster_key]["input"][i][case_id]["CNV_ranges"]).replace(",", " "))
                        type_map["rna"].append(json.dumps(cluster[project_key][cluster_key]["input"][i][case_id]["RNA"]).replace(",", " "))
                        type_map["wxs"].append(json.dumps(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"]).replace(",", " "))
    
    type_map_df  = pd.DataFrame(type_map)
    type_map_df.to_csv(data_dir+"\\"+fname+"_type_map.txt", sep=",", index=False)


# deviation found on 250321: status.json not existing anymore
# files: ["NIHMS782739-supplement-Tables_S10.txt", "NIHMS782739-supplement-Tables_S3.txt", "status.json", "SCLC_TCGA_shared_genes.txt"]
if mode == "prepare_sclc_integration":
    # mutant mapping SCLC to TCGA nomenclature
    mutant_map = {"frame_shift_del": "Frame_Shift_Del", "frame_shift_ins": "Frame_Shift_Ins", "in_frame_del": "In_Frame_Del", "in_frame_ins": "In_Frame_Ins",
    #              "intron_exon": "Intron", "missense": "Missense_Mutation", "nonsense": "Nonsense_Mutation", "nonstop": "Nonstop_Mutation", "splice": "Splice_Region", "silent": "Silent"}
                  "intron_exon": "Intron", "missense": "Missense_Mutation", "nonsense": "Nonsense_Mutation", "nonstop": "Nonstop_Mutation", "splice": "Splice_Site", "silent": "Silent"}
    
    print_file = False 

    # append status file to contain sclc fnames
    with open(data_dir+"\\"+fnames[2], "r") as file:
        status = json.load(file)

    # load shared genes file for filtering
    with open(data_dir+"\\"+fnames[3], "r") as file:
        shared_genes = pd.read_csv(data_dir+"\\"+fnames[3], delimiter=",")

    # create rna files
    with open(data_dir+"\\"+fnames[0], "r") as file:
        rna = pd.read_csv(data_dir+"\\"+fnames[0], delimiter=";")

    rna = rna[rna["gene"].isin(shared_genes["gene"])].sort_values(by="gene")

    if os.path.isdir(target_dir+"\\shared_RNA") == False:
        os.mkdir(target_dir+"\\shared_RNA")
    
    bar      = IncrementalBar("preparing expression data", max=rna.shape[1]-2)
    cols     = rna.columns.tolist()

    case_ids = {"case_id": [], "submitter_id": [], "tumor_code": [], "tumor_code_id": []}
    file_ids = {}
    showcount = 0
    for i in range(2, rna.shape[1]):
        selected_rna      = pd.DataFrame({"gene_name": rna["gene"], "fpkm_unstranded": rna[cols[i]]})
        duplicated_genes  = selected_rna[selected_rna.duplicated(subset="gene_name")].drop_duplicates(subset="gene_name")["gene_name"].tolist()
        reduced_rna       = selected_rna.drop_duplicates(subset="gene_name")
        reduced_rna.index = reduced_rna["gene_name"]

        case_ids["case_id"].append(cols[i])
        case_ids["submitter_id"].append(cols[i])
        case_ids["tumor_code"].append("SCLC")
        case_ids["tumor_code_id"].append(cols[i])
        file_ids[cols[i]] = {"CNV_genes": [], "CNV_ranges": [], "RNA": [cols[i]], "WXS": [], "CNV_genes_sample_type": [], "CNV_ranges_sample_type": [], "RNA_sample_type": [], "WXS_sample_type": []}


        if print_file == True:
            with open(target_dir+"\\shared_RNA\\"+cols[i], "w") as f:
                f.write("#" + "\n")
                f.write("gene_name" + "\t" + "fpkm_unstranded" + "\n")
                f.write("#" + "\n")
                f.write("#" + "\n")
                f.write("#" + "\n")
                f.write("#" + "\n")

                for duplicated_gene in duplicated_genes:
                    duplicated_rna = selected_rna[selected_rna["gene_name"] == duplicated_gene]
                    duplicated_rna = duplicated_rna[duplicated_rna["fpkm_unstranded"] > 0]
                    if duplicated_rna.shape[0] > 0: reduced_rna.loc[duplicated_gene].loc["fpkm_unstranded"] = duplicated_rna["fpkm_unstranded"].mean()
                    else:                           reduced_rna.loc[duplicated_gene].loc["fpkm_unstranded"] = 0

                    if i >= rna.shape[1]-1 and showcount < 30:
                        showcount += 1
                        print(duplicated_gene, selected_rna[selected_rna["gene_name"] == duplicated_gene].shape, duplicated_rna.shape, duplicated_rna["fpkm_unstranded"].mean())

                reduced_rna.to_csv(path_or_buf=f, sep="\t", header=False, index=False)

        bar.next()
    bar.finish()

    status["case_ids"]["SCLC"] = case_ids
    status["file_ids"]["SCLC"] = file_ids

    # create wxs files
    with open(data_dir+"\\"+fnames[1], "r") as file:
        wxs = pd.read_csv(data_dir+"\\"+fnames[1], delimiter=";")

    # adapt name processes
    wxs = wxs.rename(columns={"Gene_Hugo": "Hugo_Symbol", "Type_1": "Variant_Classification"})
    
    # replace mutation categories
    wxs["Variant_Classification"] = [mutant_map[wxs.iloc[i].loc["Variant_Classification"]] for i in range(wxs.shape[0])]

    if os.path.isdir(target_dir+"\\WXS") == False:
        os.mkdir(target_dir+"\\WXS")

    patient_ids = wxs.drop_duplicates(subset="PAT_ID")["PAT_ID"].tolist()

    bar = IncrementalBar("preparing wxs data", max=len(patient_ids))
    for patient_id in patient_ids:
        selected_wxs = wxs[wxs["PAT_ID"] == patient_id]
        if patient_id in status["file_ids"]["SCLC"]: status["file_ids"]["SCLC"][patient_id]["WXS"].append(patient_id)

        with open(target_dir+"\\WXS\\"+patient_id, "w") as f:
            for i in range(7):
                f.write("#" + "\n")
            
            selected_wxs.to_csv(path_or_buf=f, sep="\t", index=False, lineterminator='\n')

        bar.next()
    bar.finish()

    # save appended status to file
    status = json.dumps(status, indent=4)
    with open(data_dir+"\\status_with_sclc.json", "w") as file:
        file.write(status)


# test calculation of expression calculation for analyze_targets
if mode == "test_expressions":
    selected_projects = ["TCGA-ACC"]

    with open(status_path, "r") as file:
        status = json.load(file)

    targets = pd.read_csv(r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity\Wang_et_al_targets_non_targets.txt", sep=",")
    print("sh", targets.shape[0], targets[[False if len(targets.iloc[i].loc["gene_name"].split(".")) > 1 else True for i in range(targets.shape[0])]].shape[0])
    non_targets = targets[targets["NMD_or_not"] == "N"]
    targets     = targets[targets["NMD_or_not"] == "Y"]


    for project in status["file_ids"]:
        if project in selected_projects:
            for i, case in enumerate(status["file_ids"][project]):
                for j, rna in enumerate(status["file_ids"][project][case]["RNA"]):
                    expressions = load_rna("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\RNA_ccorr_ptc\\"+rna, {"transform_type": None})

                    if expressions.shape[0] > 0:
                        non_targets_expressions = expressions[expressions["gene_name"].isin(non_targets["gene_name"])]
                        non_targets_expressions2 = non_targets_expressions[non_targets_expressions["fpkm_unstranded"] != 0]
                        print(i, j, non_targets_expressions.shape[0])
                        non_targets_expressions = non_targets_expressions[non_targets_expressions["fpkm_unstranded"] > 0]
                        print(non_targets_expressions2[~non_targets_expressions2["gene_name"].isin(non_targets_expressions["gene_name"])])
                        print(i, j, case, "non targets", non_targets_expressions.shape[0], "/", non_targets_expressions2.shape[0], "/",
                             expressions[expressions["gene_name"].isin(targets["gene_name"])].shape[0], "/", non_targets_expressions["fpkm_unstranded"].median())
                        print(expressions[expressions["gene_name"].isin(targets["gene_name"])]["gene_name"].tolist())
                        print(expressions[expressions["gene_name"].isin(targets["gene_name"])]["fpkm_unstranded"].tolist())
                        input("x")


# test validity of gene-specific expression data generated by analyze_targets (250429)
if mode == "test_expressions2":
    tcga_dir      = r"C:\Programming\Translational_genomics\TCGA"
    path          = r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity\2025-04-30_18-18-31_TCGA_selected_genes_TPM_ccorr_HLA-C\cancer_scores_TCGA_selected_genes_TPM_ccorr_HLA-C"
    status_path   = r"C:\Programming\Translational_genomics\NMD_analysis\data\filtered_ptc_status.json"
    target_folder = "RNA"
    target_symbol = "HLA-C"
    test_projects = ["TCGA-BRCA"]
    
    cancer_scores = pd.read_csv(path, sep=",")
    cancer_scores = cancer_scores[cancer_scores["project"].isin(test_projects)]
    cancer_scores = convert_raw_data(cancer_scores)

    with open(status_path, "r") as file:
        status = json.load(file)

    tests = 0; mismatches = 0
    for i in range(cancer_scores.shape[0]):
        print(i, cancer_scores.iloc[i].loc["project"])
        sample_names = json.loads(cancer_scores.iloc[i].loc["sample_names"].replace("'", "\""))

        for j, case_id in enumerate(sample_names):
            rna_fnames = status["file_ids"][cancer_scores.iloc[i].loc["project"]][case_id]["RNA"]

            # for simplicity, only single rna files are inferred
            if len(rna_fnames) == 1:
                rna = load_rna(tcga_dir+"\\"+cancer_scores.iloc[i].loc["project"]+"\\"+target_folder+"\\"+rna_fnames[0], {"transform_type": None})
                tests += 1

                #print(tcga_dir+"\\"+cancer_scores.iloc[i].loc["project"]+"\\"+target_folder+"\\"+rna_fnames[0],
                #      rna[rna["gene_name"] == target_symbol].iloc[0].loc["fpkm_unstranded"], "/", cancer_scores.iloc[i].loc["fpkm_unstranded_raw"][j])

                if rna[rna["gene_name"] == target_symbol].iloc[0].loc["tpm_unstranded"] != cancer_scores.iloc[i].loc["tpm_unstranded_raw"][j]:
                    print("< mismatch @", cancer_scores.iloc[i].loc["project"], "/", rna_fnames[0], ":",
                          rna[rna["gene_name"] == target_symbol].iloc[0].loc["tpm_unstranded"], "/", cancer_scores.iloc[i].loc["tpm_unstranded_raw"][j])
                    mismatches += 1

    print("<", mismatches, "mismatches in", tests, "tests")


# test validity of downloaded files
if mode == "test_md5sum":
    with open(status_path, "r") as file:
        status = json.load(file)

    for project in status["file_ids"]:
        for i, case in enumerate(status["file_ids"][project]):
            for j, el in enumerate(status["file_ids"][project][case]):
                if "md5sum" not in el and "sample_type" not in el:
                    for k, fname in enumerate(status["file_ids"][project][case][el]):
                        path = "C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"\\"+fname
                        if el != "WXS" and os.path.isfile(path) == True:
                            md5 = calculate_md5(path)
                            #print(i, j, k, project, case, el, status["file_ids"][project][case][el][k], md5, status["file_ids"][project][case][el+"_md5sum"][k])

                            if md5 != status["file_ids"][project][case][el+"_md5sum"][k]:
                                print("deviation @", i, j, k, project, case, el, status["file_ids"][project][case][el][k])
                                print(" ", md5, "/", status["file_ids"][project][case][el+"_md5sum"][k])

                        if el == "WXS" and os.path.isfile(path) == True:
                            test_path = "C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed"+"\\"+fname
                            if os.path.isdir("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed") == False:
                                os.mkdir("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed")
                                
                            compress(path, test_path)
                            md5 = calculate_md5(test_path)
                            #print(i, j, k, project, case, el, status["file_ids"][project][case][el][k], md5, status["file_ids"][project][case][el+"_md5sum"][k])

                            if md5 != status["file_ids"][project][case][el+"_md5sum"][k]:
                                print("deviation @", i, j, k, project, case, el, status["file_ids"][project][case][el][k])
                                print(" ", md5, "/", status["file_ids"][project][case][el+"_md5sum"][k])


# tests the discrepancy of ptc mutations called by analyze_targets and prepare_data, respectively
if mode == "test_ptc_mutations":
    # required data
    cancer_score_path      = r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity\2025-06-23_16-06-03_TCGA_NMD_targets_analysis_FPKM_exp_ccorr_test\cancer_scores_TCGA_NMD_targets_analysis_FPKM_exp_ccorr_test" 
    genome_path            = r"C:\Programming\Translational_genomics\NMD_analysis\data\hg38_knownGene_appended.txt"
    processed_variant_path = r"C:\Programming\Translational_genomics\NMD_analysis\data\tcga.txt"
    status_path            = r"C:\Programming\Translational_genomics\NMD_analysis\data\filtered_status_ptc.json"
    variant_path           = r"C:\Programming\Translational_genomics\NMD_analysis\data\tcga_raw.txt"

    check_data_preparation = True
    feature1               = "FEATURE:ptc_mutations2"
    feature2               = "ptc_mutations"
    tcga_dir               = r"C:\Programming\Translational_genomics\TCGA"
    wxs_identifier         = "Transcript_ID" # "Gene" "Transcript_ID"
    filter_mode            = False # filters specified below are applied, smaller set of hypotheses is tested, full test should be re-run without filters

    # data loading
    cancer_scores = pd.read_csv(cancer_score_path, delimiter=",")
    cancer_scores = cancer_scores.sort_values("project", ascending=True)

    # load all lists stored as str to lists
    cancer_scores = convert_raw_data(cancer_scores)
    hg38          = pd.read_csv(genome_path, delimiter=",")
    tcga          = pd.read_csv(processed_variant_path, delimiter=",")
    tcga_raw      = pd.read_csv(variant_path, delimiter="\t")

    if filter_mode == False:
        # extended stats used to check unfiltered data
        discrepancy_stats = {"no_discrepancy": 0, "filtered_ptcs_exceed_input": {"accounted": 0, "unaccounted": 0}, "called_ptcs_deviate1": {"accounted": 0, "unaccounted": 0},
                            "called_ptcs_deviate2": {"accounted": 0, "unaccounted": 0}, "all_ptcs_deviate": {"accounted": 0, "unaccounted": 0}}

    # apply value filter to tcga and tcga_raw to get same subset with cancer_scores
    if filter_mode == True:
        # minimum stats used to check filtered data
        discrepancy_stats = {"no_discrepancy": 0, "filtered_ptcs_exceed_input": {"accounted": 0, "unaccounted": 0}}

        tcga     = tcga[tcga["ID:ptc reads"] >= 1]
        tcga     = tcga[tcga["ID:cnv total"] != 0]
        #tcga_raw = tcga_raw[tcga_raw["project"] == "TCGA-ACC"]
        tcga_raw = tcga_raw[tcga_raw["RNASEQ_ptc_fpkm_unstranded"] >= 1]

        # load clinical data and append gender info to raw file
        clinical_data = pd.read_csv(clinical_data_path, delimiter=",")
        tcga_raw      = map_data([tcga_raw, clinical_data], ["gender"], ["case_id", "case id"], "appending clinical data", tag="ID")

        # filter out variants for male patients and Y- and X-chromosomes (analyze_predictions)
        tcga_raw = tcga_raw[tcga_raw["Chromosome"] != "chrY"]
        tcga_raw = tcga_raw[[True if tcga_raw.iloc[i].loc["Chromosome"] != "chrX" or (tcga_raw.iloc[i].loc["Chromosome"] == "chrX"
                             and tcga_raw.iloc[i].loc["ID:gender"] == "FEMALE") else False for i in range(tcga_raw.shape[0])]]

        # remove CNV=0 (analyze_predictions)
        tcga_raw = tcga_raw[[False if pd.isna(tcga_raw.iloc[i].loc["RNASEQ_ptc_cnv_total"]) == False
                             and float(tcga_raw.iloc[i].loc["RNASEQ_ptc_cnv_total"]) == 0 else True
                             for i in range(tcga_raw.shape[0])]]

        print("sh", tcga[tcga["ID:cnv total"].isna()].shape[0], tcga_raw[tcga_raw["RNASEQ_ptc_cnv_total"].isna()].shape[0])

    hg38["transcript id"]     = [hg38.iloc[i].loc["transcript id"].split(".")[0] for i in range(hg38.shape[0])]
    hg38.index                = hg38["transcript id"]
    tcga_raw["Transcript_ID"] = [tcga_raw.iloc[i].loc["Transcript_ID"].split(".")[0] for i in range(tcga_raw.shape[0])]

    with open(status_path, "r") as file:
        status = json.load(file)
  
    #cancer_scores = cancer_scores[cancer_scores["project"] == "TCGA-PRAD"]
    selected_variants = {"variant_id": []}

    for i in range(cancer_scores.shape[0]):
        print("<", cancer_scores.iloc[i].loc["project"])
        values1, values2, _, sample_names = get_values(cancer_scores, {}, feature1, feature2, i, "FEATURE:ptc_mutations2", extract_ptc=True, get_sample_names=True)

        for j in range(len(values1)):
            # determine count in raw file
            selected_tcga_raw = tcga_raw[tcga_raw["case_id"] == sample_names[j]]
            #print("j", j, selected_tcga_raw["ID:variant id"].tolist(), selected_tcga_raw["ID:variant id"])

            # load wxs file(s) to conduct further tests
            wxs_fnames = status["file_ids"][cancer_scores.iloc[i].loc["project"]][sample_names[j]]["WXS"]
            
            called_ptcs_deviate1_testsizes = []; all_ptcs_deviate_testsizes = []
            wxss = []
            for wxs_fname in wxs_fnames:
                wxs     = load_wxs(tcga_dir+"//"+cancer_scores.iloc[i].loc["project"]+"//WXS//"+wxs_fname)
                ptc_wxs = create_ptc_subset(wxs, cancer_scores.iloc[i].loc["project"], wxs_fname)
                
                # remove variants with upstream duplications/insertions/deletions
                keep_variants = [True for _ in range(ptc_wxs.shape[0])]
                ptcs          = []
                for k in range(ptc_wxs.shape[0]):
                    ptc_test      = get_numbers(ptc_wxs.iloc[k].loc["HGVSp"])
                    ptc           = get_ptc(ptc_wxs.iloc[k].loc["HGVSp"])
                    ptcs.append(ptc)

                    if ptc == None or ptc_wxs.iloc[k].loc["Variant_Classification"] not in ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"]:
                        keep_variants[k] = False
                    
                    else:
                        if ptc != None and "_" in ptc_wxs.iloc[k].loc["HGVSp"] and len(ptc_test) > 1 and np.sum(ptc_test) < ptc+10:
                            print("ptc deviation @", wxs_fnames[0], ptc_wxs.iloc[k].loc["HGVSp"], ptc, "/", np.sum(ptc_test))

                        selected_wxs  = wxs[wxs[wxs_identifier] == ptc_wxs.iloc[k].loc[wxs_identifier]]
                        selected_wxs  = selected_wxs[selected_wxs["HGVSp"] != ptc_wxs.iloc[k].loc["HGVSp"]]

                        # check for additional ptcs
                        additional_ptcs = selected_wxs[[True if type(selected_wxs.iloc[l].loc["HGVSp"]) == str and "Ter" in selected_wxs.iloc[l].loc["HGVSp"] else False
                                                        for l in range(selected_wxs.shape[0])]]
                        #if additional_ptcs.shape[0] > 0:
                        #    print(additional_ptcs["HGVSp"].tolist())
                        #    print(len(test), [get_ptc(additional_ptcs.iloc[l].loc["HGVSp"]) for l in range(additional_ptcs.shape[0])])

                        selected_wxs  = selected_wxs[[True if type(selected_wxs.iloc[l].loc["HGVSp"]) == str and 
                                                     ("del" in selected_wxs.iloc[l].loc["HGVSp"] or "dup" in selected_wxs.iloc[l].loc["HGVSp"]
                                                      or "ins" in selected_wxs.iloc[l].loc["HGVSp"]) else False
                                                      for l in range(selected_wxs.shape[0])]]
                        
                        positions     = [get_numbers(selected_wxs.iloc[l].loc["HGVSp"])[0] for l in range(selected_wxs.shape[0]) if len(get_numbers(selected_wxs.iloc[l].loc["HGVSp"])) > 0]
                        test          = [position for position in positions if position < ptc]
                        test.extend([l for l in range(additional_ptcs.shape[0]) if get_ptc(additional_ptcs.iloc[l].loc["HGVSp"]) != None and get_ptc(additional_ptcs.iloc[l].loc["HGVSp"]) < ptc])
                        if len(test) > 0: keep_variants[k] = False

                # create copy for testing and conduct modifications without filtering of keep variants
                ptc_wxs_test = ptc_wxs
                ptc_wxs_test = ptc_wxs_test[[False if "ext" in ptc_wxs_test.iloc[k].loc["HGVSp"] or "?" in ptc_wxs_test.iloc[k].loc["HGVSp"] or "=" in ptc_wxs_test.iloc[k].loc["HGVSp"] else True
                                            for k in range(ptc_wxs_test.shape[0])]] # remove entries with "ext" tag (no-stop variants) or "?" tag or "=" tag
                ptc_wxs_test = ptc_wxs_test[[False if len([1 for l in range(len(ptc_wxs_test.iloc[k].loc["HGVSp"])-3) if ptc_wxs_test.iloc[k].loc["HGVSp"][l:l+3] == "Ter"]) == 2 else True
                                            for k in range(ptc_wxs_test.shape[0])]] # remove entries with two "Ter" tags (no-stop variants)
                ptc_wxs_test = ptc_wxs_test.drop_duplicates(subset=[wxs_identifier]) # remove duplicate Genes
                ptc_wxs_test = ptc_wxs_test.drop_duplicates(subset=[wxs_identifier, "HGVSp"]) # remove redundant entries due to multiple WXS files
                all_ptcs_deviate_testsizes.append(ptc_wxs_test.shape[0])

                # continue main test
                ptc_wxs.insert(0, "ptc_pos", ptcs)
                ptc_wxs = ptc_wxs[keep_variants]
                ptc_wxs = ptc_wxs[[False if "ext" in ptc_wxs.iloc[k].loc["HGVSp"] or "?" in ptc_wxs.iloc[k].loc["HGVSp"] or "=" in ptc_wxs.iloc[k].loc["HGVSp"] else True
                                   for k in range(ptc_wxs.shape[0])]] # remove entries with "ext" tag (no-stop variants) or "?" tag or "=" tag
                ptc_wxs = ptc_wxs[[False if len([1 for l in range(len(ptc_wxs.iloc[k].loc["HGVSp"])-3) if ptc_wxs.iloc[k].loc["HGVSp"][l:l+3] == "Ter"]) == 2 else True
                                   for k in range(ptc_wxs.shape[0])]] # remove entries with two "Ter" tags (no-stop variants)
                
                called_ptcs_deviate1_testsizes.append(ptc_wxs.shape[0])
                ptc_wxs            = ptc_wxs.drop_duplicates(subset=[wxs_identifier]) # remove duplicate Genes
                ptc_wxs            = ptc_wxs.drop_duplicates(subset=[wxs_identifier, "ptc_pos"]) # remove redundant entries due to multiple WXS files
                wxss.append(ptc_wxs)

            wxs_avg = np.mean([wxss[k].shape[0] for k in range(len(wxss))]) # calculate mean to compare with ptc mutations call from analyze_mutations
            
            if len(wxss) == 1: wxs = wxss[0]
            if len(wxss) > 1:  wxs = pd.concat(wxss)
            wxs = wxs.drop_duplicates(subset=[wxs_identifier, "ptc_pos"]) # remove redundant entries due to multiple WXS files

            # check for errors
            # with filters, only one error type is checked

            # ptc count from prepare_data can be smaller due to filters but never larger
            if "filtered_ptcs_exceed_input" in discrepancy_stats and values1[j] > selected_tcga_raw.shape[0]:
                discrepancy_stats["filtered_ptcs_exceed_input"]["unaccounted"] += 1; error = "filtered_ptcs_exceed_input"

            elif filter_mode == True: discrepancy_stats["no_discrepancy"] += 1; error = "no_discrepancy"

            else:
                if "no_discrepancy" in discrepancy_stats and wxs_avg == values2[j] and wxs.shape[0] == selected_tcga_raw.shape[0]:
                    discrepancy_stats["no_discrepancy"] += 1; error = "no_discrepancy"
                
                # ptc count from analyze_mutations must be identical to test value
                if "called_ptcs_deviate1" in discrepancy_stats and wxs_avg != values2[j] and wxs.shape[0] == selected_tcga_raw.shape[0]:
                    discrepancy_stats["called_ptcs_deviate1"]["unaccounted"] += 1; error = "called_ptcs_deviate1"

                # ptc count from extract_mutations must be identical to test value
                if "called_ptcs_deviate2" in discrepancy_stats and wxs_avg == values2[j] and wxs.shape[0] != selected_tcga_raw.shape[0]:
                    discrepancy_stats["called_ptcs_deviate2"]["unaccounted"] += 1; error = "called_ptcs_deviate2"

                # all values deviate
                if "all_ptcs_deviate" in discrepancy_stats and wxs_avg != values2[j] and wxs.shape[0] != selected_tcga_raw.shape[0] and wxs_avg != selected_tcga_raw.shape[0]:
                    discrepancy_stats["all_ptcs_deviate"]["unaccounted"] += 1; error = "all_ptcs_deviate"


            if error != "no_discrepancy":
                print("<", error)
                print(cancer_scores.iloc[i].loc["project"], "index", j,
                      "cancer scores/sample names", sample_names[j], "cancer scores/", feature1, values1[j], "cancer scores/", feature2, values2[j],
                      "wxs/avg", wxs_avg, "raw/count", selected_tcga_raw.shape[0], "wxs/shape", wxs.shape[0])

                if wxs.shape[0] > 0:
                    wxs = wxs.sort_values(by="Gene")
                    print(wxs_fnames)

                selected_tcga_raw = selected_tcga_raw.sort_values(by="Gene")

                if wxs.shape[0] > selected_tcga_raw.shape[0]:
                    print("deviation1")
                    deviating_wxs = wxs[~wxs["HGVSp"].isin(selected_tcga_raw["HGVSp"])]
                    for k in range(deviating_wxs.shape[0]):
                        print(" ", deviating_wxs.iloc[k].loc["Gene"], deviating_wxs.iloc[k].loc["HGVSp"])

                if wxs.shape[0] < selected_tcga_raw.shape[0]:
                    print("deviation2")
                    deviating_tcga_raw = selected_tcga_raw[~selected_tcga_raw["HGVSp"].isin(wxs["HGVSp"])]
                    for k in range(deviating_tcga_raw.shape[0]):
                        print(" ", deviating_tcga_raw.iloc[k].loc["Gene"], deviating_tcga_raw.iloc[k].loc["HGVSp"])

                if values1[j] > selected_tcga_raw.shape[0]:
                    print("deviation3")
                    selected_tcga = tcga[tcga["ID:case id"] == sample_names[j]]
                    for k in range(selected_tcga.shape[0]):
                        print(" ", selected_tcga.iloc[k].loc["ID:variant id"], selected_tcga.iloc[k].loc["ID:HGVSp"])


            if check_data_preparation == True:
                # filter steps for raw tcga data
                ptcs          = [3*(get_ptc(selected_tcga_raw.iloc[k].loc["HGVSp"])-1) for k in range(selected_tcga_raw.shape[0])]

                # filter out no-stop variants
                cds_positions = [int(selected_tcga_raw.iloc[k].loc["CDS_position"].split("/")[1].replace("\'", "").replace("]", ""))
                                 if len(selected_tcga_raw.iloc[k].loc["CDS_position"]) > 0 else None
                                 for k in range(selected_tcga_raw.shape[0])]
                filtered_ptcs = selected_tcga_raw[[True if cds_positions[k] != None and ptc < cds_positions[k]-3 else False for k, ptc in enumerate(ptcs)]]

                # filter out variants not contained in the current database (retired transcripts)
                if filtered_ptcs.shape[0] > 0: filtered_ptcs = filtered_ptcs[filtered_ptcs["Transcript_ID"].isin(hg38["transcript id"])]

                # filter out variants that have inconsistent cds sizes
                cds_sizes     = [get_cds_size(hg38, filtered_ptcs.iloc[k].loc["Transcript_ID"]) for k in range(filtered_ptcs.shape[0])]                
                filtered_ptcs = filtered_ptcs[[True if cds_sizes[k]-3 == 3*int(filtered_ptcs.iloc[k].loc["Protein_position"].split("/")[1]) else False
                                               for k in range(filtered_ptcs.shape[0])]]

                if filter_mode == True and filtered_ptcs.shape[0] > 0: filtered_ptcs = filtered_ptcs[filtered_ptcs["RNASEQ_ptc_cnv_total"] != 0]

                if filtered_ptcs.shape[0] > values1[j]:
                    print(cancer_scores.iloc[i].loc["project"], "index", j, sample_names[j])
                    print("deviation4:", filtered_ptcs.shape[0], "/", values1[j])
                    selected_tcga           = tcga[tcga["ID:case id"] == sample_names[j]]
                    deviating_filtered_ptcs = filtered_ptcs[~filtered_ptcs["HGVSp"].isin(selected_tcga["ID:HGVSp"])]
                    for k in range(deviating_filtered_ptcs.shape[0]):
                        print(" ", hg38.loc[deviating_filtered_ptcs.iloc[k].loc["Transcript_ID"]].loc["strand"],
                              deviating_filtered_ptcs.iloc[k].loc["variant_id"], deviating_filtered_ptcs.iloc[k].loc["CDS_position"],
                              deviating_filtered_ptcs.iloc[k].loc["HGVSc"], deviating_filtered_ptcs.iloc[k].loc["HGVSp"])

                        selected_variants["variant_id"].append(deviating_filtered_ptcs.iloc[k].loc["variant_id"])

                # small deviations can occur if mutations cause PTC in normal PTC site for which no prediction is calculated
                if filtered_ptcs.shape[0] < values1[j]:
                    print(cancer_scores.iloc[i].loc["project"], "index", j, sample_names[j])
                    print("deviation5:", filtered_ptcs.shape[0], "/", values1[j])
                    selected_tcga = tcga[tcga["ID:case id"] == sample_names[j]]

                    if filtered_ptcs.shape[0] > 0:
                        deviating_selected_tcga = selected_tcga[~selected_tcga["ID:HGVSp"].isin(filtered_ptcs["HGVSp"])]
                        for k in range(deviating_selected_tcga.shape[0]):
                            print(" ", deviating_selected_tcga.iloc[k].loc["ID:variant id"], deviating_selected_tcga.iloc[k].loc["ID:HGVSp"])

                    else:
                        for k in range(selected_tcga.shape[0]):
                            print(" ", selected_tcga.iloc[k].loc["ID:variant id"], selected_tcga.iloc[k].loc["ID:HGVSp"])

    selected_variants = pd.DataFrame(selected_variants)
    selected_variants.to_csv(data_dir+"\\selected_variants.txt", index=False, sep=",")
    print(json.dumps(discrepancy_stats, indent=4))

    


    
