import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
from scipy.cluster.hierarchy import linkage, dendrogram
import shutil
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from extract_mutations_utils import *
from tcga_tools_utils import *


clinical_data_path     = parent_dir+r"\data\tcga_survival_data.txt"
data_dir               = parent_dir+r"\data"
fname                  = "tcga_data_info.json"
fnames                 = ["tcga_variants.txt", "tcga_survival_data.txt"]
immune_data_path       = parent_dir+r"\data\tcga_immune_scores.txt"
# "append_clinical_data"
# "check_files" "check_sample_types" "check_sample_ids" "compare_status_files"
# "filter_sample_types"
mode                   = "append_clinical_data"
status_path            = parent_dir+r"\data\tcga_data_info.json"
target_dir             = parent_dir+r"\data"


#  used to create tcga_survival_data.txt (appending case ids based on submitter id)
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


# clinical data file: tcga_survival_data.txt
# gender data needed for exclusion of monoallelic variants (male X)
if mode == "append_clinical_data":
    features = ["gender"] # ["age_at_initial_pathologic_diagnosis", "ajcc_pathologic_tumor_stage", "gender", "OS", "OS.time", "DFI", "DFI.time", "DSS", "DSS.time", "PFI", "PFI.time"]
    data     = []
    for fname in fnames:
        data.append(pd.read_csv(data_dir+"\\"+fname, delimiter=","))

    data[0] = map_data(data, features, ["ID:case id", "case id"], "appending clinical data", tag="LABEL")
    data[0].to_csv(data_dir+"\\"+fnames[0].split(".")[0]+"_appended.txt", sep=",")


# appending tcga_immune_scores.txt
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
                        if current_case[target_col+"_sample_type"][i][j] in current_ids: current_ids[current_case[target_col+"_sample_type"][i][j]].append(current_case[target_col+"_submitter_id"][i][j].split("-")[-1])
                        else:                                                            current_ids[current_case[target_col+"_sample_type"][i][j]] = [current_case[target_col+"_submitter_id"][i][j].split("-")[-1]]
                        if current_case[target_col+"_sample_type"][i][j] in total_ids:   total_ids[current_case[target_col+"_sample_type"][i][j]].append(current_case[target_col+"_submitter_id"][i][j].split("-")[-1])
                        else:                                                            total_ids[current_case[target_col+"_sample_type"][i][j]] = [current_case[target_col+"_submitter_id"][i][j].split("-")[-1]]
            
                for key in current_ids:
                    if len(np.unique(current_ids[key])) > 1:
                        print(project, case_id, key, current_ids[key])
                        sample_stats[target_col+"_multiple"] += 1
                    
                    sample_stats[target_col+"_total"] += 1

            bar.next()
    bar.finish()
    print(json.dumps(sample_stats, indent=4))


# comparison of unfiltered and filtered status
if mode == "check_sample_types":
    target_folders = ["RNA_c_ptc"]

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


if mode == "filter_sample_types":
    extraction_size    = 50
    extraction_targets = [] # ["Solid Tissue Normal"] # define categories to be extracted to novel project
    # define categories to be preserved in the updated status
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
                    
                    for j in range(len(status["file_ids"][project][case_id][key+"_sample_type"][i])):
                        if j < len(status["file_ids"][project][case_id][key+"_sample_type"][i])-1:
                            sample_type += status["file_ids"][project][case_id][key+"_sample_type"][i][j] + ", "

                        else:
                            sample_type += status["file_ids"][project][case_id][key+"_sample_type"][i][j]

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
    type_map_df.to_csv(data_dir+"\\TCGA_type_map_filtered_ptc.txt", sep=",", index=False)


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

                            if md5 != status["file_ids"][project][case][el+"_md5sum"][k]:
                                print("deviation @", i, j, k, project, case, el, status["file_ids"][project][case][el][k])
                                print(" ", md5, "/", status["file_ids"][project][case][el+"_md5sum"][k])

                        if el == "WXS" and os.path.isfile(path) == True:
                            test_path = "C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed"+"\\"+fname
                            if os.path.isdir("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed") == False:
                                os.mkdir("C:\\Programming\\Translational_genomics\\TCGA\\"+project+"\\"+el+"_compressed")
                                
                            compress(path, test_path)
                            md5 = calculate_md5(test_path)

                            if md5 != status["file_ids"][project][case][el+"_md5sum"][k]:
                                print("deviation @", i, j, k, project, case, el, status["file_ids"][project][case][el][k])
                                print(" ", md5, "/", status["file_ids"][project][case][el+"_md5sum"][k])
