from ast import literal_eval
import collections
from datetime import datetime
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
import scipy
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from extract_mutations_utils import *


def _count_by_threshold(selected_features, score_threshold):
    return len([selected_feature for selected_feature in selected_features if selected_feature >= score_threshold[0] and selected_feature <= score_threshold[1]])


def _mean(values):
    return np.mean([value for value in values if pd.isna(value) == False])


def append_sample_scores(extended_features, sample_scores, case_ids, params):
    extended_features[params["file_tag"].replace("_", "")] = [None for _ in range(extended_features.shape[0])]

    bar = IncrementalBar(set_bar("appending sample scores"), max=len(sample_scores))
    for project in sample_scores:
        for sample in sample_scores[project][params["info"]["rna"][0]]:
            case_id = case_ids[project][sample]

            selected_index = extended_features[extended_features["ID:case id"] == case_id].index.tolist()
            for i in selected_index:
                extended_features.at[extended_features.index[i], params["appended_col"]] = sample_scores[project][params["info"]["rna"][0]][sample]["median"]
        
        bar.next()
    bar.finish()
    extended_features.to_csv(path_or_buf=params["extended_features_path"].strip(".txt")+params["file_tag"]+".txt", sep=",", index=False)
    print("< features extended.")
    return


def calculate_target_expression(targets, extended_features, params):
    bar = IncrementalBar(set_bar("project blocks are determined"), max=targets.shape[0])

    targets      = targets.sort_values(by=["project", "cluster_key"], ignore_index=True)
    #targets      = targets[targets["project"].isin(["TCGA-UCEC"])] #
    
    last_cluster = None
    last_project = None
    block_index  = []
    temp         = []

    for i in range(targets.shape[0]):
        if last_cluster == None or (last_cluster == targets.iloc[i].loc["cluster_key"] and last_project == targets.iloc[i].loc["project"]):
            temp.append(i)

        elif len(temp) > 0:
            if targets.iloc[temp].drop_duplicates(subset="project").shape[0] != 1:
                print("< error1 @calculate_target_expression. inconsistent no. of project names in block index.")
                print(targets.iloc[temp].drop_duplicates(subset="project"))

            block_index.append({"block id": last_cluster, "index": temp})
            temp = [i]

        last_cluster = targets.iloc[i].loc["cluster_key"]
        last_project = targets.iloc[i].loc["project"]

        bar.next()
    bar.finish()

    # append last block
    if len(temp) > 0: block_index.append({"block id": last_cluster, "index": temp})

    # initialize containers to store calculation results
    projects            = targets.drop_duplicates(subset="project")["project"].tolist()
    target_ids          = targets.drop_duplicates(subset=params["target_identifier"]["rna"])[params["target_identifier"]["rna"]].tolist()
    nmd_target_averages = {project: {rna_value: {target_id: [] for target_id in target_ids} for rna_value in params["info"]["rna"]} for project in projects}

    # calculate NMD target averages for filtering
    bar = IncrementalBar(set_bar("filtering is prepared"), max=len(block_index))
    for i in range(len(block_index)):
        # test for contaminations in block index
        if targets.iloc[block_index[i]["index"]].drop_duplicates(subset="project").shape[0] != 1:
            print("< error2 @calculate_target_expression. inconsistent no. of project names in block index.")
            print(targets.iloc[block_index[i]["index"]].drop_duplicates(subset="project"))
        
        for rna_value in params["info"]["rna"]:
            project = targets.iloc[block_index[i]["index"][0]].loc["project"]
            
            for j in block_index[i]["index"]:               
                # iterate over samples to calculate per-sample averages
                for k in range(len(targets.iloc[j].loc[rna_value])):
                    # apply filter to exclude values below threshold
                    target_avg, _, _ = calculate_value(json.dumps([targets.iloc[j].loc[rna_value][k]]), params["rna_calculation"],
                                                       get_stats=True, threshold=params["target_threshold"])
                    #if project == "TCGA-ACC": print(i, j, k, "target_avg", target_avg, targets.iloc[j].loc[rna_value][k])
                    nmd_target_averages[project][rna_value][targets.iloc[j].loc[params["target_identifier"]["rna"]]].append(target_avg)

        bar.next()
    bar.finish()
    

    # set up sample score container using determined sample sizes per project
    sample_scores       = {project: {**{rna_value         : {sample: [] for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][rna_value].keys())[0]]))}
                                                             for rna_value in params["info"]["rna"] if len(list(nmd_target_averages[project][rna_value].keys())) > 0},
                                     **{col               : {sample: None for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][params["info"]["rna"][0]].keys())[0]]))
                                                             if len(list(nmd_target_averages[project][params["info"]["rna"][0]].keys())) > 0} for col in params["info"]["appended_features"]},
                                     **{col               : {sample: None for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][params["info"]["rna"][0]].keys())[0]]))
                                                             if len(list(nmd_target_averages[project][params["info"]["rna"][0]].keys())) > 0} for col in params["info"]["extended_features"]},
                                     **{col               : {sample: None for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][params["info"]["rna"][0]].keys())[0]]))
                                                             if len(list(nmd_target_averages[project][params["info"]["rna"][0]].keys())) > 0} for col in params["info"]["mutations"]},
                                     **{"sample_names"    : {sample: 0 for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][params["info"]["rna"][0]].keys())[0]]))
                                                             if len(list(nmd_target_averages[project][params["info"]["rna"][0]].keys())) > 0}},
                                     **{"target_mutations": {sample: 0 for sample in range(len(nmd_target_averages[project][rna_value][list(nmd_target_averages[project][params["info"]["rna"][0]].keys())[0]]))
                                                             if len(list(nmd_target_averages[project][params["info"]["rna"][0]].keys())) > 0}}
                                    } for project in projects}
    
    nmd_target_averages = {project: {rna_value: {target_id: _mean(nmd_target_averages[project][rna_value][target_id])
                                                 for target_id in target_ids} for rna_value in params["info"]["rna"]} for project in projects}

    # apply filter
    filtered_index      = {project: {rna_value: [target_id for target_id in nmd_target_averages[project][rna_value]
                                                 if nmd_target_averages[project][rna_value][target_id] > params["target_filter"]]
                                                 for rna_value in params["info"]["rna"]} for project in projects}
    
    # initialize container to store case ids per project
    case_ids = {project: [] for project in projects}
    

    # calculate NMD sample scores
    bar = IncrementalBar(set_bar("sample scores are calculated"), max=len(block_index))
    last_project = None
    filtered_out = 0

    # iterate over cluster blocks
    for i in range(len(block_index)):
        project           = targets.iloc[block_index[i]["index"][0]].loc["project"]
        case_ids[project] = targets.iloc[block_index[i]["index"][0]].loc["case_ids"]

        if last_project == None or last_project != project:
            if last_project != None and sample_index != len(sample_scores[last_project][params["info"]["rna"][0]]):
                print("< sample index count inconsistent with input sample size @", last_project, sample_index, "/", len(sample_scores[last_project][params["info"]["rna"][0]]))

            sample_index = 0

        # iterate over samples per cluster
        for j in range(len(targets.iloc[block_index[i]["index"][0]].loc[params["info"]["rna"][0]])):
            #print("j", j, len(targets.iloc[block_index[i]["index"][0]].loc[params["info"]["rna"][0]]))
            # iterate over index per cluster block
            for k in block_index[i]["index"]:
                for rna_value in params["info"]["rna"]:
                    # register filenames
                    sample_scores[project]["sample_names"][sample_index] = case_ids[project][sample_index]

                    # filtering step
                    if targets.iloc[k].loc[params["target_identifier"]["rna"]] in filtered_index[project][rna_value]:
                        if params["target_normalization"] == False:
                            target_avg, _, _ = calculate_value(json.dumps([targets.iloc[k].loc[rna_value][j]]), params["rna_calculation"],
                                                               get_stats=True, threshold=params["target_threshold"])
                            #if j == 0:
                            #    print(project, "i1", i, "j", j, "sample_index", sample_index, "k", k, targets.iloc[k].loc[params["target_identifier"]["rna"]], targets.iloc[k].loc[rna_value][j],
                            #          len(targets.iloc[k].loc[rna_value]), target_avg)
                                
                        if params["target_normalization"] == True:
                            normalization_factor = targets.iloc[k].loc[rna_value+"_non_targets"][j]
                            target_avg, _, _     = calculate_value(json.dumps([targets.iloc[k].loc[rna_value][j]]), params["rna_calculation"],
                                                                   get_stats=True, threshold=params["target_threshold"], normalization_factor=normalization_factor)
                    
                        if target_avg != None and pd.isna(target_avg) == False:
                            if sample_index < len(sample_scores[project][rna_value]):
                                sample_scores[project][rna_value][sample_index].append(target_avg)

                            else:
                                print("< dimension error1 @calculate_target_expression: ", sample_index, "/", len(sample_scores[project][rna_value]))

                                
                        if (params["target_normalization"] == True and 
                            len(targets.iloc[k].loc[rna_value]) != len(targets.iloc[k].loc[rna_value+"_non_targets"])):
                            print("< dimension error2 @calculate_target_expression: ", len(targets.iloc[k].loc[rna_value]),
                                  "/", len(targets.iloc[k].loc[rna_value+"_non_targets"]))
                            
                    else:
                        filtered_out += 1
            
            for col in params["info"]["appended_features"]:
                sample_scores[project][col][sample_index] = targets.iloc[k].loc[col][sample_index]
                #print("k", k, project, col, sample_index, targets.iloc[k].loc[col][sample_index])

            for col in params["info"]["combined_mutations"]:
                target_sum = 0
                for mutation_col in params["info"]["combined_mutations"][col]:
                    target_avg, _, _ = calculate_value(json.dumps([targets.iloc[k].loc[mutation_col][j]]), params["rna_calculation"], get_stats=True)
                    #if len(targets.iloc[k].loc[mutation_col][j]) > 0: print(mutation_col, target_avg, targets.iloc[k].loc[mutation_col][j])
                    if target_avg != None: target_sum += target_avg

                sample_scores[project][col][sample_index] = target_sum

            for col in params["info"]["mutations"]:
                if col not in params["info"]["combined_mutations"]:
                    target_avg, _, _                          = calculate_value(json.dumps([targets.iloc[k].loc[col][j]]), params["rna_calculation"], get_stats=True)
                    sample_scores[project][col][sample_index] = target_avg

            if len(params["info"]["target_mutations"]) > 0 and "variant_classifications" in targets.columns and type(targets.iloc[k].loc["variant_classifications"]) == list:
                variant_found = False
                for variant_list in targets.iloc[k].loc["variant_classifications"][j]:
                    for variant in variant_list:
                        if variant in params["info"]["target_mutations"]:
                            variant_found = True

                if variant_found == True:        
                    sample_scores[project]["target_mutations"][sample_index] += 1

                    
            for col in params["info"]["extended_features"]:
                if col[len(col)-2:] != "_n":
                    selected_case     = extended_features[extended_features["ID:case id"] == case_ids[project][j]]

                    selected_features = [float(selected_case.iloc[l].loc[col]) for l in range(selected_case.shape[0])
                                         if selected_case.iloc[l].loc[col] != "-" and selected_case.iloc[l].loc[col] != "#NV"]

                    sample_scores[project][col+"_n"][sample_index] = len(selected_features)
                    
                    if col not in params["stats_exceptions"]:
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = np.mean(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None

                    elif params["stats_exceptions"][col] == "avg_wo_exception":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = np.mean(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = 0

                    elif params["stats_exceptions"][col] == "count_by_threshold":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = _count_by_threshold(selected_features, params["score_threshold"][col])
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None

                    elif params["stats_exceptions"][col] == "mean_by_threshold":
                        selected_features = [selected_feature for selected_feature in selected_features if selected_feature > 0]
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = np.mean(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None

                    elif params["stats_exceptions"][col] == "relative_count_by_threshold":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = _count_by_threshold(selected_features, params["score_threshold"][col]) / len(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None

                    elif params["stats_exceptions"][col] == "sum":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = np.sum(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None

                    elif params["stats_exceptions"][col] == "sum_wo_exception":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = np.sum(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = 0

                    elif params["stats_exceptions"][col] == "test":
                        if len(selected_features) >= 1: sample_scores[project][col][sample_index] = len(selected_features)
                        if len(selected_features) < 1:  sample_scores[project][col][sample_index] = None
                    
            # test for consistency of aggregated PTC variant count vs. count of PTC variant entries per sample
            if ("Start_Position" in sample_scores[project] and "Start_Position" in sample_scores[project]
                and sample_scores[project]["ptc_mutations"][sample_index] != sample_scores[project]["Start_Position"][sample_index]):
                print("< inconsistent sizes between aggregated PTC variant count and extracted PTC variant count @", project, sample_index, case_ids[project][sample_index],
                      sample_scores[project]["ptc_mutations"][sample_index], "/", sample_scores[project]["Start_Position"][sample_index])
                
            #else:
            #print(project, sample_index, case_ids[project][sample_index], sample_scores[project]["ptc_mutations"][sample_index], "/", sample_scores[project]["LABEL:NMD score"][sample_index])

            sample_index += 1

        last_project = project
        bar.next()
    bar.finish()

    if last_project != None and sample_index != len(sample_scores[project][params["info"]["rna"][0]]):
        print("< sample index count inconsistent with input sample size @", project, sample_index, "/", len(sample_scores[project][params["info"]["rna"][0]]))

    print("<", filtered_out, "items filtered out.")

    # determine samples to be filtered out
    sample_filter = {project: [] for project in projects}

    if params["filter_samples"] == True:
        keys = list(params["sample_filter"].keys())
        for i in range(len(keys)):
            temp = {project: [] for project in projects}
            [temp[project].extend([sample for sample in sample_scores[project][keys[i]]
                                          if (params["sample_filter"][keys[i]]["inclusive"] == True and sample_scores[project][keys[i]][sample] >= params["sample_filter"][keys[i]]["range"][0]
                                              and sample_scores[project][keys[i]][sample] < params["sample_filter"][keys[i]]["range"][1])
                                              or (params["sample_filter"][keys[i]]["inclusive"] == False and sample_scores[project][keys[i]][sample] > params["sample_filter"][keys[i]]["range"][0]
                                              and sample_scores[project][keys[i]][sample] < params["sample_filter"][keys[i]]["range"][1])
                                              or (params["sample_filter"][keys[i]]["include_misses"] == True and pd.isna(sample_scores[project][keys[i]][sample]) == True)])
                                              for project in projects]
            
            if i == 0:
                sample_filter = temp

            if i > 0:
                for project in projects:
                    sample_filter[project] = collections.Counter(sample_filter[project]) & collections.Counter(temp[project])

    #print(json.dumps(sample_filter["TCGA-ACC"], indent=4))
    sample_scores = {project: {**{rna_value         : {sample: {"median": np.median(sample_scores[project][rna_value][sample]),
                                                       "std":    np.std(sample_scores[project][rna_value][sample]),
                                                       "count":  len(sample_scores[project][rna_value][sample])}
                                                       for sample in sample_scores[project][rna_value] if params["filter_samples"] == False or sample in sample_filter[project]}
                                                       for rna_value in params["info"]["rna"]},
                               **{col               : {sample: {"median": sample_scores[project][col][sample]} for sample in sample_scores[project][col]
                                                       if params["filter_samples"] == False or sample in sample_filter[project]}
                                                       for col in params["info"]["appended_features"]},
                               **{col               : {sample: {"median": sample_scores[project][col][sample]} for sample in sample_scores[project][col]
                                                       if params["filter_samples"] == False or sample in sample_filter[project]}
                                                       for col in params["info"]["extended_features"]},
                               **{col               : {sample: {"median": sample_scores[project][col][sample]} for sample in sample_scores[project][col]
                                                       if params["filter_samples"] == False or sample in sample_filter[project]}
                                                       for col in params["info"]["mutations"]},
                               **{"sample_names"    : {sample: sample_scores[project]["sample_names"][sample] for sample in range(len(sample_scores[project]["sample_names"]))
                                                       if params["filter_samples"] == False or sample in sample_filter[project]}},
                               **{"target_mutations": {sample: {"median": sample_scores[project]["target_mutations"][sample]} for sample in sample_scores[project]["target_mutations"]
                                                       if params["filter_samples"] == False or sample in sample_filter[project]}}
                               } for project in projects}
    

    #print(json.dumps(sample_scores["TCGA-ACC"], indent=4))
    
    # test for inconsistent sizes of extracted targets
    for project in filtered_index:
        for rna_value in filtered_index[project]:
            for sample in sample_scores[project][rna_value]:
                if len(filtered_index[project][rna_value]) < sample_scores[project][rna_value][sample]["count"]:
                    print("< extracted target size inconsistent with filtered target size @", project, rna_value, sample,
                          "with", len(filtered_index[project][rna_value]), "/", sample_scores[project][rna_value][sample]["count"])


    if params["append_sample_scores"] == True: append_sample_scores(extended_features, sample_scores, case_ids, params)

    # calculate cancer scores
    cancer_scores_summary = {"project": [],
                             params["info"]["rna"][0]+"_mean": [],  params["info"]["rna"][0]+"_std": [],
                             params["info"]["rna"][0]+"_median": [],  params["info"]["rna"][0]+"_std": [],
                             params["info"]["rna"][0]+"_count": [], params["info"]["rna"][0]+"_raw": [],
                             **{col+"_mean":           [] for col in params["info"]["appended_features"]},
                             **{col+"_median":         [] for col in params["info"]["appended_features"]},
                             **{col+"_std":            [] for col in params["info"]["appended_features"]},
                             **{col+"_count":          [] for col in params["info"]["appended_features"]},
                             **{col+"_raw":            [] for col in params["info"]["appended_features"]},
                             **{col+"_mean":           [] for col in params["info"]["extended_features"]},
                             **{col+"_median":         [] for col in params["info"]["extended_features"]},
                             **{col+"_std":            [] for col in params["info"]["extended_features"]},
                             **{col+"_count":          [] for col in params["info"]["extended_features"]},
                             **{col+"_raw":            [] for col in params["info"]["extended_features"]},
                             **{col+"_mean":           [] for col in params["info"]["mutations"]},
                             **{col+"_median":         [] for col in params["info"]["mutations"]},
                             **{col+"_std":            [] for col in params["info"]["mutations"]},
                             **{col+"_count":          [] for col in params["info"]["mutations"]},
                             **{col+"_raw":            [] for col in params["info"]["mutations"]},
                             **{"sample_names":        []},
                             **{"target_mutations_mean":   []},
                             **{"target_mutations_median": []},
                             **{"target_mutations_std":    []},
                             **{"target_mutations_count":  []},
                             **{"target_mutations_raw":    []}}

    bar = IncrementalBar(set_bar("cancer scores are calculated"), max=len(sample_scores))

    for project in sample_scores:
        for col in sample_scores[project]:
            if col != "sample_names":
                score_medians = [float(sample_scores[project][col][sample]["median"]) for sample in sample_scores[project][col]
                                    if pd.isna(sample_scores[project][col][sample]["median"]) == False]
                cancer_scores_summary[col+"_raw"].append([sample_scores[project][col][sample]["median"] for sample in sample_scores[project][col]])

                if col == params["info"]["rna"][0]:
                    cancer_scores_summary["project"].append(project)

                if len(score_medians) > 0:
                    cancer_scores_summary[col+"_mean"].append(np.mean(score_medians))
                    cancer_scores_summary[col+"_median"].append(np.median(score_medians))
                    cancer_scores_summary[col+"_std"].append(np.std(score_medians))
                    cancer_scores_summary[col+"_count"].append(len(score_medians))

                else:
                    cancer_scores_summary[col+"_mean"].append(None)
                    cancer_scores_summary[col+"_median"].append(None)
                    cancer_scores_summary[col+"_std"].append(None)
                    cancer_scores_summary[col+"_count"].append(None)

            else: 
                sample_names = [sample_scores[project]["sample_names"][sample] for sample in sample_scores[project][col]]
                cancer_scores_summary["sample_names"].append(sample_names)

            if col == params["info"]["rna"][0]:
                length_test = [_ for sample in sample_scores[project][col] if sample_scores[project][col][sample]["count"] > len(target_ids)]
                if len(length_test) > 0: print("< error @calculate_target_expression. inconsistent number of target ids.")
        
        bar.next()
    bar.finish()

    cancer_scores_summary = pd.DataFrame(cancer_scores_summary).sort_values(by=[params["info"]["rna"][0]+"_mean"])
    cancer_scores_summary.to_csv(path_or_buf=params["newdir"]+params["os_sep"]+"cancer_scores"+params["file_tag"], sep=",", index=False)

    fig, ax = plt.subplots()
    ax.set_ylabel("cancer score")
    get_box_plot(ax, cancer_scores_summary, "project", params["info"]["rna"][0]+"_raw")
    #plt.show()
    fig.savefig(params["newdir"] + params["os_sep"] + "boxplot" + params["file_tag"] + ".png", dpi=600, bbox_inches='tight')

    mutations_mean_log10 = []
    mutations_raw_log10  = []
    for i in range(cancer_scores_summary.shape[0]):
        log10_scaled = [math.log10(mutation) for mutation in cancer_scores_summary.iloc[i].loc["total_raw"] if mutation > 0] # ptc_mutations_raw
        mutations_mean_log10.append(np.mean(log10_scaled))
        mutations_raw_log10.append(log10_scaled)
    
    cancer_scores_summary["total_mean_log10"] = mutations_mean_log10
    cancer_scores_summary["total_raw_log10"]  = mutations_raw_log10
    cancer_scores_summary                     = cancer_scores_summary.sort_values(by=["total_mean_log10"])

    fig, ax = plt.subplots()
    ax.set_ylabel("log10 mutations")
    ax = get_box_plot(ax, cancer_scores_summary, "project", "total_raw_log10")
    #plt.show()
    fig.savefig(params["newdir"] + params["os_sep"] + "mutations_plot" + params["file_tag"] + ".png", dpi=600, bbox_inches='tight')
    
    cancer_scores_summary.drop(columns=["total_mean_log10", "total_raw_log10"])
    cancer_scores_summary.to_csv(path_or_buf=params["newdir"]+params["os_sep"]+"cancer_scores"+params["file_tag"], sep=",", index=False)
    
    selected_cancer_scores_summary = create_selection(cancer_scores_summary, params["info"]["rna"][0]+"_raw")
    selected_cancer_scores_summary.to_csv(path_or_buf=params["newdir"]+params["os_sep"]+"transformed_cancer_scores"+params["file_tag"], sep=",", index=False)

    selected_cancer_scores_summary = create_selection(cancer_scores_summary, "exon_raw")
    selected_cancer_scores_summary.to_csv(path_or_buf=params["newdir"]+params["os_sep"]+"transformed_exon_mutations"+params["file_tag"], sep=",", index=False)

    selected_cancer_scores_summary = create_selection(cancer_scores_summary, "total_raw")
    selected_cancer_scores_summary.to_csv(path_or_buf=params["newdir"]+params["os_sep"]+"transformed_total_mutations"+params["file_tag"], sep=",", index=False)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', None)
    plt.close("all")
    return targets


def contains_number(text):
    numbers = [i for i in range(len(text)) if text[i].isnumeric()]
    if len(numbers) > 0: return True
    else:                return False


def create_newdir(params):
    datestr = str(datetime.now())
    index   = datestr.index(".")
    dirname = datestr[:index].replace(" ", "_").replace(":", "-")
    newdir  = params["run_dir"] + params["os_sep"] + dirname + params["file_tag"]
    print("< new directory:", newdir)
    if not os.path.exists(newdir): os.mkdir(newdir)
    params["newdir"] = newdir
    return


def create_selection(cancer_scores_summary, target_col):
    max_size = None
    for i in range(cancer_scores_summary.shape[0]):
        if max_size == None or len(cancer_scores_summary.iloc[i].loc[target_col]) > max_size:
            max_size = len(cancer_scores_summary.iloc[i].loc[target_col])

    selected_cancer_scores_summary = pd.DataFrame({**{"project": [cancer_scores_summary.iloc[i].loc["project"] for i in range(cancer_scores_summary.shape[0])]},
                                                   **{"sample score " + str(i): ["NA" for _ in range(cancer_scores_summary.shape[0])] for i in range(max_size)}})
    for i in range(cancer_scores_summary.shape[0]):
        for j in range(len(cancer_scores_summary.iloc[i].loc[target_col])):
            selected_cancer_scores_summary.at[selected_cancer_scores_summary.index[i], "sample score " + str(j)] = cancer_scores_summary.iloc[i].loc[target_col][j]

    return selected_cancer_scores_summary


def extract_variants(cluster, data, project_key, cluster_key, params, lock):
    project = cluster[project_key][cluster_key]["project"]

    for i in range(len(cluster[project_key][cluster_key]["input"])):
        for case_id in cluster[project_key][cluster_key]["input"][i]:
            if params["filter_mutations"] == True:
                rna_mean = pd.DataFrame()
                rnas = []
                for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])):
                    if params["transform_type"] == None:
                        rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"RNA"
                                       +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)
        
                    else:
                        rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+params["transform_type"]
                                       +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)
                    
                    if rna.shape[0] > 0: rnas.append(rna)

                if len(rnas) > 0:
                    temp_data = pd.DataFrame({**{params["target_identifier"]["rna"]: rnas[j][params["target_identifier"]["rna"]] for j in range(len(rnas)) if rnas[j].shape[0] > 0},
                                              **{j: rnas[j][params["info"]["rna"][0]] for j in range(len(rnas)) if rnas[j].shape[0] > 0}})

                    # if multiple rna seq files exist, calculate averages first for each gene
                    rna_mean       = temp_data[[j for j in range(len(rnas))]].mean(axis=1)

                    # remove version number if not gene symbol (added 250317)
                    if params["target_identifier"]["rna"] != "gene_name":
                        rna_mean.index = [rnas[0].iloc[j].loc[params["target_identifier"]["rna"]].split(".")[0] for j in range(rnas[0].shape[0])]

                    else:
                        rna_mean.index = rnas[params["target_identifier"]["rna"]]

                    if params["filter_mutations_type"] == "adjusted":
                        adjusted_filter = rna_mean[rna_mean > 0].sort_values()[int(rna_mean[rna_mean > 0].shape[0]*params["filter_mutation_margin"])]
                        rna_mean        = rna_mean[rna_mean > adjusted_filter]

                    if params["filter_mutations_type"] == "threshold":
                        rna_mean = rna_mean[rna_mean > params["rna_threshold"]]


            # load all whole exome sequencing files for case id 
            wxss = []
            for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"])):
                wxs = load_wxs(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"WXS"+params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][j])

                if params["filter_mutations"] == True:
                    wxs = wxs[wxs[params["target_identifier"]["wxs"]].isin(rna_mean.index.tolist())]

                wxss.append(wxs)


            wxs_dict = {**{col: [[] for _ in range(len(wxss))] for col in params["info"]["wxs"]},
                        **{col: [[] for _ in range(len(wxss))] for col in params["info"]["mutations"]}}
            
            for j in range(data.shape[0]):
                # medium layer (entries in different wxs files)
                variant_classifications = [[] for _ in range(len(wxss))]
                variant_ids             = [[] for _ in range(len(wxss))]

                # clear target-specific mutation data (as defined in info-wxs), keep non-specific mutation info
                for col in params["info"]["wxs"]:
                    wxs_dict[col] = [[] for _ in range(len(wxss))]

                for k in range(len(wxss)):
                    if len(wxss[k].columns.tolist()) > 0:
                        # extract general information once (equal for all targets)
                        if j == 0:
                            # check if mutations other than specified are present
                            #if wxss[k][wxss[k]["Variant_Classification"].isin(params["info"]["mutations"])].shape[0] != wxss[k].shape[0]:
                            #    print("< unconsidered mutational categories detetected @", project_key, cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][k])
                            #    print(wxss[k][~wxss[k]["Variant_Classification"].isin(params["info"]["mutations"])].drop_duplicates(subset=["Variant_Classification"])["Variant_Classification"].tolist())

                            # get general mutation info
                            wxs_dict["total"][k].append(wxss[k].shape[0])
                            for col in params["info"]["mutations"]:
                                if col != "total" and col != "ptc_mutations":
                                    wxs_dict[col][k].append(wxss[k][wxss[k]["Variant_Classification"] == col].shape[0])

                            # conduct check for active PTCs (not all PTC mutations must in fact be active)
                            wxs_dict = get_active_ptcs(wxs_dict, wxss[k], params, project, cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][k], k)

                        # create selection of mutations for targets
                        target_variants = wxss[k][wxss[k][params["target_identifier"]["wxs"]] == data.iloc[j].loc[params["target_identifier"]["rna"]]]
                        
                        # inner layer (multiple entries for a single wxs file)
                        for l in range(target_variants.shape[0]):
                            variant_classifications[k].append(target_variants.iloc[l].loc["Variant_Classification"])
                                                               
                            # create variant id as a unique identifier for the variant type, changed 241118
                            if "HGVSp" in target_variants.columns.tolist() and type(target_variants.iloc[l].loc["HGVSp"]) == str:                     
                                variant_ids[k].append(target_variants.iloc[l].loc[params["wxs_identifier"]] + "_" + target_variants.iloc[l].loc["HGVSp"])
                                #variant_ids[k].append(target_variants.iloc[l].loc["Gene"] + "_" + target_variants.iloc[l].loc["HGVSp"])

                            elif "HGVSp_Short" in target_variants.columns.tolist() and type(target_variants.iloc[l].loc["HGVSp_Short"]) == str:
                                variant_ids[k].append(target_variants.iloc[l].loc[params["wxs_identifier"]] + "_" + target_variants.iloc[l].loc["HGVSp_Short"])
                                #variant_ids[k].append(target_variants.iloc[l].loc["Gene"] + "_" + target_variants.iloc[l].loc["HGVSp_Short"])

                            elif "HGVSc" in target_variants.columns.tolist() and type(target_variants.iloc[l].loc["HGVSc"]) == str:
                                variant_ids[k].append(target_variants.iloc[l].loc[params["wxs_identifier"]] + "_" + target_variants.iloc[l].loc["HGVSc"])
                                #variant_ids[k].append(target_variants.iloc[l].loc["Gene"] + "_" + target_variants.iloc[l].loc["HGVSc"])

                            elif "HGVSp" in target_variants.columns.tolist(): # conditions added to exclude SCLC data which lack according data columns
                                print("< variant id could not be created for target:", target_variants.iloc[l].loc[params["wxs_identifier"]],
                                        " with HGVSp:",  target_variants.iloc[l].loc["HGVSp"], ", and HGVSp_Short:", target_variants.iloc[l].loc["HGVSp_Short"],
                                        ", and HGVSc:", target_variants.iloc[l].loc["HGVSc"])
                                pd.set_option('display.max_columns', None)
                                print(target_variants)

                            for key in params["info"]["wxs"]:
                                if key in target_variants.columns: wxs_dict[key][k].append(target_variants.iloc[l].loc[key])

                # add information to data matrix
                data.at[data.iloc[j].loc[params["target_identifier"]["rna"]], "variant_classifications"].append(variant_classifications)
                data.at[data.iloc[j].loc[params["target_identifier"]["rna"]], "variant_ids"].append(variant_ids)

                for key in wxs_dict:
                    data.at[data.iloc[j].loc[params["target_identifier"]["rna"]], key].append(wxs_dict[key])


    test_cols = ["variant_ids"]
    test_dimensions(data, test_cols, len(data.iloc[0].loc[params["info"]["rna"][0]]), description=project_key+"_"+cluster_key, function_name="extract_variants")
    return data


def extract_non_targets_rnaseq_(data, rnas, non_targets_reference, params, project, cluster_key, rna_id):
    for rnacol in params["info"]["rna"]:
        temp_data = pd.DataFrame({**{params["target_identifier"]["rna"]: rnas[i][params["target_identifier"]["rna"]] for i in range(len(rnas)) if rnas[i].shape[0] > 0},
                                  **{i: rnas[i][rnacol] for i in range(len(rnas)) if rnas[i].shape[0] > 0}})

        # if multiple rna seq files exist, calculate averages first for each gene
        if temp_data.shape[0] > 0:
            # create non-targets by excluding targets from expression genes
            if params["non_targets_explicit"] == False: non_targets = temp_data[~temp_data[params["target_identifier"]["rna"]].isin(data[params["target_identifier"]["rna"]])]
            if params["non_targets_explicit"] == True:  non_targets = temp_data[temp_data[params["target_identifier"]["rna"]].isin(non_targets_reference[params["target_identifier"]["rna"]])]

            if params["non_targets_explicit"] == False and non_targets.shape[0] != rnas[0].shape[0]-params["target_size"]:
                print("< error occurred @extract_non_targets_rnaseq_:", non_targets.shape[0], "/", rnas[0].shape[0]-params["target_size"])

            # filter out any zeros before averaging
            # changed on 250317 to exclude '-1' placeholders from apply_cnv
            non_targets = non_targets[~(non_targets == -1).any(axis=1)]
            non_targets = non_targets[~(non_targets == 0).any(axis=1)]

            mean        = non_targets[[i for i in range(len(rnas))]].mean(axis=1)
            [data.at[i, rnacol + "_non_targets"].append(mean.median()) for i in data.index]
            if temp_data.shape[1] >= 3:
                print("extract_non_targets_rnaseq_", rna_id, mean.shape, non_targets.shape, mean.median())
            
            #print("extract_non_targets_rnaseq_2", rna_id, temp_data.shape, non_targets.shape, mean.median())

            if non_targets.shape[1]-1 != len(rnas):
                print("< dimension error occurred @extract_non_targets_rnaseq_ with:", non_targets.shape[1]-1, "/", len(rnas))
                print(" ", project, cluster_key, rna_id)
        
        else:
            [data.at[i, rnacol + "_non_targets"].append(0) for i in data.index]

    return data


def extract_targets_rnaseq_(data, rnas, params, project, cluster_key, case_id, rna_id):
    # select targets from expression data
    for i in range(len(rnas)):
        if rnas[i].shape[0] > 0:
            rnas[i] = rnas[i][rnas[i][params["target_identifier"]["rna"]].isin(data[params["target_identifier"]["rna"]])]

            # remove version number if not gene symbol (added 250317)
            if params["target_identifier"]["rna"] != "gene_name":
                rnas[i].index = [rnas[i].iloc[j].loc[params["target_identifier"]["rna"]].split(".")[0] for j in range(rnas[i].shape[0])]

            else:
                rnas[i].index = rnas[i][params["target_identifier"]["rna"]]


            if rnas[i].shape[0] != params["target_size"]:
                print("< missing targets @extract_targets_rnaseq_:", rnas[i].shape[0], "/", params["target_size"])
                print(rnas[i])
                input("X")

    for i in range(data.shape[0]):
        # register case id
        data.at[data.iloc[i].loc[params["target_identifier"]["rna"]], "case_ids"].append(case_id)

        multiple_seqs    = {rnacol: [] for rnacol in params["info"]["rna"]}
        target_not_found = False

        # medium layer (if multiple expression data exist for the same case id)
        for j in range(len(rnas)):
            overlapping_seqs = {rnacol: [] for rnacol in params["info"]["rna"]}

            if data.iloc[i].loc[params["target_identifier"]["rna"]] in rnas[j].index.tolist():
                for rnacol in params["info"]["rna"]:
                    # inner layer for hypothetical ambiguous entries is not used in this context but preserved for simplicity
                    overlapping_seqs[rnacol].append(rnas[j].loc[data.iloc[i].loc[params["target_identifier"]["rna"]]].loc[rnacol])
                
            else:
                target_not_found = True

            if target_not_found == True and rnas[j].shape[0] > 0:
                print("<", data.iloc[i].loc[params["target_identifier"]["rna"]], "not found.")
                print(" ", project, cluster_key, rna_id)
                print(rnas[j])

            # write info to multiple seqs dictionary, outer layer accounts for multiple rna seq datasets
            for rnacol in params["info"]["rna"]:
                multiple_seqs[rnacol].append(overlapping_seqs[rnacol])

        for rnacol in multiple_seqs:
            data.at[data.iloc[i].loc[params["target_identifier"]["rna"]], rnacol].append(multiple_seqs[rnacol])

    return data


def extract_targets_rnaseq(cluster, data, non_targets, project_key, cluster_key, params, lock):
    project     = cluster[project_key][cluster_key]["project"]
    loaded_rnas = 0

    for i in range(len(cluster[project_key][cluster_key]["input"])):
        rna_loaded = 0
        for case_id in cluster[project_key][cluster_key]["input"][i]:
            rnas = []
            
            for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])):
                if params["transform_type"] == None:
                    rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"RNA"
                                   +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)
    
                else:
                    rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+params["transform_type"]
                                   +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)
                
                # marked (<-) added / removed on 250811
                # rnas.append(rna) # <- removed
                if rna.shape[0] > 0: # <- added
                    rnas.append(rna) # <- added

                rna_loaded = 1
            
            # marked (<-) added / removed on 250811
            #try: # <- removed
            #    if params["target_normalization"] == True: data = extract_non_targets_rnaseq_(data, rnas, non_targets, params, project, cluster_key,
            #                                                                                  cluster[project_key][cluster_key]["input"][i][case_id]["RNA"]) # <- removed
            #except: # <- removed
            #    print("normalization error @project", project_key, "rna files", cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])
            #    for rna in rnas:
            #        print(rna) # <- removed

            # if-condition added
            if params["target_normalization"] == True:
                data = extract_non_targets_rnaseq_(data, rnas, non_targets, params, project, cluster_key,
                                                   cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])

            data = extract_targets_rnaseq_(data, rnas, params, project, cluster_key, case_id, cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])

        loaded_rnas += rna_loaded
    
    test_dimensions(data, params["info"]["rna"], loaded_rnas, description=project_key+"_"+cluster_key, function_name="extract_targets_rnaseq")
    return data


def extract_targets(cluster, clusters, targets, non_targets, proc_index, params, lock, thread_id):
    for i in range(proc_index[0], proc_index[1]+1, 1):
        project_key = clusters[i]["project_key"]
        cluster_key = clusters[i]["cluster_key"]
        project     = cluster[project_key][cluster_key]["project"]

        # initialize data container
        data = pd.DataFrame({**{"project":                 [project_key for _ in range(targets.shape[0])]},
                             **{"cluster_key":             [cluster_key for _ in range(targets.shape[0])]},
                             **{"case_ids":                [[] for _ in range(targets.shape[0])]},
                             **{"variant_classifications": [[] for _ in range(targets.shape[0])]},
                             **{"variant_ids":             [[] for _ in range(targets.shape[0])]},
                             **{col:                       [targets.iloc[j].loc[col] for j in range(targets.shape[0])] for col in targets.columns},
                             **{col:                       [[] for _ in range(targets.shape[0])] for col in params["info"]["wxs"]}, # stores mutants per sample for all targets
                             **{col:                       [[] for _ in range(targets.shape[0])] for col in params["info"]["mutations"]},
                             **{col:                       [[] for _ in range(targets.shape[0])] for col in params["info"]["rna"]}, # stores values per sample for all targets
                             **{col+"_non_targets":        [[] for _ in range(targets.shape[0])] for col in params["info"]["rna"]}}) # stores median values per sample for all non-targets
        #                        index = [targets.iloc[j].loc[params["target_identifier"]["rna"]].split(".")[0] for j in range(targets.shape[0])])

        # remove version number if not gene symbol (added 250317)
        if params["target_identifier"]["rna"] != "gene_name":
            data.index = [targets.iloc[j].loc[params["target_identifier"]["rna"]].split(".")[0] for j in range(targets.shape[0])]

        else:
            data.index = targets[params["target_identifier"]["rna"]]
            
        if os.path.isdir(params["data_dir"]+params["os_sep"]+project) == True:
            if os.path.isdir(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]) == False:
                os.mkdir(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"])

            if os.path.isfile(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]+params["os_sep"]+project_key+"_"+cluster_key+".txt") == False:
                print("< thread:", thread_id, "project:", project_key, "cluster:", cluster_key)
                data = extract_targets_rnaseq(cluster, data, non_targets, project_key, cluster_key, params, lock)

                # for implicit target types, no mutational information is added (due to performance issues)
                if params["track_mutations"] == True: data = extract_variants(cluster, data, project_key, cluster_key, params, lock)

                # store data
                data.to_csv(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]+params["os_sep"]+project_key+"_"+cluster_key+".txt", index=False, sep=",")

        else:
            print("<", params["data_dir"]+params["os_sep"]+project, " does not exist.")

        del data

    return


# remove 250313
'''
def filter_variants(variants, params):
    init_shape = variants.shape[0]
    for key in params["variant_value_filter"]:
        if key not in variants.columns:
            print("< error. filter item", key, "not found in variants.")

        else:
            # if dtype is object, remove skip tag first ("-"), then covert to float
            if variants[key].dtype == "object":
                variants = variants[~variants[key].isna()]
                variants = variants.astype({key: "float32"})

            variants = variants[variants[key] >= params["variant_value_filter"][key]]
    
    print("< application of variant filter reduced variants from", init_shape, "to", variants.shape[0])
    return variants
'''


def get_box_plot(ax, data, x_col, y_col, y_label="cancer score"):
    cmap = plt.get_cmap('viridis')
    for i in range(data.shape[0]):
        #print("i", i, data.iloc[i].loc[x_col], len(data.iloc[i].loc[y_col]), np.median(data.iloc[i].loc[y_col]))
        if data.iloc[i].loc[x_col] == "TCGA-KIRP":
            for j in range(len(data.iloc[i].loc[y_col])):
                if pd.isna(data.iloc[i].loc[y_col][j]) == True:
                    print(i, j, data.iloc[i].loc["sample_names"][j])
        #print(data.iloc[i].loc[y_col])
        boxplot = ax.boxplot(data.iloc[i].loc[y_col], labels=[data.iloc[i].loc[x_col]],
                             meanline=True, notch=True, positions=[i], showmeans=True, patch_artist=True, widths=0.5)

        for patch in zip(boxplot['boxes']):
            patch[0].set_facecolor(cmap(float(i)/data.shape[0]))

        for line in zip(boxplot['caps']):
            line[0].set_linewidth(3)

        for line in zip(boxplot['medians']):
            line[0].set_linewidth(3)
            line[0].set_color("dimgray")

        for line in zip(boxplot['means']):
            line[0].set_linestyle("--")
            line[0].set_linewidth(3)
            line[0].set_color("black")

        for line in zip(boxplot['fliers']):
            line[0].set_linewidth(1)
            line[0].set_markeredgecolor(cmap(float(i)/data.shape[0]))

        for line in zip(boxplot['whiskers']):
            line[0].set_linewidth(2)

    ax.tick_params(axis='x', labelrotation=90, labelsize=14)
    ax.tick_params(axis='y', labelrotation=90, labelsize=14)
    ax.set_ylabel(y_label, fontsize=18)
    #plt.xticks(rotation=90)
    return ax


def get_active_ptcs(wxs_dict, wxs, params, project, wxs_fname, it):
    ptcs        = create_ptc_subset(wxs, project, wxs_fname)
    active_ptcs = 0
    variant_ids = []

    for i in range(ptcs.shape[0]):
        # check whether mutation leads to elongated protein ("ext")
        #if "HGVSp" in ptcs.columns.tolist() and (len(ptcs.iloc[i].loc["HGVSp"]) > 0 and "ext" not in ptcs.iloc[i].loc["HGVSp"] and "=" not in ptcs.iloc[i].loc["HGVSp"] and "?" not in ptcs.iloc[i].loc["HGVSp"]):
        if ("HGVSp" in ptcs.columns.tolist() and (len(ptcs.iloc[i].loc["HGVSp"]) > 0 and "ext" not in ptcs.iloc[i].loc["HGVSp"]
                                                  and "=" not in ptcs.iloc[i].loc["HGVSp"] and "?" not in ptcs.iloc[i].loc["HGVSp"])
            and ptcs.iloc[i].loc["Variant_Classification"] in ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"]):
            # check whether the gene is present multiple times, if so, only the most downstream ptc variant is selected
            if check_variants(wxs, ptcs, params, i) == True:
                # get ptc position
                ptc_position = get_ptc(ptcs.iloc[i].loc["HGVSp"])

                # remove duplicates (due to mutations at different positions with same position of the resulting PTC), changed 241118
                if ptc_position != None and ptcs.iloc[i].loc[params["wxs_identifier"]] + "_" + str(ptc_position) not in variant_ids:
                    active_ptcs += 1
                    variant_ids.append(ptcs.iloc[i].loc[params["wxs_identifier"]] + "_" + str(ptc_position))

    wxs_dict["ptc_mutations"][it].append(active_ptcs)
    return wxs_dict


def _append_features(df, data, params, it, step, text):
    mapped_index = append_df_with_mapping(data, "case id", "case id", "index", text)
    tumor_stages = {"Stage I": 0, "Stage II": 1, "Stage III": 2, "Stage IV": 3}

    for feature in params["info"]["appended_features"]:
        if feature in data[1].columns:
            values = []
            for i in mapped_index:
                if i != "-":
                    if data[1].iloc[int(i)].loc["case id"] != data[0].iloc[len(values)]["case id"]:
                        print("< mapping error @append_features")
                        exit()

                    if data[1].iloc[int(i)].loc[feature] != "#NV": # test for exception of immune dataset
                        if feature != "ajcc_pathologic_tumor_stage":
                            values.append(data[1].iloc[int(i)].loc[feature]) 
                        
                        if feature == "ajcc_pathologic_tumor_stage":
                            if data[1].iloc[int(i)].loc[feature] in tumor_stages: values.append(tumor_stages[data[1].iloc[int(i)].loc[feature]])
                            else:                                                 values.append(None)
                    
                    else:
                        values.append(None)
                
                # marked (<-) added / removed on 250605
                # elif "nmd mutations" not in text: # <- removed
                elif "hla mutations" not in text and "nmd mutations" not in text: # <- added
                    values.append(None)

                else: # for nmd mutations, missing values mean that no mutations were found for the respective case id
                    values.append(0) 

            if len(values) == data[0].shape[0]:
                for i in range(it, it+step):
                    df.at[df.index[i], feature] = values

            else:
                print("< dimension error @append_features")
                exit()

    return


def append_features(df, params):
    # load clinical data
    clinical_data           = pd.read_csv(params["clinical_data_path"], delimiter=",")
    clinical_data["index"]  = [i for i in range(clinical_data.shape[0])]
    check_placeholders(clinical_data, [col for col in params["info"]["appended_features"] if col in clinical_data.columns])

    # load immune data
    immune_data             = pd.read_csv(params["immune_data_path"], delimiter=",")
    immune_data["index"]    = [i for i in range(immune_data.shape[0])]
    check_placeholders(immune_data, [col for col in params["info"]["appended_features"] if col in immune_data.columns])

    # load hla mutations
    hla_mutations           = pd.read_csv(params["hla_mutations_path"], delimiter=",")
    hla_mutations["index"]  = [i for i in range(hla_mutations.shape[0])]
    check_placeholders(hla_mutations, [col for col in params["info"]["appended_features"] if col in hla_mutations.columns])

    # load nmd mutations
    nmd_mutations           = pd.read_csv(params["nmd_mutations_path"], delimiter=",")
    nmd_mutations["index"]  = [i for i in range(nmd_mutations.shape[0])]
    check_placeholders(nmd_mutations, [col for col in params["info"]["appended_features"] if col in nmd_mutations.columns])

    # load immune editing data
    immune_editing          = pd.read_csv(params["immune_editing_path"], delimiter=",")
    immune_editing["index"] = [i for i in range(immune_editing.shape[0])]
    check_placeholders(immune_editing, [col for col in params["info"]["appended_features"] if col in immune_editing.columns])


    # create containers for features to be appended
    for col in params["info"]["appended_features"]:
        df.insert(df.shape[1], col, [None for _ in range(df.shape[0])])

    projects = df.drop_duplicates(subset="project")["project"]
    step     = df[df["project"] == projects[0]].shape[0]
    bar      = IncrementalBar(set_bar("appending features"), max=len(projects))
    for i in range(0, df.shape[0], step):
        projectwise_case_ids = pd.DataFrame({"case id": df.iloc[i].loc["case_ids"]})
        #print("i", i, len(projectwise_case_ids), df.iloc[i].loc["project"])
        _append_features(df, [projectwise_case_ids, clinical_data], params, i, step, set_bar("mapping clinical data"))
        _append_features(df, [projectwise_case_ids, immune_data], params, i, step, set_bar("mapping immune data"))
        _append_features(df, [projectwise_case_ids, hla_mutations], params, i, step, set_bar("mapping hla mutations"))
        _append_features(df, [projectwise_case_ids, nmd_mutations], params, i, step, set_bar("mapping nmd mutations"))
        _append_features(df, [projectwise_case_ids, immune_editing], params, i, step, set_bar("mapping immune editing data"))
        bar.next()
    bar.finish()

    return df


def remove_misses(df, params):
    # apply changes to data
    check_df_index(df)
    last_project = None
    for i in range(df.shape[0]):
        for col in params["info"]["rna"]:
            targets                      = json.loads(df.iloc[i].loc[col])
            selected_index               = [j for j in range(len(targets)) if contains_number(json.dumps(targets[j])) == True]
            df.at[df.index[i], col]      = [targets[j] for j in selected_index]
            
            if df.iloc[i].loc["project"] != last_project:
                print("<", len(targets)-len(selected_index), "values removed for", df.iloc[i].loc["project"])

            non_targets                  = json.loads(df.iloc[i].loc[col+"_non_targets"])
            if len(non_targets) > 0:
                df.at[df.index[i], col+"_non_targets"] = [non_targets[j] for j in selected_index]
                if len(non_targets) != len(targets): print("size error1 @remove_misses:", len(non_targets), "/", len(targets))        

        for col in params["info"]["wxs"]:
            # nan values must be replaced with '-1' before json loads, then converted back to None
            try:
                samples                 = json.loads(df.iloc[i].loc[col].replace("\"\"", "").replace("5'Flank", "").replace("3'Flank", "").replace("5'UTR", "").replace("3'UTR", "").replace("'", "\"").replace("nan", "-1"))
                df.at[df.index[i], col] = [samples[j] if "-1" not in samples[j] else None for j in selected_index]
                if len(samples) != len(targets): print("size error2 @remove_misses:", len(samples), "/", len(targets))

            except:
                print("< format error occurred @remove_misses.")
                print(df.iloc[i].loc[col])
                exit()

        for col in params["info"]["mutations"]:
            if col in df.columns and col not in params["info"]["combined_mutations"]:
                if type(df.iloc[i].loc[col]) == str:
                    samples                 = json.loads(df.iloc[i].loc[col])
                    df.at[df.index[i], col] = [samples[j] for j in selected_index]
                    if len(samples) != len(targets): print("size error3 @remove_misses:", len(samples), "/", len(targets))

        for col in params["info"]["combined_mutations"]:
            for mutation_col in params["info"]["combined_mutations"][col]:
                if mutation_col in df.columns and type(df.iloc[i].loc[mutation_col]) == str:
                    samples                          = json.loads(df.iloc[i].loc[mutation_col])
                    df.at[df.index[i], mutation_col] = [samples[j] for j in selected_index]
                    if len(samples) != len(targets): print("size error3 @remove_misses:", len(samples), "/", len(targets))

        if "variant_classifications" in df.columns:
            # nan values must be replaced with '-1' before json loads, then converted back to None
            if type(df.iloc[i].loc["variant_classifications"]) == str:
                samples       = json.loads(df.iloc[i].loc["variant_classifications"].replace("\"\"", "").replace("5'Flank", "").replace("3'Flank", "").replace("5'UTR", "").replace("3'UTR", "").replace("'", "\"").replace("nan", "-1"))
                if len(samples) > 0: df.at[df.index[i], "variant_classifications"] = [samples[j] if "-1" not in samples[j] else None for j in selected_index]

        # changed on 24/11/11: remove misses from case ids
        case_ids                                            = json.loads(df.iloc[i].loc["case_ids"].replace("'", "\""))
        if len(samples) > 0: df.at[df.index[i], "case_ids"] = [case_ids[j] for j in selected_index]
        last_project                                        = df.iloc[i].loc["project"]

    # test for identical dimensions
    for i in range(df.shape[0]):
        test_dimensions(df.iloc[i:i+1], [*params["info"]["rna"], *params["info"]["wxs"]], len(df.iloc[i].loc[params["info"]["rna"][0]]), function_name="remove_misses")
    
    return df


def test_dimensions(data, target_cols, init_size, description=None, function_name="extract_targets_rnaseq"):
    failed = False
    # conduct dimension test
    last_size = None
    for i in range(data.shape[0]):
        if last_size != None and len(data.iloc[i].loc[target_cols[0]]) != last_size:
            failed = True
            if description != None: print("< size error1 occurred @", function_name, description, "with:", len(data.iloc[i].loc[target_cols[0]]), "/", last_size)
            if description == None: print("< size error1 occurred @", function_name, "with:", len(data.iloc[i].loc[target_cols[0]]), "/", last_size)

        last_current_size = None
        for target_col in target_cols:
            if last_current_size != None and len(data.iloc[i].loc[target_col]) != last_current_size:
                failed = True
                if description != None: print("< size error2 occurred @", function_name, description, "@ column: ", target_col,
                                              ", with:", len(data.iloc[i].loc[target_col]), "/", last_current_size)
                if description == None: print("< size error2 occurred @", function_name, "@ column: ", target_col, ", with:", len(data.iloc[i].loc[target_col]), "/", last_current_size)

            last_current_size = len(data.iloc[i].loc[target_col])

        last_size = last_current_size

    if last_size != init_size:
        if description != None: print("< error @", function_name, description, ", inconsistent sizes:", last_size, ", init size:", init_size); print(target_cols)
        if description == None: print("< error @", function_name, ", inconsistent sizes:", last_size, ", init size:", init_size); print(target_cols)

    return failed