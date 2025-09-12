from ast import literal_eval
import numpy as np
import os
import pandas as pd
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from extract_mutations_utils import *
from shared_utils import *


def _extract_expressions(values, stats, cluster, project_key, cluster_key, params, lock):
    project = cluster[project_key][cluster_key]["project"]

    for i in range(len(cluster[project_key][cluster_key]["input"])):
        # loop not necessary, only for convenience
        if len(cluster[project_key][cluster_key]["input"][i]) > 1:
            print("< error. inconsistent size @_extract_expressions")
            print(project_key, cluster_key)
            print(json.dumps(cluster[project_key][cluster_key]["input"][i], indent=4))
            exit()
        
        for case_id in cluster[project_key][cluster_key]["input"][i]:
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
            for j, rna_target in enumerate(params["rna_targets"]):
                temp_data = pd.DataFrame({**{j: rnas[j][rna_target].tolist() for j in range(len(rnas)) if rnas[j].shape[0] > 0}}, index=rnas[0][params["rna_identifier"]])

                # if multiple rna seq files exist, calculate averages first for each gene
                if temp_data.shape[0] > 0:
                    # filter out any zeros (or '-1' placeholders) before averaging
                    temp_data = temp_data[~(temp_data <= params["rna_threshold"]).any(axis=1)]
                    mean      = temp_data[[j for j in range(len(rnas))]].mean(axis=1)
                    [values[rna_target][mean.index[j]].append(mean.iloc[j]) for j in range(mean.shape[0])]

                    if j == 0:
                        stats[project_key+"_"+cluster_key]["averaged_cases"]  += 1
                        stats[project_key+"_"+cluster_key]["total_rna_files"] += temp_data.shape[1]

    return values, stats


def extract_expressions(cluster, clusters, template, proc_index, params, lock, thread_id):
    # initialize stats for report
    stats = {clusters[i]["project_key"]+"_"+clusters[i]["cluster_key"]: {"averaged_cases": 0, "total_rna_files": 0} for i in range(proc_index[0], proc_index[1]+1, 1)}

    for i in range(proc_index[0], proc_index[1]+1, 1):
        project_key = clusters[i]["project_key"]
        cluster_key = clusters[i]["cluster_key"]
        project     = cluster[project_key][cluster_key]["project"]

        # initialize data container
        data = pd.DataFrame({**{col: [template.iloc[j].loc[col] for j in range(template.shape[0])] for col in ["gene_id", "gene_name", "gene_type"]},
                             **{project+"_"+stats_type+"_"+col: [0 for _ in range(template.shape[0])]
                                for col in template.columns if col in params["rna_targets"]
                                for stats_type in ["mean", "median"]}},
                                index=template[params["rna_identifier"]])
        
        # added on 250213 because of dtype warning (explicit conversion required)
        data = data.astype({project+"_"+stats_type+"_"+col: 'float64'
                            for col in template.columns if col in params["rna_targets"]
                            for stats_type in ["mean", "median"]})
        
        values = {**{col: {template.iloc[j].loc[params["rna_identifier"]]: [] for j in range(template.shape[0])} for col in template.columns if col in params["rna_targets"]}}

        
        if os.path.isdir(params["data_dir"]+params["os_sep"]+project) == True:
            if os.path.isdir(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]) == False:
                os.mkdir(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"])

            if os.path.isfile(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]+params["os_sep"]+project_key+"_"+cluster_key+".txt") == False:
                print("< thread:", thread_id, "project:", project_key, "cluster:", cluster_key)
                values, stats = _extract_expressions(values, stats, cluster, project_key, cluster_key, params, lock)

                # calculate stats
                for rna_target in params["rna_targets"]:
                    for key in values[rna_target]:
                        if len(values[rna_target][key]) > 0:  data.at[key, project+"_mean_"+rna_target]   = np.mean(values[rna_target][key])
                        if len(values[rna_target][key]) == 0: data.at[key, project+"_mean_"+rna_target]   = None
                        if len(values[rna_target][key]) > 0:  data.at[key, project+"_median_"+rna_target] = np.median(values[rna_target][key])
                        if len(values[rna_target][key]) == 0: data.at[key, project+"_median_"+rna_target] = None

                # store data
                data.to_csv(params["data_dir"]+params["os_sep"]+"rnaseq"+params["file_tag"]+params["os_sep"]+project_key+"_"+cluster_key+".txt", index=False, sep=",")

        else:
            print("<", params["data_dir"]+params["os_sep"]+project, " does not exist.")
        
        del data

    print("< stats for process", thread_id)
    print(json.dumps(stats, indent=4))
    return