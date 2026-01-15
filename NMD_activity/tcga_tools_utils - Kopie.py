import gzip
import hashlib
import json
import numpy as np
import os
import pandas as pd
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from extract_mutations_utils import *
from shared_utils import *


def assign_cluster(x, y, intercept, slope):
    x1 = [x[i] for i in range(len(x)) if y[i] > x[i]*slope+intercept]
    y1 = [y[i] for i in range(len(y)) if y[i] > x[i]*slope+intercept]
    x2 = [x[i] for i in range(len(x)) if y[i] <= x[i]*slope+intercept]
    y2 = [y[i] for i in range(len(y)) if y[i] <= x[i]*slope+intercept]
    return x1, x2, y1, y2


def calculate_md5(path):
    hasher = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def compress(path, target_path):
    with open(path, "rb") as f:
        input = f.read()

    compressed_input = gzip.compress(input)

    with open(target_path, "wb") as outfile:
        outfile.write(compressed_input)


# flatten status to dataframe
def _compare_status(status):
    converted_status = []
    for key1 in status:
        if type(status[key1]) == dict or type(status[key1]) == list:
            for key2 in status[key1]:
                full_key = key1+"_"+key2
                
                if type(status[key1]) == dict and (type(status[key1][key2]) == dict or type(status[key1][key2]) == list):
                    for key3 in status[key1][key2]:
                        full_key = key1+"_"+key2+"_"+key3

                        if type(status[key1][key2][key3]) == dict or type(status[key1][key2][key3]) == list:
                            for key4 in status[key1][key2][key3]:
                                full_key = key1+"_"+key2+"_"+key3

                                if type(status[key1][key2][key3]) == dict and (type(status[key1][key2][key3][key4]) == dict or type(status[key1][key2][key3][key4]) == list):
                                    for key5 in status[key1][key2][key3][key4]:
                                        full_key = key1+"_"+key2+"_"+key3+"_"+key4
                                        if type(status[key1][key2][key3][key4]) == dict: converted_status.append(full_key+status[key1][key2][key3][key4][key5])
                                        if type(status[key1][key2][key3][key4]) == list and type(key5) != list: converted_status.append(full_key+"_"+key5)

                                else:
                                    if type(status[key1][key2][key3]) == dict: converted_status.append(full_key+status[key1][key2][key3][key4])
                                    if type(status[key1][key2][key3]) == list: converted_status.append(full_key+"_"+key4)

                        else:
                            if type(status[key1][key2][key3]) == dict: converted_status.append(full_key+"_"+status[key1][key2][key3])
                            if type(status[key1][key2][key3]) == list: converted_status.append(full_key+"_"+key3)

                else:
                    if type(status[key1]) == dict: converted_status.append(full_key+"_"+status[key1][key2])
                    if type(status[key1]) == list: converted_status.append(full_key+"_"+key2)

        else:
            if type(status[key1]) == dict: converted_status.append(key1+"_"+status[key1])

    return converted_status


def compare_status(status1, status2):
    converted_status1 = _compare_status(status1)
    converted_status2 = _compare_status(status2)

    mismatches1 = [entry for entry in converted_status1 if entry not in converted_status2]
    mismatches2 = [entry for entry in converted_status2 if entry not in converted_status1]

    print("<", len(mismatches1), "/", len(mismatches2), "were detected")
    return {"status1": sorted(mismatches1), "status2": sorted(mismatches2)}


def convert_cancer_scores(cancer_scores, feature_targets, id_targets, label_targets, datatype="raw", exclude_misses=False, randomize=False, selected_samples=[]):
    data = {**{"ID:project": []}, **{"ID:"+id.replace("FEATURE:", "").replace("ID:", ""): [] for id in id_targets if id not in label_targets},
            **{"FEATURE:"+feature.replace("FEATURE:", "").replace("ID:", ""): [] for feature in feature_targets if feature not in label_targets},
            **{"LABEL:"+label.replace("FEATURE:", "").replace("ID:", ""): [] for label in label_targets}}
    
    excluded_index = []
    for i in range(cancer_scores.shape[0]):
        for feature in feature_targets:
            if feature not in label_targets:
                values = cancer_scores.iloc[i].loc[feature+"_"+datatype]

                for j in range(len(values)):
                    data["FEATURE:"+feature.replace("FEATURE:", "").replace("ID:", "")].append(values[j])
                    
                    if pd.isna(values[j]) == True and len(data["FEATURE:"+feature.replace("FEATURE:", "").replace("ID:", "")])-1 not in excluded_index:
                        excluded_index.append(len(data["FEATURE:"+feature.replace("FEATURE:", "").replace("ID:", "")])-1)

        for id in id_targets:
            values = json.loads(cancer_scores.iloc[i].loc[id].replace("'", "\""))

            for value in values:
                data["ID:"+id.replace("FEATURE:", "").replace("ID:", "")].append(value)

        for value in values:
            data["ID:project"].append(cancer_scores.iloc[i].loc["project"])

        for label in label_targets:
            values = cancer_scores.iloc[i].loc[label+"_"+datatype]

            for j in range(len(values)):
                data["LABEL:"+label.replace("FEATURE:", "").replace("ID:", "")].append(values[j])

                if pd.isna(values[j]) == True and len(data["LABEL:"+label.replace("FEATURE:", "").replace("ID:", "")])-1 not in excluded_index:
                    excluded_index.append(len(data["LABEL:"+label.replace("FEATURE:", "").replace("ID:", "")])-1)

    data = pd.DataFrame(data)
    if exclude_misses == True:    data = data[[False if i in excluded_index else True for i in range(data.shape[0])]]
    if len(selected_samples) > 0: data = data[data["ID:sample_names"].isin(selected_samples)]
    if randomize == True:         data = get_random_cohorts(data, 5, targets=[id_targets[0]])
    return data


# load all lists stored as str to lists
def convert_raw_data(cancer_scores):
    test = [None for i in range(cancer_scores.shape[0])]
    for col in cancer_scores.columns:
        if "raw" in col:
            for i in range(cancer_scores.shape[0]):
                values                                        = json.loads(cancer_scores.iloc[i].loc[col].replace("'", "\"").replace("nan", "-1").replace("None", "-1"))
                cancer_scores.at[cancer_scores.index[i], col] = [value if value != -1 else None for value in values]

                if test[i] != None and test[i] != len(cancer_scores.loc[cancer_scores.index[i], col]):
                    print("< size error occurred @load for feature", col, "with", test[i], "/", len(cancer_scores.loc[cancer_scores.index[i]].loc[col]))
                    exit()

    return cancer_scores


def count_clusters(cluster_labels):
    unique, counts = np.unique(cluster_labels, return_counts=True)
    return dict(zip(unique, counts))


def create_case_map(status):
    case_map = {"case_id": [], "submitter_id": []}

    for project in status["case_ids"]:
        if len(status["case_ids"][project]["case_id"]) != len(status["case_ids"][project]["submitter_id"]):
            print("< case and submitter id sizes inconsistent for", project, len(status["case_ids"][project]["case_id"]), "/", len(status["case_ids"][project]["submitter_id"]))

        for i in range(len(status["case_ids"][project]["case_id"])):
            case_map["case_id"].append(status["case_ids"][project]["case_id"][i])
            case_map["submitter_id"].append(status["case_ids"][project]["submitter_id"][i])

    return pd.DataFrame(case_map)


def get_cds_size(hg, transcript_id):
    cds_size   = 0
    exonsstart = [int(i) for i in hg.loc[transcript_id].loc["exonsstart"].split(",") if len(i) > 0] # condition required because last letter is a comma
    exonsend   = [int(i) for i in hg.loc[transcript_id].loc["exonsend"].split(",") if len(i) > 0] # condition required because last letter is a comma

    for i in range(len(exonsend)):
        if ((exonsstart[i] >= hg.loc[transcript_id].loc["cdsstart"] and exonsstart[i] < hg.loc[transcript_id].loc["cdsend"])
            or (exonsend[i] >= hg.loc[transcript_id].loc["cdsstart"] and exonsend[i] < hg.loc[transcript_id].loc["cdsend"])
            or (hg.loc[transcript_id].loc["cdsstart"] >= exonsstart[i] and hg.loc[transcript_id].loc["cdsstart"] <= exonsend[i])
            or (hg.loc[transcript_id].loc["cdsend"] >= exonsstart[i] and hg.loc[transcript_id].loc["cdsend"] <= exonsend[i])):
            if hg.loc[transcript_id].loc["cdsstart"] <= exonsstart[i]: exon_cdsstart = exonsstart[i]
            if hg.loc[transcript_id].loc["cdsstart"] > exonsstart[i]:  exon_cdsstart = hg.loc[transcript_id].loc["cdsstart"]
            if hg.loc[transcript_id].loc["cdsend"] <= exonsend[i]:     exon_cdsend   = hg.loc[transcript_id].loc["cdsend"]
            if hg.loc[transcript_id].loc["cdsend"] > exonsend[i]:      exon_cdsend   = exonsend[i]
            cds_size += exon_cdsend-exon_cdsstart

    return cds_size


def get_values(cancer_scores, feature_filter, feature1, feature2, it, ptc_target,
               class_samples=None, data_type="raw", extract_ptc=False, get_sample_names=False, remove_misses=True):
    values1      = cancer_scores.iloc[it].loc[feature1+"_"+data_type]
    values2      = cancer_scores.iloc[it].loc[feature2+"_"+data_type]
    
    if extract_ptc == True:
        values3      = cancer_scores.iloc[it].loc[ptc_target+"_"+data_type]
    
    if get_sample_names == True:
        sample_names = json.loads(cancer_scores.iloc[it].loc["sample_names"].replace("'", "\""))

    # dimension check
    if len(values1) != len(values2):
        print("size error1 occurred @get_values:", len(values1), "/", len(values2))

    # apply filters
    excluded_index = []
    if len(feature_filter) > 0:
        filter_values = []
        [filter_values.append(cancer_scores.iloc[it].loc[col+"_"+data_type]) for col in feature_filter]

        for i, col in enumerate(feature_filter):
            if remove_misses == False:
                excluded_index.extend([j for j in range(len(filter_values[i])) if pd.isna(filter_values[i][j]) == False
                                       and (filter_values[i][j] >= feature_filter[col][0] or filter_values[i][j] < feature_filter[col][1])])

            if remove_misses == True:
                excluded_index.extend([j for j in range(len(filter_values[i])) if pd.isna(filter_values[i][j]) == True
                                       or filter_values[i][j] >= feature_filter[col][0] or filter_values[i][j] < feature_filter[col][1]])
            
        values1 = [values1[i] for i in range(len(values1)) if i not in excluded_index]
        values2 = [values2[i] for i in range(len(values2)) if i not in excluded_index]
        values3 = [values3[i] for i in range(len(values3)) if i not in excluded_index]

    if get_sample_names == True:
        sample_names = [sample_names[i] for i in range(len(sample_names)) if i not in excluded_index]

    # apply class filter
    if class_samples != None:
        sample_names = json.loads(cancer_scores.iloc[it].loc["sample_names"].replace("'", "\""))
        sample_names = [sample_names[i] for i in range(len(sample_names)) if i not in excluded_index]

        values1      = [values1[i] for i in range(len(values1)) if sample_names[i] in class_samples]
        values2      = [values2[i] for i in range(len(values2)) if sample_names[i] in class_samples]
        values3      = [values3[i] for i in range(len(values3)) if sample_names[i] in class_samples]
        sample_names = [sample_names[i] for i in range(len(sample_names)) if sample_names[i] in class_samples]

    if get_sample_names == False: return values1, values2, values3
    if get_sample_names == True:  return values1, values2, values3, sample_names


def map_data(data, features, target_cols, text="", tag="FEATURE"):
    if len(text) > 0: text = set_bar(text)
    data[1]["index"] = [i for i in range(data[1].shape[0])]
    mapped_index     = append_df_with_mapping(data, target_cols[0], target_cols[1], "index", text)

    mapped_features = {feature: [] for feature in features}
    for i in mapped_index:
        for feature in features:
            if i != "-":
                if data[1].iloc[int(i)].loc[feature] != "#NV": mapped_features[feature].append(data[1].iloc[int(i)].loc[feature])
                else:                                          mapped_features[feature].append(None)

            else:
                mapped_features[feature].append(None)


    for feature in features:
        data[0][tag+":"+feature] = mapped_features[feature]

    return data[0]