import gzip
import json
import os
import math
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import requests
import scipy
from sklearn import metrics
from sklearn.preprocessing import QuantileTransformer, StandardScaler
import time

import statsmodels.api as sm
import sys
from shared_utils import *

# must be imported here and not in shared_utils!
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\NMD_activity")
from extract_mutations_utils import *
sys.path.insert(0, parent_dir+"\\random_forest")
from analyze_predictions_utils import *

cohorts              = 5
data_dir             = parent_dir+r"\data"
fnames               = ["model_training_variants_full.txt"]
loading_modes        = ["pandas"]
## inventory ##
# "analyze_prediction_errors"
# "append_appris_annotation" "append_cancer_types" "append_exons" "append_expressions" "append_genotypes" "append_indel_status" "append_lindeboom_predictions" "append_msk_info" 
# "append_predictions" "append_ptc_info"
# "apply_cuomo_liftover" "apply_score_limits"
# "average_labels"
# "calculate_bulk" "change_defaults" "change_features" "change_placeholders"
# "check_cnvs" "check_cuomo_liftover" "compare_data" "compare_mskcc_chord" "compare_selection"
# "create_lindeboom" "create_shared_variants"
# "filter_by_cds_size" "filter_shared_genes" "filter_shared_variants"
# "get_block_correlation"
# "merge_mskcc_chord"
# "prepare_cuomo_liftover"
# "randomize_cohorts" "rearrange_mutations" "replace_mutation_stats"
# "select_appris" "split_by_driver_genes"
mode                 = "get_block_correlation"
os_sep               = "\\"
separators           = [","]

if len(loading_modes) != len(fnames): print("< warning. loading modes and file tags are not the same number.")


def main():
    data = []

    for i, fname in enumerate(fnames):
        if os.path.isfile(data_dir+os_sep+fname):
            if loading_modes[i] == "pandas":
                with open(data_dir+os_sep+fname, 'r') as f:
                    data.append(pd.read_csv(data_dir+os_sep+fname, delimiter=separators[i], low_memory=False))

            elif loading_modes[i] == "gzip":
                with gzip.open(data_dir+os_sep+fname, 'r') as f:
                    lines = f.readlines()

                for line in lines: print(line)
                data.append(lines)

            elif loading_modes[i] == "lines":
                with open(data_dir+os_sep+fname, 'r') as f:
                    lines = f.readlines()
                    data.append(lines)

                for i, line in enumerate(lines):
                    print(i, line)

            elif loading_modes[i] == "json":
                with open(data_dir+os_sep+fname, 'r') as openfile:
                    data.append(json.load(openfile))

        else: print("<", fname, "not found")


    # used file for appris annotation:
    # hg19_appris_principal_isoforms.txt, downloaded from https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/, most recent version: April 18th 2022
    # hg38_appris_principal_isoforms.txt, downloaded from https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/, most recent version: July 3rd 2024
    if mode == "append_appris_annotation":
        data[0]["mod. transcript id"] = [data[0].iloc[i].loc["transcript id"].split(".")[0] for i in range(data[0].shape[0])]
        data[1]["mod. transcript id"] = [data[1].iloc[i].loc["Transcript ID"].split(".")[0] for i in range(data[1].shape[0])]
        appris_annotations            = append_df_with_mapping(data, "mod. transcript id", "mod. transcript id", "APPRIS Annotation", set_bar("appending appris annotation"), non_redundant=True, reverse=True, verbose=True)

        data[0]["appris annotation"]  = appris_annotations
        data[0]                       = data[0].drop(["mod. transcript id"], axis=1)
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_appended.txt", sep=",", index=False)


    # appending scores by non coding exons (not contained in create_genome_predictions, required for correct assignment of exon position in analyze_predictions)
    if mode == "append_exons":
        stats = {"changes": 0, "no_changes": 0, "ambiguous": [], "dimension_error": [], "exon_size_error": [], "no_predictions": [], "not_found": []}
        bar   = IncrementalBar("appending exons", max=data[0].shape[0])

        for i in range(data[0].shape[0]):
            if data[0].iloc[i].loc["predictions_by_exon"] != "{}":
                data[0].at[data[0].index[i], "predictions_by_exon"] = json.loads(data[0].iloc[i].loc["predictions_by_exon"].replace("{", "{\"").replace("], ", "], \"").replace(":", "\":"))
                selected_data = data[1][data[1]["transcript id"] == data[0].iloc[i].loc["transcript id"]]
                exonsstart = [int(j) for j in selected_data.iloc[0].loc["exonsstart"].split(",") if len(j) > 0] # condition required because last letter is a comma
                exonsend   = [int(j) for j in selected_data.iloc[0].loc["exonsend"].split(",") if len(j) > 0]   # condition required because last letter is a comma
                #print("i", i, data[0].iloc[i].loc["transcript id"], data[0].iloc[i].loc["strand"], selected_data.shape[0], list(data[0].iloc[i].loc["predictions_by_exon"].keys()))

                if selected_data.shape[0] == 1:
                    if len(exonsstart) == len(data[0].iloc[i].loc["predictions_by_exon"]):
                        stats["no_changes"] += 1

                    elif len(exonsstart) < len(data[0].iloc[i].loc["predictions_by_exon"]):
                        stats["dimension_error"].append(data[0].iloc[i].loc["transcript id"])

                    elif len(exonsstart) > len(data[0].iloc[i].loc["predictions_by_exon"]):
                        exon_keys = [int(key) for key in list(data[0].iloc[i].loc["predictions_by_exon"].keys())]
                        predictions_by_exon = {}
                        if selected_data.iloc[0].loc["strand"] == "+":
                            exons = []
                            for j in range(len(exonsstart)):
                                exon = exonsstart[j]
                                # convert exon to cdsstart according to create_genome_predictions if it deviates
                                if selected_data.iloc[0].loc["cdsstart"] >= exonsstart[j] and selected_data.iloc[0].loc["cdsstart"] < exonsend[j]:
                                    exon = selected_data.iloc[0].loc["cdsstart"]

                                exons.append(exon)

                            exons = np.unique([*exons, *exon_keys]).tolist()
                            exons = sorted(exons)

                        if selected_data.iloc[0].loc["strand"] == "-":
                            exons = []
                            for j in range(len(exonsstart)):
                                exon = exonsstart[j]
                                # convert exon to cdsstart according to create_genome_predictions if it deviates
                                if selected_data.iloc[0].loc["cdsstart"] >= exonsstart[j] and selected_data.iloc[0].loc["cdsstart"] < exonsend[j]:
                                    exon = selected_data.iloc[0].loc["cdsstart"]+3

                                exons.append(exon)
                            
                            exons = np.unique([*exons, *exon_keys]).tolist()
                            exons = sorted(exons, reverse=True)

                        if len(exons) == len(exonsstart):
                            for exon in exons:
                                if str(exon) in data[0].iloc[i].loc["predictions_by_exon"]:
                                    predictions_by_exon[str(exon)] = data[0].iloc[i].loc["predictions_by_exon"][str(exon)]
                                
                                else:
                                    predictions_by_exon[str(exon)] = []

                            data[0].at[data[0].index[i], "predictions_by_exon"] = predictions_by_exon
                            stats["changes"] += 1

                            if i >= 20 and i < 40:
                                print(data[0].iloc[i].loc["transcript id"], data[0].iloc[i].loc["strand"])
                                for key in predictions_by_exon:
                                    print(key, len(predictions_by_exon[key]))


                        else:
                            stats["exon_size_error"].append(data[0].iloc[i].loc["transcript id"])

                elif selected_data.shape[0] == 0:
                    stats["not_found"].append(data[0].iloc[i].loc["transcript id"])

                elif selected_data.shape[0] > 1:
                    stats["ambiguous"].append(data[0].iloc[i].loc["transcript id"])

            else:
                stats["no_predictions"].append(data[0].iloc[i].loc["transcript id"])

            bar.next()
        bar.finish()

        print(json.dumps(stats, indent=4))
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_appended.txt", sep=",", index=False)


    # appeding expression data for TCGA-projects
    if mode == "append_expressions":
        expression_tag    = "mean_tpm_unstranded"
        expression_target = "gene_id" # TCGA: "gene_id"
        variant_target    = "ID:gene id" # TCGA: "ID:gene id

        if expression_target == "gene_id": data[1][expression_target] = [data[1].iloc[i].loc[expression_target].split(".")[0] for i in range(data[1].shape[0])]
        print("< gene duplicates", data[1].shape[0]-data[1].drop_duplicates(subset=[expression_target]).shape[0])
        if "FEATURE:expression" not in data[0]: data[0].insert(data[0].shape[1], "FEATURE:expression", [None for _ in range(data[0].shape[0])])
        else:                                   data[0]["FEATURE:expression"] = [None for _ in range(data[0].shape[0])]

        misses = []
        bar    = IncrementalBar(set_bar("appending expressions"), max=data[0].shape[0])

        for i in range(data[0].shape[0]):
            gene_id          = data[0].iloc[i].loc[variant_target]
            project          = data[0].iloc[i].loc["ID:project"]
            selected_gene_id = data[1][data[1][expression_target] == gene_id]
            expressions      = [value for value in selected_gene_id[project+"_"+expression_tag] if value > 0]

            if len(expressions) > 0:           data[0].at[data[0].index[i], "FEATURE:expression"] = np.mean(expressions)
            if selected_gene_id.shape[0] == 0: misses.append(gene_id)

            if len(expressions) > 1:
                print("i", i, gene_id, project+"_"+expression_tag, np.mean(expressions))
                print(expressions)

            if math.isinf(np.mean(expressions)) == True or pd.isna(np.mean(expressions)) == True:
                print("i", i, gene_id, project+"_"+expression_tag, np.mean(expressions))
                print(expressions)

            bar.next()
        bar.finish()

        print("<", len(misses), "missing gene ids")
        print(misses)
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_expressions.txt", sep=",", index=False)


    # maps genotypes of HPSI cells and allows filtering
    # file containing genotype info has lower size because for one cell line no data were available
    if mode == "append_genotypes":
        print(data[0])
        print(data[1])
        data[0]["combined id"] = [data[0].iloc[i].loc["ID:cell id"]+data[0].iloc[i].loc["ID:variant id"] for i in range(data[0].shape[0])]
        data[1]["combined id"] = [data[1].iloc[i].loc["ID:cell id"]+data[1].iloc[i].loc["ID:variant id"] for i in range(data[1].shape[0])]
        genotypes              = append_df_with_mapping(data, "combined id", "combined id", "ID:heterozygous", "appending genotypes",
                                                        non_redundant=True, reverse=True, verbose=True)

        data[0]["ID:heterozygous"] = genotypes
        data[0]                    = data[0][data[0]["ID:heterozygous"] == "True"]
        data[0]                    = data[0].drop(columns=["combined id"])
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].strip(".txt")+"_gt_filtered.txt", sep=",", index=False)


    if mode == "append_msk_info":
        ptc_mutations = append_df_with_mapping(data, "ID:patient id", "ID:PATIENT_ID", "FEATURE:ptc_mutations", set_bar("appending msk info"),
                                               non_redundant=True, reverse=True, verbose=True)

        ptc_mutations = [float(ptc_mutation) if ptc_mutation != "-" else None for ptc_mutation in ptc_mutations]

        if "FEATURE:ptc_mutations" not in data[0].columns.tolist(): data[0].insert(data[0].shape[1], 'FEATURE:ptc_mutations', ptc_mutations)
        else:                                                       data[0]["FEATURE:ptc_mutations"] = ptc_mutations

        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_appended.txt", sep=",", index=False)


    # appending gene ids from gProfiler_hsapiens_11.12.2023_18-37-10.csv
    if mode == "append_gene_ids_knownGene":
        data[0]["mod. transcript id"] = [data[0].iloc[i].loc["transcript id"].split(".")[0] for i in range(data[0].shape[0])]
        gene_ids                      = append_df_with_mapping(data, "mod. transcript id", "initial_alias", "converted_alias", "appending gene ids", non_redundant=False, reverse=True)
        data[0]                       = data[0].drop(["mod. transcript id"], axis=1)
        ambiguity_test               = [gene_id for gene_id in gene_ids if len(gene_id) > 1]
        print("<", len(ambiguity_test), "ambiguities found.")
        data[0].insert(1, "gene id", [gene_id[0] if len(gene_id) > 0 else "-" for gene_id in gene_ids]) # flatten gene_ids
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].strip(".txt")+"_appended.txt", sep=",", index=False)


    # appending gene ids from gene names.txt
    if mode == "append_gene_symbols_knownGene":
        print(data[0])
        data[0]["mod. transcript id"] = [data[0].iloc[i].loc["transcript id"].split(".")[0] for i in range(data[0].shape[0])]
        data[1]["mod. transcript id"] = [data[1].iloc[i].loc["MANE Select Ensembl transcript ID (supplied by NCBI)"].split(".")[0]
                                         if pd.isna(data[1].iloc[i].loc["MANE Select Ensembl transcript ID (supplied by NCBI)"]) == False
                                         else None
                                         for i in range(data[1].shape[0])]
        gene_symbols                  = append_df_with_mapping(data, "mod. transcript id", "mod. transcript id", "Approved symbol", set_bar("appending gene symbols"),
                                                               non_redundant=False, reverse=True)
        data[0]                       = data[0].drop(["mod. transcript id"], axis=1)
        ambiguity_test               = [gene_symbol for gene_symbol in gene_symbols if len(gene_symbol) > 1]
        print("<", len(ambiguity_test), "ambiguities found.")
        data[0].insert(1, "gene symbol", [gene_symbol[0] if len(gene_symbol) > 0 else "-" for gene_symbol in gene_symbols]) # flatten gene_ids
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_appended.txt", sep=",", index=False)


    if mode == "append_indel_status":
        datatype = "tcga" # "tcga" "msk"

        if datatype == "msk":
            hgvscs = data[0]["ID:HGVSc"].tolist()

        show = True
        report  = {"nonsense": 0, "del-1": 0, "ins+1": 0, "del<-1": 0, "ins>+1": 0, "delins": 0}
        multiple_indels = []
        for i in range(data[0].shape[0]):
            if datatype == "msk":
                data[0].at[data[0].index[i], "ID:HGVSc"] = data[0].iloc[i].loc["ID:HGVSc"].split(":")[1]
                
            # second condition due to "c.1909-6_1911dup"
            if "delins" in data[0].iloc[i].loc["ID:HGVSc"] in data[0].iloc[i].loc["ID:HGVSc"] or "-" in data[0].iloc[i].loc["ID:HGVSc"]:
                multiple_indels.append(True)
                report["delins"] += 1

            elif "del" in data[0].iloc[i].loc["ID:HGVSc"]: #"ins" in data[0].iloc[i].loc["ID:HGVSp"] or "del" in data[0].iloc[i].loc["ID:HGVSp"]:
                counts = get_numbers(data[0].iloc[i].loc["ID:HGVSc"])

                # differentiating datatypes required as single deletions in msk are formatted as c.1991delA and c.6594_6595delAC or c.881_890del (!)
                if datatype == "msk":
                    if len(counts) > 1: indels = counts[1]-counts[0]+1
                    else:               indels = 1
                    if "ins" in data[0].iloc[i].loc["ID:HGVSc"]: indels -= len(data[0].iloc[i].loc["ID:HGVSc"].split("ins")[1])

                # differentiating datatypes required as single deletions in tcga are formatted as e.g. c.1454del or c.342_343del
                if datatype == "tcga":
                    if len(counts) > 1: indels = counts[1]-counts[0]+1
                    else:               indels = 1
                    if "ins" in data[0].iloc[i].loc["ID:HGVSc"]: indels -= len(data[0].iloc[i].loc["ID:HGVSc"].split("ins")[1])

                if indels == 1:
                    multiple_indels.append(False)
                    report["del-1"] += 1
                    if i < 100: print("i1", i, data[0].iloc[i].loc["ID:HGVSc"])

                if indels > 1:
                    multiple_indels.append(True)
                    report["del<-1"] += 1
                    if i < 100: print("i2", i, data[0].iloc[i].loc["ID:HGVSc"])

            elif "dup" in data[0].iloc[i].loc["ID:HGVSc"]: #"ins" in data[0].iloc[i].loc["ID:HGVSp"] or "del" in data[0].iloc[i].loc["ID:HGVSp"]:
                counts = get_numbers(data[0].iloc[i].loc["ID:HGVSc"])
                if len(counts) > 1: indels = counts[1]-counts[0]+1
                else:               indels = 1

                if indels == 1:
                    multiple_indels.append(False)
                    report["ins+1"] += 1
                    if i < 100: print("i3", i, data[0].iloc[i].loc["ID:HGVSc"])

                if indels > 1:
                    multiple_indels.append(True)
                    report["ins>+1"] += 1
                    if i < 100: print("i4", i, data[0].iloc[i].loc["ID:HGVSc"])

            elif "ins" in data[0].iloc[i].loc["ID:HGVSc"]: #"ins" in data[0].iloc[i].loc["ID:HGVSp"] or "del" in data[0].iloc[i].loc["ID:HGVSp"]:
                counts = get_numbers(data[0].iloc[i].loc["ID:HGVSc"])
                indels = len(data[0].iloc[i].loc["ID:HGVSc"].split("ins")[1])

                if indels == 1:
                    multiple_indels.append(False)
                    report["ins+1"] += 1
                    if i < 100: print("i5", i, data[0].iloc[i].loc["ID:HGVSc"])

                if indels > 1:
                    multiple_indels.append(True)
                    report["ins>+1"] += 1
                    if i < 100: print("i6", i, data[0].iloc[i].loc["ID:HGVSc"])

            else:
               report["nonsense"] += 1
               multiple_indels.append(False)

            if len(multiple_indels) != i+1:
                print(i, data[0].iloc[i].loc["ID:HGVSc"], len(multiple_indels))
                multiple_indels.append(True)

        print(report)
        if datatype == "msk":
            data[0]["ID:HGVSc"] = hgvscs

        data[0]["ID:multiple indels"] = multiple_indels
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_appended.txt", sep=",", index=False)


    # appending of isoform info from appris_principal_isoforms.txt, applied to UCSC knownGene
    # knownGene hg19 downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/, version from: June 30th 2013
    # knownGene hg38 downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/, version from: April 2nd 2024
    if mode == "append_isoforms_knownGene":
        data[0]["mod. transcript id"] = [data[0].iloc[i].loc["transcript id"].split(".")[0] for i in range(data[0].shape[0])]
        appris_annotations            = append_df_with_mapping(data, "mod. transcript id", "Transcript ID", "APPRIS Annotation", "appending APPRIS annotation", non_redundant=True, reverse=True, verbose=False)
        data[0]                       = data[0].drop(["mod. transcript id"], axis=1)
        data[0]["appris annotation"]  = appris_annotations
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].strip(".txt")+"_appended.txt", sep=",", index=False)


    if mode == "append_lindeboom_predictions":
        # extended variant id required for Cuomo data
        if "ID:variant transcript pos id" in data[1].columns:
            print("< Cuomo data detected.")
            data[0]["ID:variant transcript pos id"] = [data[0].iloc[i].loc["ID:variant id"] + "_" + data[0].iloc[i].loc["ID:transcript id"] + "_" + str(data[0].iloc[i].loc["FEATURE:ptc cds position"])
                                                       for i in range(data[0].shape[0])]
            lindeboom_predictions = append_df_with_mapping(data, "ID:variant transcript pos id", "ID:variant transcript pos id", "FEATURE:lindeboom prediction",
                                                           "appending lindeboom predictions", non_redundant=True, reverse=True, show_progress=True)
            #lindeboom_predictions = append_df_with_mapping2(data, "ID:variant transcript pos id", "ID:variant transcript pos id", "FEATURE:lindeboom prediction", text="progress", verbose=False)
            #

        else:
            lindeboom_predictions = append_df_with_mapping(data, "ID:variant id", "ID:variant id", "FEATURE:lindeboom prediction",
                                                           "appending lindeboom predictions", non_redundant=True, reverse=True, show_progress=True)

        lindeboom_predictions = [pred if pred != "-" else None for pred in lindeboom_predictions]
        if 'FEATURE:lindeboom prediction' not in data[0].columns.tolist(): data[0].insert(data[0].shape[1], 'FEATURE:lindeboom prediction', lindeboom_predictions)
        else:                                                              data[0]["FEATURE:lindeboom prediction"] = lindeboom_predictions

        test = False

        if test == True:
            last_variant_id = ""
            for i in range(data[1].shape[0]):
                if data[1].iloc[i].loc["ID:variant id"] != last_variant_id:
                    selected_df = data[0][data[0]["ID:variant id"] == data[1].iloc[i].loc["ID:variant id"]]
                    if i % 1000 == 0: print("i", i, "/", data[1].shape[0], ":", selected_df.shape[0])

                    tests = [selected_df.iloc[j].loc["FEATURE:lindeboom prediction"] for j in range(selected_df.shape[0])
                            if str(data[1].iloc[i].loc["FEATURE:lindeboom prediction"]) != str(selected_df.iloc[j].loc["FEATURE:lindeboom prediction"])]
                    if len(tests) > 0:
                        print("i2", i, data[1].iloc[i].loc["ID:variant id"], tests, "/", data[1].iloc[i].loc["FEATURE:lindeboom prediction"])

                        last_variant_id = data[1].iloc[i].loc["ID:variant id"]

        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_lindeboom.txt", sep=",", index=False)


    # not to use for Cuomo data due to ambiguous variant ids (instead, produce full output with read_predictions)
    if mode == "append_predictions":
        target      = "ID:variant id" # "ID:variant id"
        predictions = append_df_with_mapping(data, target, target, "FEATURE:prediction", set_bar("appending predictions"), non_redundant=True, reverse=True, verbose=True)
        predictions = [float(prediction) if prediction != "-" else None for prediction in predictions]

        if "FEATURE:prediction" not in data[0].columns.tolist(): data[0].insert(data[0].shape[1], 'FEATURE:prediction', predictions)
        else:                                                    data[0]["FEATURE:prediction"] = predictions

        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_preds.txt", sep=",", index=False)


    # tool to pre-calculate per-sample information (e.g. ptc mutations), used in analyze_predictions
    # gender information is required (mapped using tcga_tools)
    if mode == "append_ptc_info":
        targets = ["FEATURE:ptc_mutations", "FEATURE:frameshift"]

        if "FEATURE:frameshift" in targets:
            if "FEATURE:frameshift" not in data[0].columns:           data[0].insert(data[0].shape[1], "FEATURE:frameshift", [None for _ in range(data[0].shape[0])])
            else:                                                     data[0]["FEATURE:frameshift"] = [None for _ in range(data[0].shape[0])]
            if "FEATURE:frameshift_mutations" not in data[0].columns: data[0].insert(data[0].shape[1], "FEATURE:frameshift_mutations", [None for _ in range(data[0].shape[0])])
            else:                                                     data[0]["FEATURE:frameshift_mutations"] = [None for _ in range(data[0].shape[0])]
        
        if "FEATURE:ptc_mutations" in targets:
            if "FEATURE:ptc_mutations" not in data[0].columns:        data[0].insert(data[0].shape[1], "FEATURE:ptc_mutations", [None for _ in range(data[0].shape[0])])
            else:                                                     data[0]["FEATURE:ptc_mutations"] = [None for _ in range(data[0].shape[0])]

        bar = IncrementalBar("appending ptc info", max=data[0].shape[0])
        for i in range(data[0].shape[0]):
            if "FEATURE:frameshift" in targets:
                data[0].at[data[0].index[i], "FEATURE:frameshift"] = get_frameshift(data[0].iloc[i].loc["ID:HGVSp"])
                
                if data[0].loc[data[0].index[i]].loc["FEATURE:frameshift"] > 0: data[0].at[data[0].index[i], "FEATURE:frameshift_mutations"] = 1
                else:                                                           data[0].at[data[0].index[i], "FEATURE:frameshift_mutations"] = 0
           
            if "FEATURE:ptc_mutations" in targets:
                # apply gender and chromosome-based selection prior to calcuation
                selected_data = data[0][data[0]["ID:case id"] == data[0].iloc[i].loc["ID:case id"]]
                selected_data = selected_data[selected_data["ID:chromosome"] != "chrY"]
                selected_data = selected_data[[True if selected_data.iloc[j].loc["ID:chromosome"] != "chrX"
                                               or (selected_data.iloc[j].loc["ID:chromosome"] == "chrX" and selected_data.iloc[j].loc["LABEL:gender"] == "FEMALE") else False
                                               for j in range(selected_data.shape[0])]]

                data[0].at[data[0].index[i], "FEATURE:ptc_mutations"] = selected_data.shape[0]
                if i < 10 or "delins" in data[0].iloc[i].loc["ID:HGVSp"]:
                    print("i", i, data[0].iloc[i].loc["ID:HGVSp"], data[0].loc[data[0].index[i]].loc["FEATURE:frameshift"],
                        data[0].loc[data[0].index[i]].loc["FEATURE:frameshift_mutations"], data[0].loc[data[0].index[i]].loc["FEATURE:ptc_mutations"])

            bar.next()
        bar.finish()
        data[0] = data[0].rename(columns={col: col if col != "LABEL:gender" else "ID:gender" for col in data[0].columns})
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_ptc_info.txt", sep=",", index=False)


    # as in most cases, transcript id was not sufficient to find uniprot id, gene ids were retrieved in a second step
    if mode == "append_uniprot_ids":
        data[0]["transcript id_"] = [data[0].iloc[i].loc["transcript id"].split(".")[0] for i in range(data[0].shape[0])]
        data[1]["transcript id_"] = [data[1].iloc[i].loc["MANE Select Ensembl transcript ID (supplied by NCBI)"].split(".")[0]
                                     if type(data[1].iloc[i].loc["MANE Select Ensembl transcript ID (supplied by NCBI)"]) == str else "-" for i in range(data[1].shape[0])]

        uniprot_ids1 = append_df_with_mapping(data, "transcript id_", "transcript id_", "UniProt ID(supplied by UniProt)", "appending uniprot ids", non_redundant=True, reverse=True, verbose=False)
        data[0] = data[0].drop(columns="transcript id_")

        uniprot_ids2 = append_df_with_mapping(data, "gene id", "Ensembl gene ID", "UniProt ID(supplied by UniProt)", "appending uniprot ids", non_redundant=True, reverse=True, verbose=True)
        if "FEATURE:prediction" not in data[0].columns.tolist(): data[0].insert(data[0].shape[1], 'uniprot id', [uniprot_ids2[i] if uniprot_ids1[i] == "-" else uniprot_ids1[i] for i in range(len(uniprot_ids1))])
        else:                                                    data[0]["uniprot id"] = [uniprot_ids2[i] if uniprot_ids1[i] == "-" else uniprot_ids1[i] for i in range(len(uniprot_ids1))]

        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].strip(".txt")+"_appended.txt", sep=",", index=False)


    # appending ucids from hg19_knownToEnsembl.txt (not required for hg38)
    # downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/, version from: April 6th 2014
    if mode == "append_transcript_ids_knownGene":
        data[0]["uc id"] = [data[0].iloc[i].loc["uc id"].split(".")[0] for i in range(data[0].shape[0])]
        data[1]["uc id"] = [data[1].iloc[i].loc["uc id"].split(".")[0] for i in range(data[1].shape[0])]

        transcript_ids   = append_df_with_mapping(data, "uc id", "uc id", "transcript id", "appending knownGene entries", non_redundant=True, reverse=True)
        data[0].insert(0, 'transcript id', transcript_ids)
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".txt")[0]+"_appended.txt", sep=",", index=False)


    if mode == "apply_wes_liftover":
        start_hg38 = [0 for _ in range(data[0].shape[0])]
        end_hg38   = [0 for _ in range(data[0].shape[0])]
        for i in range(data[1].shape[0]):
            start_hg38[data[1].iloc[i].loc["Pos"]] = data[1].iloc[i].loc["Start"]
            end_hg38[data[1].iloc[i].loc["Pos"]]   = data[1].iloc[i].loc["End"]

        data[0]["Start_hg38"] = start_hg38
        data[0]["End_hg38"]   = end_hg38
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_hg38.csv", sep="\t", index=False)


    # averages and reduces data for forestNMD
    if mode == "average_labels":
        target_col = "ID:transcript pos id"

        if "ID:transcript pos id" not in data[0].columns:
            data[0]["ID:transcript pos id"] = [data[0].iloc[i].loc["ID:transcript id"] + "_" + str(data[0].iloc[i].loc["FEATURE:ptc cds position"]) for i in range(data[0].shape[0])]

        if "LABEL:NMD score" not in data[0].columns:
            print("< error. label not found.")
            exit()

        data[0]       = data[0].sort_values(by=target_col)
        block_index   = create_blocks(data[0], [target_col])
        averaged_data = {col: [] for col in data[0].columns}

        bar = IncrementalBar(set_bar("averaging labels"), max=len(block_index))
        for i in range(len(block_index)):
            for col in data[0].columns:
                if col == "LABEL:alt counts":     averaged_data[col].append(data[0].iloc[block_index[i]["index"]][col].sum())
                elif col == "LABEL:total counts": averaged_data[col].append(data[0].iloc[block_index[i]["index"]][col].sum())
                elif col == "LABEL:NMD score":    averaged_data[col].append(data[0].iloc[block_index[i]["index"]][col].mean())
                else:                             averaged_data[col].append(data[0].iloc[block_index[i]["index"][0]].loc[col])
            
            bar.next()
        bar.finish()

        averaged_data = pd.DataFrame(averaged_data)
        print("< averaging labels reduced data from", data[0].shape[0], "to", averaged_data.shape[0])
        averaged_data.to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_reduced_avg.txt", sep=",", index=False)


    if mode == "change_defaults":
        test = data[0].copy()
        data[0]["FEATURE:dist. from last EJC"] = [-data[0].iloc[i].loc["FEATURE:ptc upstream distance"] if pd.isna(data[0].iloc[i].loc["FEATURE:dist. from last EJC"]) == True
                                                  else data[0].iloc[i].loc["FEATURE:dist. from last EJC"]
                                                  for i in range(data[0].shape[0])]

        for i in range(data[0].shape[0]):
            if data[0].iloc[i].loc["FEATURE:dist. from last EJC"] != test.iloc[i].loc["FEATURE:dist. from last EJC"]:
                print(i, data[0].iloc[i].loc["FEATURE:ptc upstream distance"], data[0].iloc[i].loc["FEATURE:dist. from last EJC"], test.iloc[i].loc["FEATURE:dist. from last EJC"])

        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_changed.txt", sep=",", index=False)


    if mode == "change_features":
        updated_cols = data[0].columns.tolist()

        print("features:")
        for i in range(len(updated_cols)):
            print(i, " : ", updated_cols[i])
        
        print("< selection by index:on/off or indexA-indexB:on/off, example(1:on). end selection typing 'x'.")
        print("  alternatively, selection-based profile: A-H")
        console_input = ""
        while console_input != "x" and console_input not in ["A", "B", "C", "D", "E", "F", "G", "H"]:
            console_input = input()
            if console_input != "x" and console_input not in ["A", "B", "C", "D", "E", "F", "G", "H"]:
                try:
                    if "-" not in console_input:
                        index = int(console_input.split(":")[0])

                        if console_input.split(":")[1] == "on":
                            updated_cols[index] = "FEATURE:" + updated_cols[index].split(":")[1]

                        if console_input.split(":")[1] == "off": updated_cols[index] = "FEATURxE:" + updated_cols[index].split(":")[1]

                    else:                   
                        indexstr = console_input.split(":")[0]
                        index1   = int(indexstr.split("-")[0])
                        index2   = int(indexstr.split("-")[1])

                        for i in range(index1, index2+1, 1):
                            if console_input.split(":")[1] == "on":  updated_cols[i] = "FEATURE:" + updated_cols[i].split(":")[1]
                            if console_input.split(":")[1] == "off": updated_cols[i] = "FEATURxE:" + updated_cols[i].split(":")[1]
     
                except:
                    print("< error occurred. try again.")

            if console_input in ["A", "B", "C", "D", "E", "F", "G", "H"]:
                if console_input == "A":
                    selected_features = ["3'utr-size", "downstream exons", "downstream EJC density", "ptc downstream distance",
                                         "dist. from last EJC", "ptc exon position", "ptc-wt stop codon distance",
                                         "downstream GC ptc count ptc exon"]

                if console_input == "B":
                    selected_features = ["3'utr-size", "downstream exons", "downstream EJC density", "ptc downstream distance",
                                         "dist. from last EJC", "ptc exon position", "ptc-wt stop codon distance",
                                         "downstream GC ptc count ptc exon", "lindeboom prediction"]

                if console_input == "C":
                    selected_features = ["downstream exons", "last exon", "ptc cds position", "ptc-wt stop codon distance", "ptc exon size",
                                        "ptc downstream distance", "upstream exons", 
                                        "rna half-life", "50 nt to last EJC", "AF", "5'utr-size", "3'utr-size", "total exon size"]

                if console_input == "D":
                    selected_features = ["dist. from last EJC", "last exon", "ptc cds position", "ptc-wt stop codon distance", "ptc exon size",
                                        "rna half-life", "50 nt to last EJC", "AF"]

                if console_input == "E":
                    selected_features = ["3'utr-size", "downstream exons", "downstream EJC density", "ptc downstream distance",
                                         "dist. from last EJC", "ptc exon position", "ptc-wt stop codon distance",
                                         "downstream GC ptc count ptc exon", "5'utr-size", "ptc exon size", "ptc cds position", "rna half-life",
                                         "last exon", "50 nt to last EJC", "AF", "total exon size", "upstream exons"]

                if console_input == "F":
                    selected_features = ["5'utr-size", "3'utr-size", "upstream cds exons", "downstream cds exons", "cds EJC density",
                                         "downstream cds EJC density", "ptc downstream distance", "dist. from last EJC", "ptc cds position",
                                         "ptc exon position", "ptc-wt stop codon distance", "total GC count exons", "downstream GC ptc count ptc exon",
                                         "50 nt to last cds EJC", "last cds exon", "ptc exon size", "mean expression", "median expression",
                                         "AF"]

                if console_input == "G":
                    selected_features = [col.split(":")[1] for col in updated_cols if "prediction" in col]

                if console_input == "H":
                    selected_features = ["5'utr-size", "3'utr-size", "upstream cds exons", "downstream cds exons", "cds EJC density", "downstream cds EJC density", "ptc downstream distance",
                                         "dist. from last EJC", "ptc cds position", "ptc exon position", "ptc-wt stop codon distance", "last cds exon", "50 nt to last cds EJC", "ptc exon size",
                                         "total GC count exons", "downstream GC ptc count ptc exon", "lindeboom prediction"]
                    

                updated_cols = ["FEATURE:"+updated_col.split(":")[1]
                                if "FEATURxE" in updated_col and updated_col.split(":")[1] in selected_features else updated_col
                                for updated_col in updated_cols]
                updated_cols = ["FEATURxE:"+updated_col.split(":")[1]
                                if "FEATURE" in updated_col and updated_col.split(":")[1] not in selected_features else updated_col
                                for updated_col in updated_cols]
        
        data[0].columns = updated_cols
        print(data[0].columns.tolist())
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0], sep=",", index=False)


    # for genome build 19, downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/, version: August 21st 2018
    # for genome build 38, downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p14/, version: October 27th 2022 
    if mode == "convert_genomes":
        hg_build_dir = "hg19"
        genome       = load_whole_genome(data_dir+"\\hg19.fa")

        if not os.path.isdir(data_dir+"\\"+hg_build_dir): os.mkdir(data_dir+"\\"+hg_build_dir)
        for chr_tag in genome:
            with open(data_dir+"\\"+hg_build_dir+"\\"+chr_tag, 'w') as f:
                f.write(">" + chr_tag + "\n")
                f.write(genome[chr_tag])
        

    if mode == "count_isoforms":
        data[0]        = data[0].sort_values(by="gene id")
        # calling create_blocks must be corrected (3rd argument removed)
        block_index    = create_blocks(data[0], "gene id", "gene id")
        isoform_counts = [0 for _ in range(data[0].shape[0])]

        for i in range(len(block_index)):
            if block_index[i]["block id"] != "None" and block_index[i]["block id"] != None and block_index[i]["block id"] != "-":
                for j in block_index[i]["index"]:
                    isoform_counts[j] = len(block_index[i]["index"])

        data[0]["isoform counts"] = isoform_counts
        data[0].to_csv(path_or_buf=data_dir+os_sep+fnames[0].strip(".txt") + "_isoform_counts.txt", sep=",", index=False)


    if mode == "create_genome_selection":
        # genome_build
        genome_build = "hg19" # "hg19" "h38"

        if genome_build == "hg19":
            data[0]["transcript id"] = [i.split(".")[0] for i in data[0]["uc id"]]
            transcript_ids           = [i.split(".")[0] for i in data[1]["ID:uc id"]]

        if genome_build == "hg38":
            data[0]["transcript id"] = [i.split(".")[0] for i in data[0]["transcript id"]]
            transcript_ids           = [i.split(".")[0] for i in data[1]["ID:transcript id"]]

        genome_selection = data[0][data[0]["transcript id"].isin(transcript_ids)]
        genome_selection.to_csv(path_or_buf=data_dir+"\\"+fnames[0].split(".")[0]+"_selection.txt", sep=",", index=False)


    # custom pre-filter for Cuomo data to exclude data with inconsistent cds size (due to usage of different genome builds)
    if mode == "filter_by_cds_size":
        apu = Analyze_predictions_utils({"appris_selection"  : False,
                                         "create_newdir"     : False,
                                         "data_dir"          : r"C:\Programming\Translational_genomics\NMD_analysis\data",
                                         "error_filter"      : ["chromosome_not_found", "total_exon_size_unknown"],
                                         "os_sep"            : "\\",
                                         "prediction_fnames" : [fnames[1]],
                                         "use_variant_filter": False})
        apu.load([fnames[1]], mode="predictions")

        init_shape = data[0].shape[0]
        data[0]    = data[0][[True if data[0].iloc[i].loc["ID:transcript id"] in apu.predictions.index
                              and data[0].iloc[i].loc["FEATURE:total cds size"]-3 == len(apu.predictions.loc[data[0].iloc[i].loc["ID:transcript id"]].loc["predictions"])
                              else False for i in range(data[0].shape[0])]]

        print("< filtering of inconsistent cds sizes reduced data from", init_shape, "to", data[0].shape[0])
        data[0].to_csv(path_or_buf=data_dir+"\\"+fnames[0].split(".")[0]+"_hg_filtered.txt", sep=",", index=False)


    if mode == "filter_shared_genes":
        shared_genes_total    = data[0][data[0]["ID:gene id"].isin(data[1]["ID:gene id"])]
        data[0]               = data[0][~data[0]["ID:gene id"].isin(data[1]["ID:gene id"])]
        test1 = data[0].drop_duplicates(subset=["ID:cell id"]).shape[0]
        test2 = data[0].drop_duplicates(subset=["ID:sample id"]).shape[0]
        test3 = data[0].drop_duplicates(subset=["ID:gene id"]).shape[0]

        print("< shared genes", shared_genes_total.shape[0], "unshared genes", data[0].shape[0], "tests", test1, "/", test2, "/", test3)
        data[0].to_csv(path_or_buf=data_dir+"\\"+fnames[0].split(".")[0]+"_filtered_genes.txt", sep=",", index=False)


    if mode == "filter_shared_variants":
        data[0]["varname"]    = [data[0].iloc[i].loc["ID:gene id"] + "_" + str(data[0].iloc[i].loc["FEATURE:ptc cds position"]) for i in range(data[0].shape[0])]
        data[1]["varname"]    = [data[1].iloc[i].loc["ID:gene id"] + "_" + str(data[1].iloc[i].loc["FEATURE:ptc cds position"]) for i in range(data[1].shape[0])]
        shared_variants_total = data[0][data[0]["varname"].isin(data[1]["varname"])]
        shared_variants       = shared_variants_total.drop_duplicates(subset="varname")["varname"].tolist()
        data[0]               = data[0][~data[0]["varname"].isin(data[1]["varname"])]
        test1 = data[0].drop_duplicates(subset=["ID:cell id"]).shape[0]
        test2 = data[0].drop_duplicates(subset=["ID:sample id"]).shape[0]
        test3 = data[0].drop_duplicates(subset=["ID:variant id"]).shape[0]

        print("< shared variants", shared_variants_total.shape[0], "/", len(shared_variants), "unshared variants", data[0].shape[0], "tests", test1, "/", test2, "/", test3)
        data[0].to_csv(path_or_buf=data_dir+"\\"+fnames[0].split(".")[0]+"_filtered.txt", sep=",", index=False)


    # conduct correlation test to estimate maximum achievable prediction performance for NMD prediction using averaged NMD scores as predictors
    if mode == "get_block_correlation":
        blocks      = {"ID:variant id": [], "LABEL:NMD score": [], "LABEL:avg. NMD score": []}
        data[0]     = data[0].sort_values(by=["ID:variant id"])
        block_index = create_blocks(data[0], ["ID:variant id"])

        test_count  = 0
        bar = IncrementalBar(set_bar("calculating block correlation"), max=len(block_index))
        for i in range(len(block_index)):
            outliers  = [data[0].iloc[j].loc["LABEL:NMD score"] for j in block_index[i]["index"]
                         if data[0].iloc[j].loc["LABEL:NMD score"] < 0 or data[0].iloc[j].loc["LABEL:NMD score"] > 1]

            if len(outliers) > 0:
                print("< outliers found @get_block_correlation for variant", data[0].iloc[block_index[i]["index"][0]].loc["ID:variant id"])
                print(outliers)

            else:
                avg_score = np.mean([data[0].iloc[j].loc["LABEL:NMD score"] for j in block_index[i]["index"]])

                for j in block_index[i]["index"]:
                    # print("i", i, "j", j, data[0].iloc[j].loc["ID:variant id"], data[0].iloc[j].loc["LABEL:NMD score"], avg_score)
                    blocks["ID:variant id"].append(data[0].iloc[j].loc["ID:variant id"])
                    blocks["LABEL:NMD score"].append(data[0].iloc[j].loc["LABEL:NMD score"])
                    blocks["LABEL:avg. NMD score"].append(avg_score)
                    test_count += 1

                    if data[0].iloc[j].loc["ID:variant id"] != data[0].iloc[block_index[i]["index"][0]].loc["ID:variant id"]:
                        print("< block creation error1 with", data[0].iloc[j].loc["ID:variant id"], "and", data[0].iloc[block_index[i]["index"][0]].loc["ID:variant id"])
            bar.next()
        bar.finish()

        if test_count != data[0].shape[0]: print("< block creation error2 with", test_count, "/", data[0].shape[0])

        print("<", len(block_index), "blocks detected. pearson", scipy.stats.pearsonr(blocks["LABEL:NMD score"], blocks["LABEL:avg. NMD score"]))
        blocks = pd.DataFrame(blocks)
        blocks.to_csv(path_or_buf=data_dir+"\\"+fnames[0].split(".")[0]+"_block_correlation.txt", sep=",", index=False)


    if mode == "merge_df":
        test = 0
        for i in range(len(data)):
            test += data[i].shape[0]

        df = pd.concat(data, ignore_index=True)
        print("< shapes", df.shape[0], "/", test)
        df.to_csv(path_or_buf=data_dir+os_sep+"hg38_NMD_scores.txt", sep=",", index=False)


    if mode == "merge_expression_data":
        merged_df = pd.concat(data, axis=0, ignore_index=True)
        merged_df.to_csv(path_or_buf=data_dir+os_sep+"merged_expression_data.txt", sep=",", index=False)


    if mode == "randomize_cohorts":
        targets = ["ID:gene symbol"] # ["ID:gene id"]
        rand_df = get_random_cohorts(data[0], cohorts, targets)
        # test randomization
        for i in range(cohorts):
            df1 = rand_df[rand_df["ID:cohort"] == i+1]
            df2 = rand_df[rand_df["ID:cohort"] != i+1]
            if df1[df1[targets[0]].isin(df2[targets[0]].tolist())].shape[0] > 0:
                print("< gene symbol contamination @cohort", i+1)

        rand_df.to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_randomized.txt", sep=",", index=False)


    # program to rearrange the mutations (e.g. nmd_mutations) output from extract_mutations to generate input for analyze_targets (checked 250623)
    if mode == "rearrange_mutations":
        key = "nmd" # "hla" "nmd"

        if data[0][~data[0]["Gene"].isin(data[1]["gene_id"])].shape[0] > 0:
            print("< error occurred. unknown gene ids detected.")
        
        # statistics containers
        case_ids   = data[0].drop_duplicates(subset="case_id")["case_id"].tolist()
        case_stats = pd.DataFrame({"case id": [case_id for case_id in case_ids], "total_nmd": [0 for _ in case_ids], "frameshift_nmd": [0 for _ in case_ids],
                                   "missense_nmd": [0 for _ in case_ids], "nonsense_nmd": [0 for _ in case_ids]}, index=[case_id for case_id in case_ids])
        
        incidences       = 0
        bar              = IncrementalBar(set_bar("selecting mutations"), max=data[0].shape[0])
        for i in range(data[0].shape[0]):
            case_id = data[0].iloc[i].loc["case_id"]
            test    = 0

            # nonsense or frameshift mutations are always considered deleterious
            if data[0].iloc[i].loc["Variant_Classification"] in ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"]:
                case_stats.at[case_id, "total_"+key]          += 1
                test                                          += 1

                if data[0].iloc[i].loc["Variant_Classification"] == "Frame_Shift_Del" or data[0].iloc[i].loc["Variant_Classification"] == "Frame_Shift_Ins":
                    case_stats.at[case_id, "frameshift_"+key] += 1

                elif data[0].iloc[i].loc["Variant_Classification"] == "Nonsense_Mutation":
                    case_stats.at[case_id, "nonsense_"+key]   += 1

            # missense mutations are considered deleterious if mareked as deleterious or probably damaging
            elif (data[0].iloc[i].loc["Variant_Classification"] == "Missense_Mutation"
                and data[0].iloc[i].loc["SIFT"].split("(")[0] == "deleterious" and data[0].iloc[i].loc["PolyPhen"].split("(")[0] == "probably_damaging"):
                case_stats.at[case_id, "missense_"+key]       += 1
                case_stats.at[case_id, "total_"+key]          += 1
                test                                          += 1

            elif (type(data[0].iloc[i].loc["SIFT"]) == str and type(data[0].iloc[i].loc["PolyPhen"]) == str
                    and data[0].iloc[i].loc["SIFT"].split("(")[0] == "deleterious" and data[0].iloc[i].loc["PolyPhen"].split("(")[0] == "probably_damaging"):
                case_stats.at[case_id, "total_"+key]          += 1
                test                                          += 1

                if data[0].iloc[i].loc["Variant_Classification"] == "Silent":
                    print("< warning. mutation ("+data[0].iloc[i].loc["HGVSp"]+") defined as deleterious")

            if test == 1:
                incidences += 1

            if test > 1:
                print("< error @", case_id)
                exit()

            bar.next()
        bar.finish()
        print("<", case_stats.shape[0], "selected cases with", incidences, "incidences")
        case_stats.to_csv(path_or_buf=data_dir+os_sep+fnames[0].split(".")[0]+"_rearranged.txt", sep=",", index=False)


    # tool to create mutation stats with gene factors from msk dataset
    # data: "tcga_mutation_stats.json", "msk_chord_data_mutations.txt", "msk_chord_survival_analysis.txt", "hg38_seqs.txt"
    if mode == "replace_mutation_stats":
        exon_mutations = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                          "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Site"]
        
        data[1]        = data[1][data[1]["Tumor_Sample_Barcode"] == "P-0010783-T03-IM6"]
        
        cancer_types   = np.unique(data[2]["ID:CANCER_TYPE"])
        init_shape     = data[2].shape[0]
        data[2]        = data[2][~data[2]["ID:specimen time"].isna()]
        print("< removal of data entries without sequence info reduced dataset from", init_shape, "to", data[2].shape[0])

        # initialize case container
        cases = {cancer_type: {} for cancer_type in [*cancer_types, "total"]}
        
        for key in cases:
            if key != "total": patients = np.unique(data[2][data[2]["ID:CANCER_TYPE"] == key]["ID:PATIENT_ID"])
            else:              patients = np.unique(data[2]["ID:PATIENT_ID"])
            cases[key] = {patient: 0 for patient in patients}

        # initialize gene container
        gene_factors   = {cancer_type: {} for cancer_type in [*cancer_types, "total"]}

        data[1]["patient id"] = [data[1].iloc[i].loc["Tumor_Sample_Barcode"].split("-")[0]+"-"+data[1].iloc[i].loc["Tumor_Sample_Barcode"].split("-")[1]
                                 for i in range(data[1].shape[0])]
        
        data[3]["transcript id"] = [transcript_id.split(".")[0] for transcript_id in data[3]["transcript id"]]
        data[3].index            = data[3]["transcript id"]

        # remove duplicates
        if data[1].drop_duplicates(subset=["patient id", "Transcript_ID", "HGVSc"]).shape[0] != data[1].shape[0]:
            print("<", data[1].shape[0]-data[1].drop_duplicates(subset=["patient id", "Transcript_ID", "HGVSc"]).shape[0], "duplicates found.")
            data[1] = data[1].drop_duplicates(subset=["patient id", "Transcript_ID", "HGVSc"])

        # filter by mutation entries
        data[1] = data[1][data[1]["Variant_Classification"].isin(exon_mutations)]

        bar = IncrementalBar(set_bar("replacing gene factors"), max=data[1].shape[0])
        for i in range(data[1].shape[0]):
            patient_id = data[1].iloc[i].loc["Tumor_Sample_Barcode"].split("-")[0]+"-"+data[1].iloc[i].loc["Tumor_Sample_Barcode"].split("-")[1]
            selected   = data[2][data[2]["ID:PATIENT_ID"] == patient_id]
            
            if selected.shape[0] > 0:
                cancer_type   = selected.iloc[0].loc["ID:CANCER_TYPE"]
                transcript_id = data[1].iloc[i].loc["Transcript_ID"]

                if transcript_id in data[3].index:
                    seq_size = len(data[3].loc[transcript_id].loc["cds"])

                    if transcript_id in gene_factors[cancer_type]: gene_factors[cancer_type][transcript_id] += 1/seq_size
                    else:                                          gene_factors[cancer_type][transcript_id]  = 1/seq_size

                    if transcript_id in gene_factors["total"]:     gene_factors["total"][transcript_id]     += 1/seq_size
                    else:                                          gene_factors["total"][transcript_id]      = 1/seq_size

                cases[cancer_type][patient_id] += 1

            else:
                print("<", patient_id, "not found.")

            bar.next()
        bar.finish()

        data[0]["cases"] = cases
        data[0]["genes"] = gene_factors

        for key in data[0]["genes"]:
            data[0]["avg. gene factor"][key] = np.mean([data[0]["genes"][key][gene] for gene in data[0]["genes"][key]])
        
        with open(data_dir+os_sep+"msk_mutation_stats.json", "w") as f:
            f.write(json.dumps(data[0], indent=4))

    
    # for cuomo data extracted in appris=redundant mode, variants that are represented with multiple transcript ids (with the same appris level) are reduced to only one
    if mode == "select_appris":
        data[0]     = data[0].sort_values(by="ID:variant id")
        block_index = create_blocks(data[0], ["ID:variant id"])

        selected_data = []

        for i in range(len(block_index)):
            block_data            = data[0].iloc[block_index[i]["index"]]
            unique_transcript_ids = np.unique(block_data["ID:transcript id"])
            selected_data.append(block_data[block_data["ID:transcript id"] == unique_transcript_ids[0]])
            #print("<", i, block_index[i]["block id"], len(unique_transcript_ids), block_data[block_data["ID:transcript id"] == unique_transcript_ids[0]].shape[0])

            if len(unique_transcript_ids) == 1 and len(block_index[i]["index"]) != block_data[block_data["ID:transcript id"] == unique_transcript_ids[0]].shape[0]:
                print("< inconsistent size.", len(block_index[i]["index"]), "/", block_data[block_data["ID:transcript id"] == unique_transcript_ids[0]].shape[0])

        selected_data = pd.concat(selected_data)
        print("< selecting for single transcript reduced data from", data[0].shape[0], "to", selected_data.shape[0])
        selected_data.to_csv(data_dir+os_sep+fnames[0].split(".")[0]+"_appris_selected.txt")


    # input data e.g. tcga_variants.txt, Census_allThu Jul 10 11_18_52 2025.csv (datatype=census) or TCGA.PanCancer.onco.genes.OncoVar.tsv (datatype=pan_cancer)
    if mode == "split_by_driver_genes":
        datatype  = "pan_cancer" # "pan_cancer" "census"

        if datatype == "census":
            data[1]   = data[1][data[1]["Molecular Genetics"] == "Dom"]
            oncogenes = data[1][data[1]["Role in Cancer"] == "oncogene"]["Gene Symbol"]
            tsgs      = data[1][data[1]["Role in Cancer"] == "TSG"]["Gene Symbol"]

        if datatype == "pan_cancer":
            oncogenes = data[1][data[1]["OncoKB"] == "Oncogene"]["Gene_symbol"]
            tsgs      = data[1][data[1]["OncoKB"] == "TSG"]["Gene_symbol"]

        oncogene_data = data[0][data[0]["ID:gene symbol"].isin(oncogenes)]
        oncogene_data.to_csv(data_dir+os_sep+fnames[0].split(".")[0]+"_oncogenes.txt")
        tsgs_data     = data[0][data[0]["ID:gene symbol"].isin(tsgs)]
        tsgs_data.to_csv(data_dir+os_sep+fnames[0].split(".")[0]+"_tsgs.txt")


if __name__ == '__main__':
    main()