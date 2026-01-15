import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from progress.bar import IncrementalBar
# marked (<-) added on 250429 to allow integration of hg38 for Cuomo data in combination with liftover
from pyliftover import LiftOver # <-
import scipy
from sklearn import metrics
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


chrs = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11,
        "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "20": 19, "21": 20, "22": 21, "X": 22, "Y": 23}


def create_lindeboom_dictionary(lines):
    gene_dictionary = {}

    exception = "numeric" # exception found in Lindeboom preditions for hg19
    bar   = IncrementalBar(set_bar("creating Lindeboom dictionary"), max=len(lines)/1000)
    steps = 0

    for line in lines:
        if "CDS" in line:
            residual = line.split("\t")
            #if residual[0][3::] in chrs: chr = chrs[residual[0][3::]]
            #else:                        chr = -1
            chr      = residual[0]
            cdsstart = residual[3]
            cdsend   = residual[4]
            strand   = residual[6]
            gene_id  = residual[8].split(" ")[1].strip(";")
            gene_id  = gene_id.split(".")[0] # remove version number
            scores   = residual[8].split(" ")[9].strip(";\n").split(",")
            scores   = [float(score) for score in scores if exception not in score]

            if strand != "-" and strand != "+":
                print("strand error occurred @line:")
                print(line)

            else:
                if gene_id not in gene_dictionary.keys():  gene_dictionary[gene_id] = [{"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                        "strand": strand, "scores": scores}]
                else:                                      gene_dictionary[gene_id].append({"cdsstart": int(cdsstart), "cdsend": int(cdsend), "chr": chr,
                                                                                            "strand": strand, "scores": scores})

        if steps % 1000 == 0: bar.next()
        steps += 1

    bar.finish()
    return gene_dictionary

    
def evaluate(targets, params):
    if params["datatype"] == "custom" or params["datatype"] == "tcga":
        # marked (<-) added / removed on 250509 due to conversion error
        #labels = [1/(2*math.exp(-targets.iloc[i].loc["LABEL:NMD score"]*math.log(2)))
        #          if targets.iloc[i].loc["LABEL:NMD score"] != None and targets.iloc[i].loc["LABEL:NMD score"] <= 0 else None
        #          for i in range(targets.shape[0])] # <- removed
        labels = convert_labels(targets["LABEL:NMD score"].tolist()) # <- added

    # marked (<-) newly added on 250410 to integrate data from Kim et al.
    if params["datatype"] == "kim": # <-
        labels = targets["LABEL:NMD score"].tolist() # <-
    
    # marked (<-) added/removed on 250429 to allow integration of hg38 in combination with liftover
    # if params["datatype"] == "cuomo": # <- removed
    if "cuomo" in params["datatype"]: # <- added
        labels = targets["LABEL:NMD score"].tolist()

    if params["datatype"] == "teran":
        labels = targets["REF_RATIO"].tolist()

    predictions = targets["score"].tolist()

    test = [label for label in labels if label < 0] # <- added on 250616
    if len(test) > 0: # <- added on 250616
        print("< error. negative label value(s) detected.") # <- added on 250616
        print(test) # <- added on 250616
        exit() # <- added on 250616

    results = {
              # -log2( [A+B] / [2A]) as in Lindeboom et al. 2016
              "nonASE (log)" : {"labels":      [-math.log2(1 / (2*labels[i])) for i in range(len(labels))
                                                # marked (<-) was edited on 250409 as previous version caused error
                                                # if np.isnan(labels[i]) == False and labels[i] != 0 and np.isnan(predictions[i]) == False], # <- replaced
                                                if pd.isna(labels[i]) == False and labels[i] > 0 and pd.isna(predictions[i]) == False], # <- new
                                "predictions": [predictions[i] for i in range(len(predictions))
                                                # marked (<-) was edited on 250409 as previous version caused error
                                                # if np.isnan(labels[i]) == False and labels[i] != 0 and np.isnan(predictions[i]) == False], # <- replaced
                                                if pd.isna(labels[i]) == False and labels[i] > 0 and pd.isna(predictions[i]) == False], # <- new
                                "log":          True},

              # A / [A+B], as in Teran et al. 2021
              # r = -log2( 1 / 2A )
              # r = -ln(1 / 2A) / ln(2)
              # -r*ln(2) = ln(1 / 2A)
              # exp(-r*ln(2)) = 1/2A
              # A = 1/(2*exp(-r*ln(2))
              # changed 240917
              "ASE"         :   {"labels":      [labels[i] for i in range(len(labels))
                                                 # marked (<-) was edited on 250409 as previous version caused error
                                                 # if np.isnan(labels[i]) == False and np.isnan(predictions[i]) == False], # <- replaced
                                                 if pd.isna(labels[i]) == False and pd.isna(predictions[i]) == False], # <- new
                                 "predictions": [1/(2*math.exp(-predictions[i]*math.log(2))) for i in range(len(predictions))
                                                 # marked (<-) was edited on 250409 as previous version caused error
                                                 # if np.isnan(labels[i]) == False and np.isnan(predictions[i]) == False], # <- replaced
                                                 if pd.isna(labels[i]) == False and pd.isna(predictions[i]) == False], # <- new
                                 "log":         False}

              }
    
    fig, ax = plt.subplots(int(len(results)/2)+1, 2)
    col = 0; row = 0
    for result in results:
        if results[result]["log"] == False: decision_threshold = params["decision_threshold"]
        else:                               decision_threshold = -math.log2(params["decision_threshold"])
        r              = scipy.stats.pearsonr(results[result]["labels"], results[result]["predictions"]).statistic
        r2             = r**2
        rmse           = metrics.mean_squared_error(results[result]["labels"], results[result]["predictions"])

        balanced_df, _ = balance(pd.DataFrame({"label": results[result]["labels"], "prediction": results[result]["predictions"]}),
                                               "label", threshold=decision_threshold, randomize=True)
        # could theoritically lead to errors if labels contain None values (excluded for Teran and Cuomo data)
        cat_labels     = [1 if balanced_df.iloc[i].loc["label"] >= decision_threshold else 0 for i in range(balanced_df.shape[0])]

        auroc          = metrics.roc_auc_score(cat_labels, balanced_df["prediction"].tolist())

        print("<", result, "r:", round(r, 4), ", r2:", round(r2, 4), ", rmse:", round(rmse, 4),
              "balanced auroc", round(auroc, 4), "no. of values", len(results[result]["predictions"]))
        
        ax[row][col].scatter(results[result]["labels"], results[result]["predictions"])
        ax[row][col].set_title(result)
        col += 1
        if col == 2: col = 0; row += 1
    
    plt.show()


def get_absolute_index_(lindeboom_dictionary, ucid, abs_ptc_index):
    cds_length       = 0
    ptc_match_index  = -1
    ptc_match_length = 0
    strand           = lindeboom_dictionary[ucid][0]["strand"]

    for i in range(len(lindeboom_dictionary[ucid])):
        if (abs_ptc_index >= lindeboom_dictionary[ucid][i]["cdsstart"]
            and abs_ptc_index < lindeboom_dictionary[ucid][i]["cdsend"]):
            ptc_match_index  = i
            ptc_match_length = cds_length

        cds_length += lindeboom_dictionary[ucid][i]["cdsend"]-lindeboom_dictionary[ucid][i]["cdsstart"]+1

    if strand == "+":
        cds_index = (abs_ptc_index-(lindeboom_dictionary[ucid][ptc_match_index]["cdsstart"]-1))+ptc_match_length
    
    if strand == "-":
        cds_index = (cds_length-ptc_match_length)-1-(abs_ptc_index-lindeboom_dictionary[ucid][ptc_match_index]["cdsstart"])-1

    return cds_index, cds_length, ptc_match_index


def get_relative_index_(lindeboom_dictionary, ucid, ptc_index):
    abs_mutation_index = -1
    cds_length         = 0
    ptc_match_index    = -1
    ptc_cds_exon_index = -1
    strand             = lindeboom_dictionary[ucid][0]["strand"]
    
    if strand == "+":
        for i in range(len(lindeboom_dictionary[ucid])):
            exon_cds_length = lindeboom_dictionary[ucid][i]["cdsend"]-lindeboom_dictionary[ucid][i]["cdsstart"]+1
            
            if ptc_index >= cds_length and ptc_index < cds_length+exon_cds_length:
                ptc_match_index    = i
                ptc_cds_exon_index = ptc_index-cds_length
                abs_mutation_index = lindeboom_dictionary[ucid][len(lindeboom_dictionary[ucid])-1-i]["cdsstart"]-1+ptc_cds_exon_index

            cds_length += exon_cds_length

    if strand == "-":
        for i in range(len(lindeboom_dictionary[ucid])):
            exon_cds_length = lindeboom_dictionary[ucid][len(lindeboom_dictionary[ucid])-1-i]["cdsend"]-lindeboom_dictionary[ucid][len(lindeboom_dictionary[ucid])-1-i]["cdsstart"]+1
            
            if ptc_index >= cds_length and ptc_index < cds_length+exon_cds_length:
                ptc_match_index    = len(lindeboom_dictionary[ucid])-1-i
                ptc_cds_exon_index = exon_cds_length-(ptc_index-cds_length)-1
                abs_mutation_index = lindeboom_dictionary[ucid][len(lindeboom_dictionary[ucid])-1-i]["cdsstart"]-1+ptc_cds_exon_index

            cds_length += exon_cds_length
    
    return abs_mutation_index, ptc_cds_exon_index, cds_length, ptc_match_index


def get_lindeboom_predictions_(scores, report, targets, lindeboom_dictionary, genome, params, it):
    # with insertion, no wt base is given (marked as "-")
    # marked (<-) newly added/removed on 250410 to integrate data from Kim et al.
    # if params["datatype"] == "cuomo" or params["datatype"] == "custom" or params["datatype"] == "tcga": # <- removed
    if params["datatype"] != "teran": # <- added
        # for non-frameshifting variants, abs. mutation index and abs. ptc index are identical
        # marked (<-) added/removed on 250429 to allow calculation of cuomo hg38 along with hg19
        # if params["datatype"] == "cuomo": # <- removed
        if params["datatype"] == "cuomo19": # <- added
            abs_mutation_index = int(targets.iloc[it].loc["ID:variant id"].split("_")[1])-1
            # marked (<-) changed on 250429
            # abs_ptc_index      = int(targets.iloc[it].loc["ID:variant id"].split("_")[1])-1 <- removed
            abs_ptc_index      = targets.iloc[it].loc["FEATURE:abs. ptc index"] # <- added
            ref_base           = targets.iloc[it].loc["ID:variant id"].split("_")[2]

        # marked (<-) added on 250429 to allow calculation of cuomo hg38 along with hg19
        if params["datatype"] == "cuomo38": # <- added
            abs_mutation_index = int(targets.iloc[it].loc["ID:abs. mutation index"]) # <- added
            abs_ptc_index      = targets.iloc[it].loc["FEATURE:abs. ptc index"] # <- added
            ref_base           = targets.iloc[it].loc["ID:variant id"].split("_")[2] # <- added

        if params["datatype"] == "custom":
            abs_mutation_index = int(targets.iloc[it].loc["ID:abs. mutation index"])
            abs_ptc_index      = int(targets.iloc[it].loc["FEATURE:abs. ptc index"])
            ref_base           = targets.iloc[it].loc["ID:wt base"][0]

        # marked (<-) newly added/removed on 250410 to integrate data from Kim et al.
        #if params["datatype"] == "tcga": # <- removed
        if params["datatype"] == "kim" or params["datatype"] == "tcga": # <- added
            abs_mutation_index = int(targets.iloc[it].loc["ID:abs. mutation index"])
            abs_ptc_index      = int(targets.iloc[it].loc["FEATURE:abs. ptc index"])
            ref_base           = targets.iloc[it].loc["ID:wt base"][0]

        protein_length = targets.iloc[it].loc["FEATURE:total cds size"]/3
        ptc_index      = targets.iloc[it].loc["FEATURE:ptc cds position"]

    # for non-frameshifting variants, abs. mutation index and abs. ptc index are identical
    if params["datatype"] == "teran":
        abs_mutation_index    = targets.iloc[it].loc["POS"]-1
        if targets.iloc[it].loc["STRAND"] == 1:  frameshift_correction = (3*int((targets.iloc[it].loc["CDS_position_only"]-1)/3))-(targets.iloc[it].loc["CDS_position_only"]-1)
        if targets.iloc[it].loc["STRAND"] == -1: frameshift_correction = (targets.iloc[it].loc["CDS_position_only"]-1)-(3*int((targets.iloc[it].loc["CDS_position_only"]-1)/3))

        abs_ptc_index         = (targets.iloc[it].loc["POS"]-1)+frameshift_correction # targets.iloc[it].loc["POS"]-1
        protein_length        = targets.iloc[it].loc["Protein_length"]
        ptc_index             = 3*int((targets.iloc[it].loc["CDS_position_only"]-1)/3) # targets.iloc[it].loc["CDS_position_only"]-1
        ref_base              = targets.iloc[it].loc["REF_ALLELE"]

    # iterate through ucids, check for which splice variant both mutation index and protein length fit
    temp_scores = []
    info        = [str(targets.iloc[it].loc[params["block_identifier"]])]
    for i in range(len(targets.iloc[it].loc["uc id"])):
        info.append("it " + str(it) + " i " + str(i) + " " + str(targets.iloc[it].loc["uc id"][i]))
        ucid = targets.iloc[it].loc["uc id"][i].split(".")[0] # remove version number

        if ucid in lindeboom_dictionary.keys():
            chr              = lindeboom_dictionary[ucid][0]["chr"]
            strand           = lindeboom_dictionary[ucid][0]["strand"]

            if chr not in genome:
                print("< error @get_lindeboom_predictions_. Chromosome", chr, "not found in genome.")

            if params["mode"] == "absolute": cds_index, cds_length, ptc_match_index                              = get_absolute_index_(lindeboom_dictionary, ucid, abs_ptc_index)
            if params["mode"] == "relative": abs_mutation_index, ptc_cds_exon_index, cds_length, ptc_match_index = get_relative_index_(lindeboom_dictionary, ucid, ptc_index)
            # marked (<-) newly added/removed on 250410 to integrate data from Kim et al.
            #if params["mode"] == "absolute" and (params["datatype"] == "cuomo" or params["datatype"] == "custom" or params["datatype"] == "tcga"): cds_index -= 1 # <- removed
            if params["mode"] == "absolute" and params["datatype"] != "teran": cds_index -= 1 # <- added

            info.append("  cds_length " + str(cds_length) + " " + str(3*protein_length) + " " + str(ptc_match_index))
            
            # marked (<-) newly added/removed on 250410 to integrate data from Kim et al.
            #if ((((params["datatype"] == "cuomo" or params["datatype"] == "custom" or params["datatype"] == "tcga") and cds_length == 3*protein_length-3) # <- removed
            if (((params["datatype"] != "teran" and cds_length == 3*protein_length-3) # <- added
                 or (params["datatype"] == "teran" and cds_length == 3*protein_length))
                and ptc_match_index != -1):
                # test section
                tests_passed = 0

                # check whether cds index is identical with ptc index (works only with absolute index)
                if "absolute" in params["mode"]:
                    if cds_index == ptc_index: tests_passed += 1
                    else:                      info.append("< reconstructed cds index and ptc index do not match: " + str(cds_index) + " / " + str(ptc_index) + ", strand: " + strand)

                if "relative" in params["mode"]:
                    tests_passed += 1

                # check whether the correct wt base is found in genome (test not possible for insertion as ref_base = "-")
                if ref_base != "-" and ref_base == genome[chr][abs_mutation_index]:
                    tests_passed += 1

                elif ref_base == "-":
                    tests_passed += 1

                else:
                    info.append("wt base is not found in genome. wt base: " + ref_base + ", context: " + genome[chr][abs_mutation_index-4:abs_mutation_index+4])

                # check whether the correct start codon is found in genome
                if strand == "+":
                    cdsstart = lindeboom_dictionary[ucid][0]["cdsstart"]-1

                    info.append("  start, + " + genome[chr][cdsstart-2:cdsstart+5] + genome[chr][cdsstart:cdsstart+3])
                    if "ATG" == genome[chr][cdsstart:cdsstart+3]: tests_passed += 1
                    else:                                         info.append("start codon is not found in genome.")

                if strand == "-":
                    cdsend   = lindeboom_dictionary[ucid][-1]["cdsend"]
                    info.append("  start, - " + genome[chr][cdsend-5:cdsend+2] + genome[chr][cdsend-3:cdsend])
                    if "CAT" == genome[chr][cdsend-3:cdsend]:     tests_passed += 1
                    else:                                         info.append("start codon is not found in genome.")

                # check whether the correct stop codon is found in genome
                if strand == "+":
                    cdsend   = lindeboom_dictionary[ucid][-1]["cdsend"]+3
                    info.append("  stop, + " + genome[chr][cdsend-5:cdsend+2] + genome[chr][cdsend-3:cdsend])
                    if genome[chr][cdsend-3:cdsend] in ["TAA", "TAG", "TGA"]:     tests_passed += 1
                    else:                                                         info.append("stop codon is not found in genome.")

                if strand == "-":
                    cdsstart = lindeboom_dictionary[ucid][0]["cdsstart"]-4
                    info.append("  stop, - " + genome[chr][cdsstart-5:cdsstart+7] + genome[chr][cdsstart:cdsstart+3])
                    if genome[chr][cdsstart:cdsstart+3] in ["TTA", "CTA", "TCA"]: tests_passed += 1
                    else:                                                         info.append("stop codon is not found in genome.")

                if params["mode"] == "absolute" and tests_passed != 4:
                    print("  error occured @", targets.iloc[it].loc[params["block_identifier"]], "associated ucid entry:", targets.iloc[it].loc["uc id"], "tests_passed", tests_passed)
                
                # determine the relative score position and check whether the position equals the wt base
                else:
                    if strand == "+":
                        if params["mode"] == "absolute": exon_pos = abs_ptc_index-(lindeboom_dictionary[ucid][ptc_match_index]["cdsstart"]-1)
                        if params["mode"] == "relative": exon_pos = ptc_cds_exon_index

                        if exon_pos < len(lindeboom_dictionary[ucid][ptc_match_index]["scores"]):
                            temp_scores.append(lindeboom_dictionary[ucid][ptc_match_index]["scores"][exon_pos])

                        else:
                            print("< error occured @", targets.iloc[it].loc[params["block_identifier"]], "exon position exceeds score dimension",
                                  exon_pos, "/", len(lindeboom_dictionary[ucid][ptc_match_index]["scores"]), "strand:", strand)
                            info.append("< error occured @" + targets.iloc[it].loc[params["block_identifier"]] + ". exon position exceeds score dimension "
                                        + str(exon_pos), " / ", str(len(lindeboom_dictionary[ucid][ptc_match_index]["scores"])) + ", strand: " + strand)
                        
                    if strand == "-":
                        if "absolute" in params["mode"]: exon_pos = abs_ptc_index-(lindeboom_dictionary[ucid][ptc_match_index]["cdsstart"]-1) #lindeboom_dictionary[ucid][ptc_match_index]["cdsend"]-abs_ptc_index
                        if "relative" in params["mode"]: exon_pos = ptc_cds_exon_index

                        if exon_pos < len(lindeboom_dictionary[ucid][ptc_match_index]["scores"]):
                            temp_scores.append(lindeboom_dictionary[ucid][ptc_match_index]["scores"][exon_pos])

                        else:
                            print("< error occured @", targets.iloc[it].loc[params["block_identifier"]], "exon position exceeds score dimension",
                                  exon_pos, "/", len(lindeboom_dictionary[ucid][ptc_match_index]["scores"]), "strand:", strand)
                            info.append("< error occured @" + targets.iloc[it].loc[params["block_identifier"]] + ". exon position exceeds score dimension "
                                        + str(exon_pos), " / ", str(len(lindeboom_dictionary[ucid][ptc_match_index]["scores"])) + ", strand: " + strand)

        # marked (<-) added on 250429 to register missing ucids
        else: # <- added
            print("<", ucid, "not found") # <- added


    if len(temp_scores) == 0:
        report["missing scores"] += 1

        # print error report if ucid was found
        if len(info) > 0:
            print("< error report:")
            for i in info:
                print("  " + i)

    if len(temp_scores) > 1: report["multiple scores"] += 1
    if len(temp_scores) > 0: scores[it] = np.mean(np.array(temp_scores))
    return scores, report


def get_lindeboom_predictions(targets, lindeboom_dictionary, genome, params):
    report  = {"missing scores": 0, "multiple scores": 0}
    targets = targets.sort_values(by=params["block_identifier"])
    scores  = [None for _ in range(targets.shape[0])]

    bar = IncrementalBar("retrieving Lindeboom predictions", max=targets.shape[0])

    if params["blockwise"] == True:
        last_score       = None
        last_variant_id  = None

        for i in range(targets.shape[0]):
            if last_variant_id == None or targets.iloc[i].loc[params["block_identifier"]] != last_variant_id:
                scores, report = get_lindeboom_predictions_(scores, report, targets, lindeboom_dictionary, genome, params, i)

            else:
                scores[i] = last_score

            last_score      = scores[i]
            last_variant_id = targets.iloc[i].loc[params["block_identifier"]]

    else:
        for i in range(targets.shape[0]):
            scores, report = get_lindeboom_predictions_(scores, report, targets, lindeboom_dictionary, genome, params, i)

        #bar.next()
    #bar.finish()
    targets["score"] = scores
    return targets, report


# marked (<-) added on 250429 to allow integration of hg38 in combination with liftover
def _make_liftover(liftover, chr, pos): # <- from here
    liftover_results = liftover.convert_coordinate(chr, pos)
                
    if len(liftover_results) == 1: return liftover_results[0][1] # check for valid, unequivocal results
    else:                          return None # <- until here


# marked (<-) added on 250429 to allow integration of hg38 in combination with liftover
def make_liftover(targets): # <- from here
    liftover = LiftOver('hg19', 'hg38')
    targets["ID:abs. mutation index"] = [_make_liftover(liftover, "chr"+targets.iloc[i].loc["ID:variant id"].split("_")[0],
                                                        int(targets.iloc[i].loc["ID:variant id"].split("_")[1])) for i in range(targets.shape[0])]

    return targets # <- until here


def print_results(targets, params):
    # additional id needed for cuomo data
    # marked (<-) added on 250429 to allow integration of hg38 in combination with liftover
    # if params["datatype"] == "cuomo": # <- removed
    if "cuomo" in params["datatype"]:
        df = pd.DataFrame({"ID:variant id": targets[params["block_identifier"]],
                           "ID:variant transcript pos id": [targets.iloc[i].loc[params["block_identifier"]] + "_" + targets.iloc[i].loc["ID:transcript id"] + "_" + str(targets.iloc[i].loc["FEATURE:ptc cds position"])
                                                            for i in range(targets.shape[0])],
                           "LABEL:NMD score": targets["LABEL:NMD score"],
                           "FEATURE:lindeboom prediction": targets["score"]})

    # marked (<-) newly added/removed on 250410 to integrate data from Kim et al.
    #if params["datatype"] == "tcga": # <- removed
    if params["datatype"] == "kim" or params["datatype"] == "tcga":
        df = pd.DataFrame({"ID:variant id": targets[params["block_identifier"]],
                           "LABEL:NMD score": targets["LABEL:NMD score"],
                           "FEATURE:lindeboom prediction": targets["score"]})
        
    if params["datatype"] == "teran":
        df = pd.DataFrame({"ID:variant id": targets[params["block_identifier"]],
                           "LABEL:NMD score": targets["REF_RATIO"],
                           "FEATURE:lindeboom prediction": targets["score"]})
    
    df = df[~df["FEATURE:lindeboom prediction"].isna()]
    df.to_csv(path_or_buf=params["data_dir"]+params["os_sep"]+params["target_file"].split(".txt")[0]+"_lindeboom_predictions.txt",
              sep=",", index=False)