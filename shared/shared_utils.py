import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
import random
import scipy
from scipy.optimize import curve_fit
import statsmodels.api as sm
import threading

aa_groups       = {
                   "by_aa":
                    {
                    "A":        ["aliphatic"],
                    "C":        ["polar"],
                    "D":        ["acidic"],
                    "E":        ["acidic"],
                    "F":        ["aromatic"],
                    "G":        ["aliphatic"],
                    "H":        ["basic"],
                    "I":        ["aliphatic"],
                    "K":        ["basic"],
                    "L":        ["aliphatic"],
                    "M":        ["polar"],
                    "N":        ["polar"],
                    "V":        ["aliphatic"],
                    "Q":        ["polar"],
                    "P":        ["proline"],
                    "R":        ["basic"],
                    "S":        ["polar"],
                    "T":        ["polar"],
                    "W":        ["aromatic"],
                    "Y":        ["aromatic"],
                    },
                   "by_group":
                    {
                    "acidic":    ["D", "E"],
                    "aliphatic": ["A", "G", "I", "L", "V"],
                    "aromatic":  ["F", "W", "Y"],
                    "basic":     ["H", "K", "R"],
                    "polar":     ["C", "M", "N", "Q", "S", "T"],
                    "proline":   ["P"]
                    }
                  }

alphanumericals = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H",
                   "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

amino_acids     = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

chrs            = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11,
                   "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "20": 19, "21": 20, "22": 21, "X": 22, "Y": 23}

genetic_code    = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                   "TGT": "C", "TGC": "C",
                   "GAT": "D", "GAC": "D",
                   "GAA": "E", "GAG": "E",
                   "TTT": "F", "TTC": "F",
                   "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
                   "CAT": "H", "CAC": "H",
                   "ATT": "I", "ATC": "I", "ATA": "I",
                   "AAA": "K", "AAG": "K",
                   "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                   "ATG": "M",
                   "AAT": "N", "AAC": "N",
                   "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                   "CAA": "Q", "CAG": "Q",
                   "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
                   "AGT": "S", "AGC": "S", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                   "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                   "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                   "TGG": "W",
                   "TAT": "Y", "TAC": "Y",
                   "TAA": "X", "TGA": "X", "TAG": "X"}


def append_df_with_mapping(data, mapping_col1, mapping_col2, target_col, text="", non_redundant=True, pos=4, reverse=False, show_progress=False, verbose=False):
    if len(text) > 0: bar = IncrementalBar(set_bar(text), max=data[0].shape[0])
    if non_redundant == True:  appended_list = ["-" for _ in range(data[0].shape[0])]
    if non_redundant == False: appended_list = [[] for _ in range(data[0].shape[0])]

    variants2      = data[1][mapping_col2].tolist()
    variants2_dict = [[] for _ in range(pow(36, 4))]
    create_dictionary(variants2_dict, variants2, digits=4, alphanumerical=True, non_redundant=non_redundant, reverse=reverse, pos=pos, show_progress=show_progress)
    variants1      = data[0][mapping_col1].tolist()

    no_matches = 0; no_values = 0
    for i in range(len(variants1)):
        match_index = find_match_index(variants2_dict, variants1[i], digits=4, alphanumerical=True, non_redundant=non_redundant, reverse=reverse, pos=pos)

        if type(match_index) == int:
            #if type(data[1].iloc[match_index].loc[target_col]) != str and np.isnan(data[1].iloc[match_index].loc[target_col]) == False:
            if type(data[1].iloc[match_index].loc[target_col]) != str and pd.isna(data[1].iloc[match_index].loc[target_col]) == False:
                appended_list[i] = str(data[1].iloc[match_index].loc[target_col])

            if type(data[1].iloc[match_index].loc[target_col]) == str:
                appended_list[i] = data[1].iloc[match_index].loc[target_col]

            #elif type(data[1].iloc[match_index].loc[target_col]) != str and np.isnan(data[1].iloc[match_index].loc[target_col]) == True:
            elif type(data[1].iloc[match_index].loc[target_col]) != str and pd.isna(data[1].iloc[match_index].loc[target_col]) == True:
                appended_list[i] = "-"
                no_values += 1

        elif type(match_index) == list and len(match_index) > 0:
            for j in match_index:
                #if type(data[1].iloc[j].loc[target_col]) != str and np.isnan(data[1].iloc[j].loc[target_col]) == False:
                if type(data[1].iloc[j].loc[target_col]) != str and pd.isna(data[1].iloc[j].loc[target_col]) == False:
                    appended_list[i].append(str(data[1].iloc[j].loc[target_col]))

                if type(data[1].iloc[j].loc[target_col]) == str:
                    appended_list[i].append(data[1].iloc[j].loc[target_col])

                #elif type(data[1].iloc[j].loc[target_col]) != str and np.isnan(data[1].iloc[j].loc[target_col]) == True:
                elif type(data[1].iloc[j].loc[target_col]) != str and pd.isna(data[1].iloc[j].loc[target_col]) == True:   
                    no_values += 1

        else:
            no_matches += 1
            if verbose == True: print(variants1[i], "not found.")

        if len(text) > 0: bar.next()
        
    if len(text) > 0:
        bar.finish()
        print("< no matches found in", no_matches, "occasions, no values found in", no_values, "occasions.")
        
    return appended_list


def append_df_with_mapping2(data, mapping_col1, mapping_col2, target_col, text="progress", verbose=False):
    bar = IncrementalBar(set_bar(text), max=data[0].shape[0])
    appended_list = [[] for _ in range(data[0].shape[0])]

    variants2      = data[1][mapping_col2].tolist()
    variants2_dict = [[] for _ in range(pow(36, 4))]
    create_dictionary(variants2_dict, variants2, digits=4, alphanumerical=True, non_redundant=False, reverse=True, pos=4)
    
    variants1 = data[0][mapping_col1].tolist()

    no_matches = 0
    for i in range(len(variants1)):
        match_index = find_match_index(variants2_dict, variants1[i], digits=4, alphanumerical=True, non_redundant=False, reverse=True, pos=4)

        if len(match_index) > 0:
            for j in match_index:
                appended_list[i].append(data[1].iloc[j].loc[target_col])

        else:
            no_matches += 1
            if verbose == True: print(variants1[i], "not found.")

        bar.next()
    bar.finish()
    
    print("< no matches found in", no_matches, "of", len(variants1), "occasions.")
    return appended_list


def apply_value_filter(data, value_filter, text):
    if len(value_filter) > 0:
        init_shape = data.shape[0]
        for key in value_filter:
            # marked (<-) added / removed on 250625
            # if key.replace("_EXCLUDED", "").replace("_IDENTICAL", "") not in data.columns: # <- removed
            if key.replace("_EXCLUDED", "").replace("_IDENTICAL", "").replace("_LESS", "").replace("_GREATER", "") not in data.columns: # <- added
                print("< error. filter item", key, "not found in variants.")

            elif type(value_filter[key]) == float or type(value_filter[key]) == int:
                # changed 250103
                if "_IDENTICAL" in key: data = data[data[key.replace("_IDENTICAL", "")] == value_filter[key]]; print("=", key, data[key.replace("_IDENTICAL", "")].min(), data[key.replace("_IDENTICAL", "")].max())
                #if "_IDENTICAL" in key: data = data[data[key.replace("_IDENTICAL", "")] >= value_filter[key]]; print("id", key)
                elif "_LESS" in key:    data = data[data[key.replace("_LESS", "")] < value_filter[key]]; print("<", key, data[key.replace("_LESS", "")].min(), data[key.replace("_LESS", "")].max()) # <- added on 250625
                elif "_GREATER" in key: data = data[data[key.replace("_GREATER", "")] > value_filter[key]]; print(">", key, data[key.replace("_GREATER", "")].min(), data[key.replace("_GREATER", "")].max()) # <- added on 250625
                else:                   data = data[data[key] >= value_filter[key]]; print(">=", key, data[key].min(), data[key].max())

            elif type(value_filter[key]) == list:
                if "_EXCLUDED" in key: data = data[~data[key.replace("_EXCLUDED", "")].isin(value_filter[key])]; print("excluded", key, data[key.replace("_EXCLUDED", "")].min(), data[key.replace("_EXCLUDED", "")].max())
                else:                  data = data[data[key].isin(value_filter[key])]; print("in", key, data[key].min(), data[key].max())

        print("< applying variant value filter to", text, "reduced variants from", init_shape, "to", data.shape[0])
    
    return data


# calculate the probability of drawing a random sub-sample of size sample_size with mean sample_mean out of distribution values
# MAYBE randint is not the correct implementation as no index can occur twice
def bootstrapping_test(values, sample_mean, sample_size, simulation_steps=100):
    means = []
    [means.append(np.mean([values[i] for i in np.random.randint(low=0, high=len(values), size=sample_size)])) for _ in range(simulation_steps)]

    norm = scipy.stats.norm(loc=np.mean(means), scale=np.std(means))
    if sample_mean <= np.mean(means): return norm.cdf(x=sample_mean), np.mean(means), np.std(means)
    if sample_mean > np.mean(means):  return norm.sf(x=sample_mean), np.mean(means), np.std(means)


def balance(df, target, threshold, randomize=True, verbose=True):
    df = df.reset_index(drop=True)
    lower_label_index = df[df[target] < threshold].index.tolist()
    upper_label_index = df[df[target] >= threshold].index.tolist()

    if randomize == True:
        random.shuffle(lower_label_index)
        random.shuffle(upper_label_index)

    min_size = min(len(lower_label_index), len(upper_label_index))

    label_index = lower_label_index[0:min_size] + upper_label_index[0:min_size]
    updated_df  = df.iloc[label_index]

    if verbose == True: print("< balancing reduced size from", df.shape[0], "to", len(label_index))
    return updated_df, label_index


def categorize(data):
    for col in data.columns:
        unique_categories = data.drop_duplicates(subset=col)[col].tolist()
        data[col]         = [unique_categories.index(data.iloc[i].loc[col]) for i in range(data.shape[0])]

    return data


# checks whether pandas dataframe index contains discontinuities
def check_df_index(df):
    test = [i for i in range(df.shape[0]) if df.index[i] != i]
    if len(test) > 0:
        print("< warning. discontinuity found in dataframe index.")
        print(test)
        print(df)
        input("type any character to proceed.")


def check_placeholders(df, cols):
     for col in cols:
        init_shape = df.shape[0]

        if col not in df.columns:
            print("< error.", col, "not found in variants.")
            exit()

        if df[df[col] == -1].shape[0] > 0:
            print("< expired placeholder '-1' detected @", col)
            exit()

        if df[df[col] == "-"].shape[0] > 0:
            print("< expired placeholder '-' detected @", col)
            exit()

def combine_pvalues(pvalues):
    conv_pvalues    = -2 * np.sum(np.log(np.array(pvalues)))
    df              = 2 * len(pvalues)
    return 1-scipy.stats.chi2.cdf(conv_pvalues, df)
            

def conduct_class_test(x, y, threshold, show_hist=False):
    test_matrix = [[] for _ in range(2)]
    for i in range(len(x)):
        if x[i] <= threshold: test_matrix[0].append(y[i]) 
        else:                 test_matrix[1].append(y[i])

    if show_hist == True:
        plt.hist(test_matrix[0], bins=40, density=True, histtype='step', label="category 1")
        plt.hist(test_matrix[1], bins=40, density=True, histtype='step', label="category 2")
        plt.legend()
        plt.show()

    if len(test_matrix[0]) > 0 and len(test_matrix[1]) > 0:
        avg1   = np.mean(test_matrix[0])
        sem1   = np.std(test_matrix[0]) / math.sqrt(len(test_matrix[0]))
        avg2   = np.mean(test_matrix[1])
        sem2   = np.std(test_matrix[1]) / math.sqrt(len(test_matrix[1]))
        pvalue = get_ttest([avg1], [sem1], [len(test_matrix[0])], [avg2], [sem2], [len(test_matrix[1])])[0]
        return test_matrix, [round(avg1/avg2, 4), round(pvalue, 4)]
    
    else:
        return test_matrix, [None, None]


def convert_labels(labels, input_format="nonASE", output_format="ASE"):
    converted_labels = []

    if input_format == "nonASE" and output_format == "ASE":
        converted_labels = [1/(2*math.exp(-math.log(2)*labels[i])) for i in range(len(labels))]
        
    if input_format == "ASE" and output_format == "nonASE":
        converted_labels = [-math.log2(1/(2*labels[i])) if labels[i] > 0 else 0 for i in range(len(labels))]

    return converted_labels


def convert_pos(pos):
    if pos == None:   return 0
    elif pos != None: return pos


def create_blocks(df, block_targets):
    if type(block_targets) != list: block_targets = [block_targets]
    block_index = []

    # added on 250312 to check whether data are not sorted
    if df[block_targets[0]].is_monotonic_increasing == False:
        print("< warning. input blocks are not sorted.")
        exit()

    df = df.sort_values(by=block_targets)
    
    bar            = IncrementalBar(set_bar("identifying block index"), max=df.shape[0])
    last_block_ids = {block_target: None for block_target in block_targets}
    temp           = []

    for i in range(df.shape[0]):
        checks = [block_target for block_target in block_targets if df.iloc[i].loc[block_target] == last_block_ids[block_target]]

        if len(temp) == 0 or len(checks) == len(block_targets):
            temp.append(i)

        elif len(temp) > 0:
            block_index.append({"block id": last_block_ids, "index": temp})
            temp = [i]

        last_block_ids = {block_target: df.iloc[i].loc[block_target] for block_target in block_targets}

        bar.next()
    bar.finish()

    # append last block
    if len(temp) > 0: block_index.append({"block id": last_block_ids, "index": temp})
    print("<", len(block_index), "blocks found")
    return block_index


def create_dictionary(dictionary, ids, digits=4, alphanumerical=False, non_redundant=False, reverse=True, pos=None, show_progress=False):
    if show_progress == True: bar = IncrementalBar(set_bar("creating dictionary"), max=len(ids))
    pos = convert_pos(pos)

    for i in range(len(ids)):
        id = str(ids[i])
        if alphanumerical == False:
            if reverse == False:
                selected_id = int(id[pos:pos+digits])

            else:
                selected_id = int((id[(len(id)-pos-digits):len(id)-pos]))

        else:
            if reverse == False:
                selected_id = get_index(id[pos:pos+digits])

            else:
                try:
                    selected_id = get_index((id[(len(id)-pos-digits):len(id)-pos]))

                except:
                    print("< error occurred for id:", id, "@ index", i)

        if non_redundant == False:
            id_dict = {"id": id, "index": i}
            dictionary[selected_id].append(id_dict)
            #print(ids[i], selected_id, id_dict)
        
        else:
            match = False
            j = 0
            while j < len(dictionary[selected_id]) and match == False:
                if dictionary[selected_id][j]["id"] == id:
                    match = True
                
                j += 1

            if match == False:
                id_dict = {"id": id, "index": i}
                dictionary[selected_id].append(id_dict)
        
        if show_progress == True: bar.next()
    if show_progress == True: bar.finish()

    return dictionary


# marked (<-) added / removed on 250427 to speed up gene detection
# def create_knowngene_dict_by_position(knowngene, stepsize=1000000): # <- removed
def create_knowngene_dict_by_position(knowngene, redundant=False, stepsize=1000000): # <- added
    knowngene_dict = [[] for _ in range(len(chrs))]
    
    bar = IncrementalBar(set_bar("creating KnownGene dictionary"), max=knowngene.shape[0])

    if redundant == False: # <- added (if-statement content not new)
        for i in range(knowngene.shape[0]):
            if knowngene.iloc[i].loc["chr"][3::] in chrs:
                chr     = chrs[knowngene.iloc[i].loc["chr"][3::]]
                exonend = int(knowngene.iloc[i].loc["exonend"])

                converted_exonend = int(exonend/stepsize)
                # marked (<-) added on 250427 to test whether stepsize assures findability of all genes
                if converted_exonend-int(int(knowngene.iloc[i].loc["exonstart"])/stepsize) > 1: # <-
                    print("< warning. with stepsize", stepsize, knowngene.iloc[i].loc["transcript id"], "might not be detected.") # <-

                if converted_exonend < len(knowngene_dict[chr]):
                    knowngene_dict[chr][converted_exonend].append(i)

                else:
                    for _ in range(converted_exonend-len(knowngene_dict[chr])+1):
                        knowngene_dict[chr].append([])

                    knowngene_dict[chr][converted_exonend].append(i)
            else:
                pass

            bar.next()
        bar.finish()

    if redundant == True: # <- added from here
        for i in range(knowngene.shape[0]):
            if knowngene.iloc[i].loc["chr"][3::] in chrs:
                chr = chrs[knowngene.iloc[i].loc["chr"][3::]]

                exonstart = int(knowngene.iloc[i].loc["exonstart"])
                exonend   = int(knowngene.iloc[i].loc["exonend"])

                converted_exonstart = int(exonstart/stepsize)
                converted_exonend   = int(exonend/stepsize)

                for j in range(converted_exonend-converted_exonstart+1):
                    converted_pos = converted_exonstart + j

                    if converted_pos < len(knowngene_dict[chr]):
                        knowngene_dict[chr][converted_pos].append(i)

                    else:
                        for _ in range(converted_pos-len(knowngene_dict[chr])+1):
                            knowngene_dict[chr].append([])

                        knowngene_dict[chr][converted_pos].append(i)

            else:
                pass

            bar.next()
        bar.finish()
        # <- until here

    return knowngene_dict


def find_match_index(dictionary, id, digits=4, alphanumerical=False, non_redundant=True, reverse=True, pos=None):
    id  = str(id)
    pos = convert_pos(pos)
    if non_redundant == True: match_index = None
    else:                     match_index = []

    if alphanumerical == False:
        if reverse == False:
            selected_id = int(id[pos:pos+digits])

        else:
            selected_id = int((id[(len(id)-pos-digits):len(id[i])-pos]))

    else:
        if reverse == False:
            selected_id = get_index(id[pos:pos+digits])

        else:
            selected_id = get_index((id[(len(id)-pos-digits):len(id)-pos]))
        

        match = False
        i = 0
        while i < len(dictionary[selected_id]) and match == False:
            if id == dictionary[selected_id][i]["id"]:
                if non_redundant == True:
                    match_index = dictionary[selected_id][i]["index"]
                    match       = True

                else:
                    match_index.append(dictionary[selected_id][i]["index"])

            i += 1

    return match_index


def fit(x, y, function="linear", components=1):
    fit_y = []

    if function == "exponential":
        parameters, _ = curve_fit(fit_exp, x, y)
        for i in range(len(x)): fit_y.append(fit_exp(x[i], parameters[0], parameters[1]))

    if function == "gaussian":
        init_params = []
        for _ in range(components):
            init_params.append(np.mean(x))
            init_params.append(np.max(x))
            init_params.append(np.std(x)) 

        i_pk  = scipy.signal.find_peaks_cwt(y, widths=range(3, int(len(x)//0.5*components)))
        i_pk = [7, 27, 60] # two-peaks in escape region fitted as single peak
        i_pk = [7, 37, 60] # two-peaks in escape region fitted as single peak
        DX    = (np.max(x)-np.min(x))/float(components) # starting guess for component width
        guess = np.ravel([[x[i], y[i], DX] for i in i_pk]) # starting guess for (x, amp, width) for each component
        print("< guess:", guess)

        parameters, _ = curve_fit(fit_gaussian, x, y, p0=guess, bounds=(0, np.inf))#, p0=init_params)
        fit_y         = [[] for _ in range(int(len(parameters)/3)+1)]

        for i in range(len(fit_y)-1):
            for j in range(len(x)):
                fit_y[i].append(get_gaussian(x[j], parameters[3*i:3*i+3]))
                if i == 0: fit_y[-1].append(get_gaussian(x[j], parameters[3*i:3*i+3]))
                else:      fit_y[-1][j] += get_gaussian(x[j], parameters[3*i:3*i+3])

    if function == "linear":
        parameters, _ = curve_fit(fit_linear, x, y)
        for i in range(len(x)): fit_y.append(fit_linear(x[i], parameters[0], parameters[1]))

    if function == "poisson":
        parameters, _ = curve_fit(fit_poisson, x, y)
        for i in range(len(x)): fit_y.append(fit_poisson(x[i], parameters[0]))

    return parameters, fit_y


def fit_exp(x, A, B):
    y = A*np.exp(-B*x)
    return y


# returns the smallest category into which a value fits in a list of categories (e.g. value 0.007 and [0.01, 0.05, 0.1] -> index 0)
def get_category_index(values, categories):
    category_index      = []
    return_single_value = False
    if type(values) != list: values = [values]; return_single_value = True

    sorted_index = np.argsort(categories)

    for value in values:
        current_category_index = [sorted_index[len(sorted_index)-1-i] for i in range(len(sorted_index)) if value < categories[sorted_index[len(sorted_index)-1-i]]]
        if len(current_category_index) == 0: category_index.append(None)
        else:                                category_index.append(current_category_index[-1])

    if return_single_value == False: return category_index
    if return_single_value == True:  return category_index[0]
            

def fit_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)

    return y


def fit_linear(x, A, B):
    y = A+B*x
    return y


def get_frameshift(info):
    frameshift = 0
    try:
        if "Ter" in info and len(info.split("Ter")[1]) > 0:
            frameshift = float(info.split("Ter")[1])

        elif "Ter" in info and "delins" in info:
            if len(info.split("delins")[1]) % 3 != 0 and len(info.split("delins")[1]) >= 3:
                print("< warning. unknown format @shared_utils, get_frameshift with", info)
            
            else:
                frameshift = len(info.split("delins")[1])/3-1

    except:
        print("frameshift-fail", info, frameshift)
        frameshift = None

    return frameshift


def get_gaussian(x, params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)

    return y


def get_mutation_outcome(protein_change):
    mutation_outcome = "consequential"
    wt_aa      = protein_change[0]
    mutated_aa = protein_change[-1]

    if wt_aa in aa_groups["by_aa"] and mutated_aa in aa_groups["by_aa"] and aa_groups["by_aa"][wt_aa][0] == aa_groups["by_aa"][mutated_aa][0]:
        mutation_outcome = "non_consequential"
        
    return mutation_outcome


def get_normal_transform(x, y):
    mean_x = x.mean()
    std_x  = x.std()

    mean_y = y.mean()
    std_y  = y.std()

    return [((value-mean_y)/std_y)*std_x+mean_x for value in y]


# found on 240827 here: https://stackoverflow.com/questions/22579434/python-finding-the-intersection-point-of-two-gaussian-curves
def get_gaussian_intersection(m1, m2, std1, std2, s1, s2):
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log((std2*s1)/(std1*s2))
    return np.roots([a, b, c])


def get_gaussian_intersection2(m1, m2, std1, std2, s1, s2, steps=10000):
    x = [1/steps*i for i in range(steps)]
    mean_a = scipy.stats.norm.pdf(x, m1, std1)[int(m1*steps)]
    mean_b = scipy.stats.norm.pdf(x, m2, std2)[int(m2*steps)]
    intersection_index = np.argwhere(np.diff(np.sign(s1/mean_a*scipy.stats.norm.pdf(x, m1, std1)-s2/mean_b*scipy.stats.norm.pdf(x, m2, std2)))).flatten()
    return [x[i] for i in intersection_index]


def get_interval_index(data, intervals):
    interval_index = [[0, 0] for _ in range(len(intervals)+1)]
    intervals      = [*intervals, *[*[data[-1]+1]]]
    
    interval_count = 0
    i = 0
    while i < len(data) and interval_count < len(intervals):
        if interval_count > 0 and data[i] >= intervals[interval_count]:
            interval_index[interval_count-1][1] = i-1
            interval_index[interval_count][0]   = i
            interval_count += 1
            #print("i1", i, data[i], interval_count, intervals[interval_count])

        if interval_count == 0:
            interval_index[interval_count][0] = i
            interval_count += 1
            #print("i2", i, data[i], interval_count, intervals[interval_count])

        i += 1

    interval_index        = interval_index[0:len(interval_index)-1]
    interval_index[-1][1] = len(data)-1
    return interval_index


def fit_poisson(x, A):
    return scipy.stats.poisson.pmf(x, A)


def get_fdr(pvalues):
    fdr_pvalues = []
    # exlude nan or None values
    #selected_pvalues  = [pvalue for pvalue in pvalues if pvalue != None]
    selected_pvalues  = [pvalue for pvalue in pvalues if pd.isna(pvalue) == False]
    corrected_pvalues = sm.stats.fdrcorrection(selected_pvalues, alpha=0.05)[1]

    i = 0
    for pvalue in pvalues:
        #if pvalue != None:
        if pd.isna(pvalue) == False:
            fdr_pvalues.append(corrected_pvalues[i])
            i += 1
        
        else:
            #fdr_pvalues.append(np.nan)
            fdr_pvalues.append(None)

    return fdr_pvalues


def get_index(selected_str):
    index = 0
    for i in range(len(selected_str)):
        match = False
        j = 0
        while j < len(alphanumericals) and match == False:
            if selected_str[i] == alphanumericals[j]:
                match = True
                index += pow(36, len(selected_str)-1-i)*j
                #print("i", i, "j", j, selected_str[i], index)

            j += 1

    return index


def get_numbers(text):
    if type(text) != str: text = str(text)
    numbers    = []
    number_str = ""

    for i in range(len(text)):
        if text[i].isdigit() == True:
            number_str += text[i]
        
        elif len(number_str) > 0:
            numbers.append(int(number_str))
            number_str = ""

    if len(number_str) > 0: numbers.append(int(number_str))
    return numbers


# targets specifies the values not to be present in multiple cohorts
def get_random_cohorts(df, cohorts, targets=None):
    if targets == None: targets = ["ID:gene id"]
    df      = df.sort_values(by=targets)
    blocks  = []
    i = 0
    while i < df.shape[0]:
        steps = 0
        j     = i
        variants        = []
        current_variant = df[targets[0]].iloc[j]
        while j < df.shape[0] and df[targets[0]].iloc[j] == current_variant:
            variants.append(j)
            current_variant = df[targets[0]].iloc[j]
            j     += 1
            steps += 1

        blocks.append(variants)
        if steps > 0: i += steps
        else:         i += 1

    bar      = IncrementalBar(set_bar("randomize variant blocks"), max=df.shape[0])
    rand_ids = np.random.permutation(len(blocks))

    if "ID:cohort" not in df.columns: cols = df.columns.insert(0, "ID:cohort")
    else:                             cols = df.columns
    rand_dict = {col: [] for col in cols}

    for i in range(cohorts):
        if i < cohorts-1:
            start_cohortsize = i*int(len(rand_ids)/cohorts)
            end_cohortsize   = (i+1)*int(len(rand_ids)/cohorts)

        else:
            start_cohortsize = i*int(len(rand_ids)/cohorts)
            end_cohortsize   = len(rand_ids)

        for j in range(start_cohortsize, end_cohortsize, 1):
            for k in range(len(blocks[rand_ids[j]])):
                for key in rand_dict:
                    if key == "ID:cohort": rand_dict[key].append(i+1)
                    else:                  rand_dict[key].append(df.iloc[blocks[rand_ids[j]][k]].loc[key])

                bar.next()
    bar.finish()

    rand_df = pd.DataFrame(rand_dict)
    return rand_df


# expanded 250113 to guarantee representation of certain items (i.e. labels) in each single cohort
# targets specifies the values not to be present in multiple cohorts
def _get_random_cohorts2(df, cohorts, cohort_exclusives=None, verbose=True):
    if cohort_exclusives == None: cohort_exclusives = ["ID:gene id"]
    df      = df.sort_values(by=cohort_exclusives)
    blocks  = []
    i = 0
    while i < df.shape[0]:
        steps = 0
        j     = i
        variants        = []
        current_variant = df[cohort_exclusives[0]].iloc[j]
        while j < df.shape[0] and df[cohort_exclusives[0]].iloc[j] == current_variant:
            variants.append(j)
            current_variant = df[cohort_exclusives[0]].iloc[j]
            j     += 1
            steps += 1

        blocks.append(variants)
        if steps > 0: i += steps
        else:         i += 1

    bar      = IncrementalBar(set_bar("randomize variant blocks"), max=df.shape[0])
    rand_ids = np.random.permutation(len(blocks))

    if "ID:cohort" not in df.columns: cols = df.columns.insert(0, "ID:cohort")
    else:                             cols = df.columns
    rand_dict = {col: [] for col in cols}

    for i in range(cohorts):
        if i < cohorts-1:
            start_cohortsize = i*int(len(rand_ids)/cohorts)
            end_cohortsize   = (i+1)*int(len(rand_ids)/cohorts)

        else:
            start_cohortsize = i*int(len(rand_ids)/cohorts)
            end_cohortsize   = len(rand_ids)

        for j in range(start_cohortsize, end_cohortsize, 1):
            for k in range(len(blocks[rand_ids[j]])):
                for key in rand_dict:
                    if key == "ID:cohort": rand_dict[key].append(i+1)
                    else:                  rand_dict[key].append(df.iloc[blocks[rand_ids[j]][k]].loc[key])

                if verbose == True: bar.next()
    if verbose == True: bar.finish()

    rand_df = pd.DataFrame(rand_dict)
    return rand_df


# expanded 250113 to guarantee representation of certain items (i.e. labels) in each single cohort (cohort_inclusives)
# cohort_exclusives specifies the values not to be present in multiple cohorts
def get_random_cohorts2(df, cohorts, cohort_exclusives=None, cohort_inclusives=None):
    if cohort_inclusives == None:
        return _get_random_cohorts2(df, cohorts, cohort_exclusives=cohort_exclusives)
    
    else:
        rand_df     = []
        df          = df.sort_values(by=cohort_inclusives)
        block_index = create_blocks(df, cohort_inclusives)
        bar = IncrementalBar("randomizing variant blocks", max=len(block_index))
        for i in range(len(block_index)):
            rand_df.append(_get_random_cohorts2(df.iloc[block_index[i]["index"]], cohorts, cohort_exclusives=cohort_exclusives, verbose=False))
            bar.next()
        bar.finish()

        return pd.concat(rand_df).sort_values(by="ID:cohort")


def get_ttest(avgs1, sems1, counts1, avgs2, sems2, counts2):
    pvalues = []
    for i in range(len(avgs1)):
        avg1  = avgs1[i]
        avg2  = avgs2[i]
        size1 = counts1[i]
        size2 = counts2[i]
        var1  = pow((sems1[i]*math.sqrt(size1)), 2)
        var2  = pow((sems2[i]*math.sqrt(size2)), 2)
        #p = -1; t = 0
        t = 0

        if size1+size2-2 > 0:
            std12 = math.sqrt(((size1-1)*var1 + (size2-1)*var2) / (size1+size2-2))

        if size1 > 1 and size2 > 1 and std12 > 0:
            t = (avg1-avg2) / (std12*math.sqrt(1/size1 + 1/size2))

        if t != 0:
            p = scipy.stats.t.sf(abs(t), df=size1+size2-2)*2
        
        pvalues.append(p)
    
    return pvalues


def get_ttest2(x : list, y : list):
    pvalue = None

    avg1  = np.mean(x)
    avg2  = np.mean(y)
    size1 = len(x)
    size2 = len(y)
    var1  = np.var(x)
    var2  = np.var(y)
    t = 0

    if size1+size2-2 > 0:
        std12 = math.sqrt(((size1-1)*var1 + (size2-1)*var2) / (size1+size2-2))

    if size1 > 1 and size2 > 1 and std12 > 0:
        t = (avg1-avg2) / (std12*math.sqrt(1/size1 + 1/size2))

    if t != 0:
        pvalue = scipy.stats.t.sf(abs(t), df=size1+size2-2)*2
       
    return pvalue


# marked (<-) added / removed on 250627
#def invert(values): # <- removed
def invert(values: list): # <- added
    return [values[len(values)-1-i] for i in range(len(values))]


def invert_sequence(seq):
    inverse_seq = ""
    bases         = ["A", "C", "G", "T"]
    inverse_bases = ["T", "G", "C", "A"]
    for i in range(len(seq)):
        j = 0
        match = False
        while j < len(bases) and match == False:
            if seq[len(seq)-1-i] == bases[j]:
                inverse_seq += inverse_bases[j]
                match = True
            j += 1

    return inverse_seq


def is_nan(values):
    if type(values) != list:
        if is_number(values) == True: return False
        else:                         return True
    
    else:
        # questionable, seems to be usable for narrow case only
        test = [value for value in values if is_number(value) == False]
        if len(test) == 0: return False
        else:              return True


def is_number(value):
    value_str = str(value).replace(".", "")
    if value_str.isdigit() == True: return True
    else:                           return False


def load_by_key(path, target, delimiter=","):
    data = pd.DataFrame()
    
    if os.path.isfile(path) == True:
        with open(path, 'r') as f:
            lines = f.readlines()

        cols           = lines[0].split(delimiter)
        target_index   = cols.index(list(target.keys())[0])

        excluded_index = [i for i, line in enumerate(lines) if i > 0 and line.split(delimiter)[target_index] not in target[list(target.keys())[0]]]
        data           = pd.read_csv(path, delimiter=",", index_col=False, skiprows=excluded_index)
        #data = pd.DataFrame({**{col: [lines[j].split(delimiter)[i].strip("\n") for j in selected_index] for i, col in enumerate(cols)}})

    else:
        print("<", path, "not found.")

    return data


def load_split_genome_(fname, path, genome, bar):
    with open(path, 'r') as f:
        next(f)
        sequence = f.read()
        sequence = sequence.replace("\n", "")

    genome[fname.split(".")[0]] = sequence
    bar.next()
    return


def load_split_genome(dir, os_sep="\\"):
    genome = {}
    fnames = os.listdir(dir)
    genome = {fname.split(".")[0]: "" for fname in fnames}
    paths  = [dir+os_sep+fname for fname in fnames]

    threads = []
    bar = IncrementalBar(set_bar("loading genome"), max=len(fnames))
    for i in range(len(fnames)):
        thread = threading.Thread(target=load_split_genome_, args=(fnames[i], paths[i], genome, bar))
        threads.append(thread)
        thread.start()

    for i, thread in enumerate(threads):
        thread.join()
    
    bar.finish()
    return genome


def load_whole_genome(path):
    genome = {} 
    with open(path, "r") as f:
        lines = f.readlines()

    bar     = IncrementalBar(set_bar("loading genome"), max=len(lines)/1000)
    chr_key = None

    for i in range(len(lines)):
        if ">" in lines[i]:
            if chr_key != None:
                if chr_key in genome: print("< warning. chromosome key", chr_key, "already exists.")
                else:                 genome[chr_key] = sequence

            chr_key  = lines[i][1::].strip("\n")
            sequence = ""

        if ">" not in lines[i]:
            sequence += lines[i].strip("\n").upper()
        
        if i % 1000 == 0: bar.next()
    bar.finish()

    if chr_key != None:
        if chr_key in genome: print("< warning. chromosome key", chr_key, "already exists.")
        else:                 genome[chr_key] = sequence
        
    return genome


def remove_trailing_digits(input):
    if type(input) == list: return [str(element).split(".")[0] for element in input]
    else:                   return str(input).split(".")[0]


def round_digitwise(number):
    number_str  = str(number)
    is_negative = False
    if number_str[0] == "-": is_negative = True

    digit_index               = None
    first_significance_index  = None
    second_significance_index = None
    i = 0
    while i < len(number_str) and (digit_index == None or second_significance_index == None):
        if first_significance_index != None and number_str[i] != "." and second_significance_index == None:
            second_significance_index = i

        if number_str[i] != "0" and number_str[i] != "-" and number_str[i] != "." and first_significance_index == None:
            first_significance_index = i

        if number_str[i] == ".":
            digit_index = i
        
        #print(i, number_str[i], number_str, first_significance_index, second_significance_index, digit_index)
        i += 1
    
    if first_significance_index != None:
        if second_significance_index != None and int(number_str[second_significance_index]) >= 5: first_significance_digit = int(number_str[first_significance_index])+1
        else:                                                                                     first_significance_digit = int(number_str[first_significance_index])

        if digit_index != None and first_significance_index < digit_index:   round_number = first_significance_digit*math.pow(10, digit_index-first_significance_index-1)
        elif digit_index != None and first_significance_index > digit_index: round_number = first_significance_digit*math.pow(10, digit_index-first_significance_index)
        else:                                                                round_number = first_significance_digit*math.pow(10, len(number_str)-first_significance_index-1)

    else:
        round_number = 0

    if number_str[0] == "-": round_number *= -1
    return round_number
    

def search_knowngene_dict_by_position(knowngene_dict, knowngene, chr, pos, stepsize=1000000):
    uc_ids = []
    # find matching uc ids
    converted_pos = int(pos/stepsize)
    if converted_pos >= len(knowngene_dict[chr]):
        print("< error occurred @search_known_gene_by_position. converted position exceeds dictionary dimension")

    else:
        search_space = [converted_pos-1, converted_pos, converted_pos+1]
        for converted_pos in search_space:
            if converted_pos < len(knowngene_dict[chr]):
                i = 0
                while i in range(len(knowngene_dict[chr][converted_pos])):
                    knowngene_index = knowngene_dict[chr][converted_pos][i]
                    #print("i", i, converted_pos, pos, knowngene_index, knowngene.iloc[knowngene_index].loc["cdsstart"], "/", knowngene.iloc[knowngene_index].loc["cdsend"])

                    if (pos >= knowngene.iloc[knowngene_index].loc["cdsstart"]
                        and pos <= knowngene.iloc[knowngene_index].loc["cdsend"]):
                        uc_ids.append(knowngene.iloc[knowngene_index].loc["uc id"])
                    
                    i += 1
    
    # return uc_ids # <- removed
    # marked (<-) added on 250427 to speed up gene detection (that requires removal of redundant entries in dictionary)
    return np.unique(uc_ids) # <- added


def set_bar(text, bar_pos=40):
    if len(text) < bar_pos:
        spaces     = bar_pos-len(text)
        for _ in range(spaces):
            text += " " 
        
    return text


def smooth_data(x, y, intervals=10):
    smoothed_x = []
    smoothed_y = []

    sorted_index = np.argsort(x)
    x            = [x[i] for i in sorted_index]
    y            = [y[i] for i in sorted_index]

    if type(intervals) == int:
        index = split_index(len(x), intervals)

    if type(intervals) == list:
        index = get_interval_index(x, intervals)

    for i in range(len(index)):
        #print("i", i, index[i][0], "/", index[i][1])
        smoothed_x.append(np.mean(x[index[i][0]:index[i][1]+1]))
        smoothed_y.append(np.mean(y[index[i][0]:index[i][1]+1]))
        #print("i", i, len((x[index[i][0]:index[i][1]+1])))

    return smoothed_x, smoothed_y


def split_index(entries, splits):
    usable_splits = min(int(splits), int(entries))

    split_size   = int(entries/usable_splits)
    thread_index = [[0 for _ in range(2)] for _ in range(usable_splits)]
    last_thread_index = -1

    for i in range(usable_splits-1):
        if last_thread_index == i*split_size: thread_index[i][0] = i*split_size+1
        else:                                 thread_index[i][0] = i*split_size
        thread_index[i][1] = (i+1)*split_size-1
        last_thread_index  = thread_index[i][1]

    if last_thread_index != -1 and last_thread_index == usable_splits-1*split_size: thread_index[usable_splits-1][0] = (usable_splits-1)*split_size+1
    else:                                                                           thread_index[usable_splits-1][0] = (usable_splits-1)*split_size
    thread_index[usable_splits-1][1] = entries-1

    total_size = 0
    for i in range(len(thread_index)):
        total_size += thread_index[i][1]-thread_index[i][0]+1
    
    if total_size != entries:
        print("< error @split_index. total size doesn't match input size:", total_size, "/", entries)

    return thread_index
