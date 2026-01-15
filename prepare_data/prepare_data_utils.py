import itertools
import json
from unittest.result import failfast
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")
sys.path.insert(0, parent_dir+"\\NMD_activity")
import threading

from extract_mutations_utils import *
from shared_utils import *


bases          = ["A", "C", "G", "T"]
chrs           = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11,
                  "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "20": 19, "21": 20, "22": 21, "X": 22, "Y": 23}
stop_codons    = ["TAA", "TAG", "TGA"]
read_through_1 = []
read_through_2 = []

for stop_codon in stop_codons:
    for base1 in bases:
        read_through_1.append(stop_codon+base1)
        for base2 in bases:
            read_through_2.append(stop_codon+base1+base2)


# append rna half-lifes: nar-00991-met-g-2009-File003.txt
def append_half_lifes(selection, df):
    half_lifes = append_df_with_mapping([selection, df], "ID:gene symbol", "Gene symbol", "Half-life [min]", "append half-life measures")
    selection.insert(selection.shape[1], 'FEATURE:RNA half-life', half_lifes)
    return selection


def append_lindeboom_predictions(selection, df):
    lindeboom_predictions = append_df_with_mapping([selection, df], "ID:variant id", "ID:variant id", "FEATURE:lindeboom prediction", "appending lindeboom predictions")
    selection.insert(selection.shape[1], 'FEATURE:lindeboom prediction', [pred if pred != "-" else None for pred in lindeboom_predictions])
    return selection


# creates n-mers of bases for pattern search
def append_features(params):
    bases          = ["A", "C", "G", "T"]
    codon_features = []
    for combination in itertools.product(bases, repeat=3):
        codon = combination[0]+combination[1]+combination[2]
        codon_features.append("FEATURE:total codon " + codon + " count cds")
        codon_features.append("FEATURE:upstream codon " + codon + " ptc count cds")
        codon_features.append("FEATURE:downstream codon " + codon + " ptc count cds")
        codon_features.append("FEATURE:total codon " + codon + " count ptc cds")
        codon_features.append("FEATURE:upstream codon " + codon + " ptc count ptc cds")
        codon_features.append("FEATURE:downstream codon " + codon + " ptc count ptc cds")

        params["motifs"]["codon " + codon]         = [codon]
        params["motif_distance"]["codon " + codon] = False
        params["motif_inframe"]["codon " + codon]  = True
        params["motif_relative"]["codon " + codon] = True

    params["features"] = [*params["features"], *codon_features]
    return params


def check_fails(fails, params):
    failed = False
    if len(fails) > 0:
        for fail in fails:
            if fail in params["error_handling"] and params["error_handling"][fail] == True:
                failed = True

    return failed


def create_knowngene_dict(knowngene, mapping_target="uc id"):
    if mapping_target not in ["uc id", "transcript id", "gene id"]:
        print("< error @create_knowngene_dict.", mapping_target, "not found.")

    knowngene_dict = {knowngene.iloc[i].loc[mapping_target].split(".")[0]: i for i in range(knowngene.shape[0])}
    return knowngene_dict


def get_chr_id(path, with_modification=True):
    chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
            "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    
    if with_modification == True:
        chr_str = path.split("chr")[1]
        chr_str = chr_str.split(".fasta")[0]

    else:
        chr_str = path

    chr_id = None
    i = 0
    while i < len(chrs) and chr_id == None:
        if chr_str == chrs[i]:
            chr_id = i

        i += 1
    
    return chr_id


def invert(seq):
    inverse_seq = ""
    bases         = ["A", "C", "G", "T", "N", "-"]
    inverse_bases = ["T", "G", "C", "A", "N", "-"]
    
    for i in range(len(seq)):
        j = 0
        match = False
        while j < len(bases) and match == False:
            if seq[len(seq)-1-i] == bases[j]:
                inverse_seq += inverse_bases[j]
                match = True
            j += 1

    return inverse_seq


def revert(seq):
    reverse_seq = ""
    for i in range(len(seq)):
        reverse_seq += seq[len(seq)-1-i]

    return reverse_seq


def fill_error_report(error_report, selection, fails, strand, it):
    indel_tag = "1"
    if "ID:HGVSc" in selection.columns and "_" in selection.iloc[it].loc["ID:HGVSc"] and ("del" in selection.iloc[it].loc["ID:HGVSc"] or "ins" in selection.iloc[it].loc["ID:HGVSc"]):
        indel_tag =">1"

    if len(fails) == 0:
        error_report["total"]["no_error"].append(it)
        error_report[strand+"total"]["no_error"].append(it) # added on 250616

        if "ID:HGVSc" in selection.columns and "dup" in selection.iloc[it].loc["ID:HGVSc"]:
            error_report[strand+"dup"+indel_tag]["no_error"].append(it)

        elif "ID:HGVSc" in selection.columns and "del" in selection.iloc[it].loc["ID:HGVSc"] and "delins" not in selection.iloc[it].loc["ID:HGVSc"]:
            error_report[strand+"del"+indel_tag]["no_error"].append(it)

        elif "ID:HGVSc" in selection.columns and "ins" in selection.iloc[it].loc["ID:HGVSc"] and "delins" not in selection.iloc[it].loc["ID:HGVSc"]:
            error_report[strand+"ins"+indel_tag]["no_error"].append(it)

        elif "ID:HGVSc" in selection.columns and "delins" in selection.iloc[it].loc["ID:HGVSc"]:
            error_report[strand+"delins"+indel_tag]["no_error"].append(it)

        elif "ID:HGVSc" in selection.columns:
            error_report[strand+"nonsense"]["no_error"].append(it)

    else:
        for error in fails:
            if "ID:HGVSc" in selection.columns and "dup" in selection.iloc[it].loc["ID:HGVSc"]:
                error_report[strand+"dup"+indel_tag][error].append(it)
                # error_report["total"][error].append(it) # removed on 250616

            elif "ID:HGVSc" in selection.columns and "del" in selection.iloc[it].loc["ID:HGVSc"] and "delins" not in selection.iloc[it].loc["ID:HGVSc"]:
                error_report[strand+"del"+indel_tag][error].append(it)
                # error_report["total"][error].append(it) # removed on 250616

            elif "ID:HGVSc" in selection.columns and "ins" in selection.iloc[it].loc["ID:HGVSc"] and "delins" not in selection.iloc[it].loc["ID:HGVSc"]:
                error_report[strand+"ins"+indel_tag][error].append(it)
                # error_report["total"][error].append(it) # removed on 250616

            elif "ID:HGVSc" in selection.columns and "delins" in selection.iloc[it].loc["ID:HGVSc"]:
                error_report[strand+"delins"+indel_tag][error].append(it)
                # error_report["total"][error].append(it) # removed on 250616

            elif "ID:HGVSc" in selection.columns:
                error_report[strand+"nonsense"][error].append(it)
                # error_report["total"][error].append(it) # removed on 250616

            error_report["total"][error].append(it) # added on 250616
            error_report[strand+"total"][error].append(it) # added on 250616

    return error_report


def get_info(df, variant_id, match_index, ptc_index, mutation_index, wt_base, mutated_base, cds_size, genome, params, ptc_mode="absolute"):
    info  = {}
    fails = []

    for feature in params["features"]:
        if "FEATURE" in feature: info[feature] = info.get(feature, 0)

    extracted_sequence = []
    all_cds = ""; all_exons = ""; utr5 = ""; utr3 = ""
    cds_start          = False
    ptc_index_found    = False
    ptc_in_exon        = False
    ptc_cds_index      = None
    ptc_exon_index     = None
    upstream_cds_size  = None
    upstream_exon_size = None

    transcript         = df.iloc[match_index]
    exonsstart         = [int(i) for i in transcript.loc["exonsstart"].split(",") if len(i) > 0] # condition required because last letter is a comma
    exonsend           = [int(i) for i in transcript.loc["exonsend"].split(",") if len(i) > 0] # condition required because last letter is a comma

    # check whether chromosome is present
    chr = transcript.loc["chr"]
    if ptc_mode == "absolute" and "_" in chr: fails.append("chromosome_incompatible")
    if chr not in genome:                     fails.append("chromosome_not_found"); print(variant_id, transcript.loc["chr"], chr)
    show = False
    #if transcript.loc["strand"] == "-": show = True
    if show == True: print("variant_id", variant_id, "ptc_index", ptc_index, "mutation_index", mutation_index, "wt_base", wt_base)
    if show == True: print("transcript", transcript)
    if transcript.loc["strand"] == "+": i = 0
    if transcript.loc["strand"] == "-": i = transcript.loc["exons"]-1
    
    while ((transcript.loc["strand"] == "+" and i < transcript.loc["exons"]) or (transcript.loc["strand"] == "-" and i >= 0)):
        info["FEATURE:exons"]           += 1
        info["FEATURE:total exon size"] += exonsend[i]-exonsstart[i]
                        
        # register position of the last exon junction complex (penultimate exon)
        if transcript.loc["strand"] == "+" and i == transcript.loc["exons"]-2: info["FEATURE:dist. from last EJC"] = info["FEATURE:total exon size"]
        if transcript.loc["strand"] == "-" and i == 1:                         info["FEATURE:dist. from last EJC"] = info["FEATURE:total exon size"]

        if transcript.loc["strand"] == "+":   all_exons += genome[chr][exonsstart[i]:exonsend[i]]
        elif transcript.loc["strand"] == "-":
            all_exons += invert(genome[chr][exonsstart[i]:exonsend[i]])

        if show == True: print("i1", i, "exonsstart", exonsstart[i], "exonsend", exonsend[i], "cdsend", transcript.loc["cdsend"])

    	# calculations for exons with no CDS (UTR-only)
        if (transcript.loc["cdsstart"] == transcript.loc["cdsend"] # as defined in the knownGene database
            or (transcript.loc["cdsstart"] != transcript.loc["cdsend"]
            and (exonsend[i] < transcript.loc["cdsstart"] or exonsstart[i] > transcript.loc["cdsend"]))):
            if cds_start == False:
                if transcript.loc["strand"] == "+":   info["FEATURE:5'utr-size"] += exonsend[i]-exonsstart[i]
                elif transcript.loc["strand"] == "-": info["FEATURE:5'utr-size"] += exonsend[i]-exonsstart[i]
                if transcript.loc["strand"] == "+":   utr5 += genome[chr][exonsstart[i]:exonsend[i]]
                elif transcript.loc["strand"] == "-": utr5 += invert(genome[chr][exonsstart[i]:exonsend[i]])

            elif cds_start == True:
                if transcript.loc["strand"] == "+":   info["FEATURE:3'utr-size"] += exonsend[i]-exonsstart[i]
                elif transcript.loc["strand"] == "-": info["FEATURE:3'utr-size"] += exonsend[i]-exonsstart[i]
                if transcript.loc["strand"] == "+":   utr3 += genome[chr][exonsstart[i]:exonsend[i]]
                elif transcript.loc["strand"] == "-": utr3 += invert(genome[chr][exonsstart[i]:exonsend[i]])
            
            if show == True: print("i2", i, "5utr-size", info["FEATURE:5'utr-size"], "/", len(utr5), "3utr-size", info["FEATURE:3'utr-size"], len(utr3))
            if ptc_index >= exonsstart[i] and ptc_index <= exonsend[i]: fails.append("ptc_in_utr")

        # calculations for exons with CDS
        else:
            # marked (<-) added on 250425 to account for cds exons only
            info["FEATURE:cds exons"] += 1 # <-
            cds_start = True

            # check if exon contains non-cds
            if exonsstart[i] < transcript.loc["cdsstart"]: exon_cds_start = transcript.loc["cdsstart"]
            else:                                          exon_cds_start = exonsstart[i]
            if exonsend[i] > transcript.loc["cdsend"]:     exon_cds_end   = transcript.loc["cdsend"]
            else:                                          exon_cds_end   = exonsend[i]
            if show == True: print("i3", i, "exon_cds_start", exon_cds_start, "exon_cds_end", exon_cds_end)

            # register position of the second next exon junction complex from the ptc
            info["FEATURE:total cds size"] += exon_cds_end-exon_cds_start
            info["FEATURE:cds"]            += 1
        
            # calculate absolute ptc index from relative ptc index if datatype=custom or tcga
            if ptc_mode == "relative" and ptc_index < info["FEATURE:total cds size"] and ptc_index_found == False:
                ptc_index_found  = True
                test_index       = ptc_index
                ptc_cds_exon_pos = ptc_index-(info["FEATURE:total cds size"]-(exon_cds_end-exon_cds_start))
                if transcript.loc["strand"] == "+":   ptc_index = exon_cds_start + ptc_cds_exon_pos                
                elif transcript.loc["strand"] == "-": ptc_index = exon_cds_end-1 - ptc_cds_exon_pos
                info["FEATURE:abs. ptc index"] = ptc_index
                if show == True: print("i4", i, "ptc_index", ptc_index, "ptc_index_before", test_index, "ptc_cds_pos", ptc_cds_exon_pos)
                
            if transcript.loc["strand"] == "+":
                all_cds                    += genome[chr][exon_cds_start:exon_cds_end]
                info["FEATURE:5'utr-size"] += exon_cds_start-exonsstart[i]
                utr5                       += genome[chr][exonsstart[i]:exon_cds_start] 
                info["FEATURE:3'utr-size"] += exonsend[i]-exon_cds_end
                utr3                       += genome[chr][exon_cds_end:exonsend[i]]

                if show == True: print("i5", i, "utr5+", len(utr5), "/", info["FEATURE:5'utr-size"], "utr3", len(utr3), "/", info["FEATURE:3'utr-size"], "cds length", len(all_cds))

            elif transcript.loc["strand"] == "-":
                all_cds                    += invert(genome[chr][exon_cds_start:exon_cds_end])
                info["FEATURE:3'utr-size"] += exon_cds_start-exonsstart[i]
                if exon_cds_start-exonsstart[i] > 0: utr3 += invert(genome[chr][exonsstart[i]:exon_cds_start])
                info["FEATURE:5'utr-size"] += exonsend[i]-exon_cds_end
                if exonsend[i]-exon_cds_end > 0:     utr5 += invert(genome[chr][exon_cds_end:exonsend[i]])
                if show == True: print("i6", i, "utr5-size", info["FEATURE:5'utr-size"], "/", len(utr5), "utr3-size", info["FEATURE:3'utr-size"], "/", len(utr3), "cds length", len(all_cds))

            if ptc_index >= exon_cds_start and ptc_index < exonsend[i] and ptc_index < exon_cds_end:
                if show == True: print("ptc_in_exon", ptc_in_exon)
                ptc_in_exon        = True
                upstream_cds_size  = info["FEATURE:total cds size"] - (exon_cds_end-exon_cds_start)
                upstream_exon_size = info["FEATURE:total exon size"] - (exonsend[i]-exonsstart[i])

                if transcript.loc["strand"] == "+":
                    info["FEATURE:ptc upstream distance"]   = ptc_index-exonsstart[i]
                    info["FEATURE:ptc downstream distance"] = exonsend[i]-1-ptc_index
                    ptc_cds_index                           = (info["FEATURE:total cds size"] + ptc_index-max(exonsstart[i], transcript.loc["cdsstart"])
                                                               - (exon_cds_end-exon_cds_start))
                    ptc_exon_index                          = info["FEATURE:total exon size"] - (exonsend[i]-exonsstart[i]) + (ptc_index-exonsstart[i])
                    if show == True: print("i7", i, "ptc-up", info["FEATURE:ptc upstream distance"], "ptc-down", info["FEATURE:ptc downstream distance"])
                    try:
                        extracted_sequence.append({"info": "ptc cds", "sequence": genome[chr][exon_cds_start:exon_cds_end]})
                        extracted_sequence.append({"info": "ptc", "sequence": genome[chr][exonsstart[i]:exonsend[i]]})
                        extracted_sequence.append({"info": "5'ptc", "sequence": genome[chr][exonsstart[i]:ptc_index]})
                        extracted_sequence.append({"info": "3'ptc", "sequence": genome[chr][ptc_index:exonsend[i]]})

                    except:
                        print("error", chr, exon_cds_start, exon_cds_end, ptc_index)                        
                        exit()

                elif transcript.loc["strand"] == "-":
                    info["FEATURE:ptc upstream distance"]   = exonsend[i]-1-ptc_index
                    info["FEATURE:ptc downstream distance"] = ptc_index-exonsstart[i]
                    ptc_cds_index                           = (info["FEATURE:total cds size"] + min(exonsend[i], transcript.loc["cdsend"])-1-ptc_index
                                                                - (exon_cds_end-exon_cds_start))
                    ptc_exon_index                          = info["FEATURE:total exon size"] - (exonsend[i]-exonsstart[i]) + (exonsend[i]-1-ptc_index)
                    
                    if show == True: print("i8", i, "ptc-up", info["FEATURE:ptc upstream distance"], "ptc-down", info["FEATURE:ptc downstream distance"])
                    extracted_sequence.append({"info": "ptc cds", "sequence": invert(genome[chr][exon_cds_start:exon_cds_end])})
                    extracted_sequence.append({"info": "ptc", "sequence": invert(genome[chr][exonsstart[i]:exonsend[i]])})
                    extracted_sequence.append({"info": "5'ptc", "sequence": invert(genome[chr][ptc_index:exonsend[i]])})
                    extracted_sequence.append({"info": "3'ptc", "sequence": invert(genome[chr][exonsstart[i]:ptc_index])})

                if ptc_cds_index != None:  info["FEATURE:ptc cds position"]  = ptc_cds_index
                else:                      fails.append("ptc_cds_index_not_found")
                if ptc_exon_index != None: info["FEATURE:ptc exon position"] = ptc_exon_index
                else:                      fails.append("ptc_exon_index_not_found")

            # marked (<-) added on 250425 to account for cds exons only
            if ptc_in_exon == True: # <-
                info["FEATURE:downstream cds exons"] += 1 # <-
            
        if ptc_in_exon == True:
            info["FEATURE:downstream exons"] += 1
            if show == True: print("i9", i)
        
        if transcript.loc["strand"] == "+": i += 1
        if transcript.loc["strand"] == "-": i -= 1


    # check if PTC was not found
    if ptc_in_exon == False: fails.append("ptc_not_in_exon")

    # check whether mutation is in cds
    if ptc_index < transcript.loc["cdsstart"] or ptc_index > transcript.loc["cdsend"]: fails.append("ptc_not_in_cds")
    
    # calculate the number of upstream exons
    info["FEATURE:upstream exons"] = info["FEATURE:exons"]-info["FEATURE:downstream exons"]
    # marked (<-) added on 250425 to account for cds exons only
    info["FEATURE:upstream cds exons"] = info["FEATURE:cds exons"]-info["FEATURE:downstream cds exons"] # <-

    # marked (<-) added on 250425 to integrate additional features motivated by Kim et al.
    if info["FEATURE:downstream exons"] == 1: # <-
        info["FEATURE:last exon"] = 1 # <-

    if info["FEATURE:downstream cds exons"] == 1: # <-
        info["FEATURE:last cds exon"] = 1 # <- 

    if info["FEATURE:downstream exons"] == 1 or (info["FEATURE:downstream exons"] == 2 and info["FEATURE:ptc downstream distance"] <= 50): # <-
        info["FEATURE:50 nt to last EJC"] = 1 # <-

    if info["FEATURE:downstream cds exons"] == 1 or (info["FEATURE:downstream cds exons"] == 2 and info["FEATURE:ptc downstream distance"] <= 50): # <-
        info["FEATURE:50 nt to last cds EJC"] = 1 # <-

    info["FEATURE:ptc exon size"] = info["FEATURE:ptc upstream distance"]+info["FEATURE:ptc downstream distance"] # <-

    # calculate the EJC density (complete and downstream)
    if info["FEATURE:total exon size"] > 0:
        info["FEATURE:EJC density"]            = (info["FEATURE:exons"]-1) / info["FEATURE:total exon size"]
        info["FEATURE:downstream EJC density"] = (info["FEATURE:downstream exons"]-1) / info["FEATURE:total exon size"]

    else:
        fails.append("no_exons")

    # marked (<-) added on 250425 to account for cds exons only
    if info["FEATURE:total cds size"] > 0: # <-
        info["FEATURE:cds EJC density"]            = (info["FEATURE:cds exons"]-1) / info["FEATURE:total cds size"] # <-
        info["FEATURE:downstream cds EJC density"] = (info["FEATURE:downstream cds exons"]-1) / info["FEATURE:total cds size"] # <-

    else: # <-
        fails.append("no_cds") # <-

    # must be done here as relative ptc index is converted to absolute index during iteration
    if ptc_exon_index != None and info["FEATURE:dist. from last EJC"] != 0:
        info["FEATURE:dist. from last EJC"] = str(info["FEATURE:dist. from last EJC"]-ptc_exon_index-1)

    # add skip tag if item is zero (could not be determined)
    elif info["FEATURE:dist. from last EJC"] == 0:
        info["FEATURE:dist. from last EJC"] = None

    # write sequences to dictionaries and append sequence report
    if ptc_cds_index != None: info["FEATURE:ptc-wt stop codon distance"] = info["FEATURE:total cds size"]-ptc_cds_index
    extracted_sequence.insert(0, {"info": "5'utr", "sequence": utr5})
    extracted_sequence.append({"info": "3'utr", "sequence": utr3})
    extracted_sequence.append({"info": "all cds", "sequence": all_cds})
    extracted_sequence.append({"info": "all exons", "sequence": all_exons})

    # check for correct start and stop codons
    if transcript.loc["cdsstart"] != transcript.loc["cdsend"]: # exclude transcripts without CDS
        if all_cds[0:3] != "ATG":                                                                                            fails.append("no_start")
        if all_cds[len(all_cds)-3::] != "TAG" and all_cds[len(all_cds)-3::] != "TGA" and all_cds[len(all_cds)-3::] != "TAA": fails.append("no_stop")

    for section in extracted_sequence:
        if section["info"] == "5'utr":
            if len(section["sequence"]) > 0: info["FEATURE:5'utr"] = section["sequence"]
            else:                            info["FEATURE:5'utr"] = "X"

        elif section["info"] == "ptc cds":
            if len(section["sequence"]) > 0: info["FEATURE:ptc cds"] = section["sequence"]
            else:                            info["FEATURE:ptc cds"] = "X"

        elif section["info"] == "ptc":
            if len(section["sequence"]) > 0: info["FEATURE:ptc"] = section["sequence"]
            else:                            info["FEATURE:ptc"] = "X"
            
        elif section["info"] == "5'ptc":
            if len(section["sequence"]) > 0: info["FEATURE:5'ptc"] = section["sequence"]
            else:                            info["FEATURE:5'ptc"] = "X"
            if len(section["sequence"]) > 0: info["FEATURE:5'ejc"] = revert(section["sequence"])
            else:                            info["FEATURE:5'ejc"] = "X"

        elif section["info"] == "3'ptc":
            if len(section["sequence"]) > 0: info["FEATURE:3'ptc"] = section["sequence"]
            else:                            info["FEATURE:3'ptc"] = "X"
            if len(section["sequence"]) > 0: info["FEATURE:3'ejc"] = revert(section["sequence"])
            else:                            info["FEATURE:3'ejc"] = "X"

        elif section["info"] == "3'utr":
            if len(section["sequence"]) > 0: info["FEATURE:3'utr"] = section["sequence"]
            else:                            info["FEATURE:3'utr"] = "X"

        elif section["info"] == "all cds":
            if len(section["sequence"]) > 0: info["FEATURE:all cds"] = section["sequence"]
            else:                            info["FEATURE:all cds"] = "X"

        elif section["info"] == "all exons":
            if len(section["sequence"]) > 0: info["FEATURE:all exons"] = section["sequence"]
            else:                            info["FEATURE:all exons"] = "X"

    # search for patterns in the cds
    if type(info["FEATURE:all cds"]) == str: 
        info, fails = get_motif(info, fails, params, info["FEATURE:all cds"], "cds", framestart=0, selected_motifs=["start"],
                                relative_ptc_index=ptc_cds_index, variant_id=variant_id)
    
    # search for patterns in the ptc cds
    if type(info["FEATURE:ptc cds"]) == str and upstream_cds_size != None:
        framestart = (ptc_cds_index-upstream_cds_size) - int((ptc_cds_index-upstream_cds_size)/3)*3
        if framestart != 0: framestart = 3-framestart
        info, fails = get_motif(info, fails, params, info["FEATURE:ptc cds"], "ptc cds", framestart=framestart, relative_ptc_index=ptc_cds_index-upstream_cds_size,
                                selected_motifs=["start"], variant_id=variant_id)
    
    # search for codon patterns in the cds
    #print([motif for motif in params["motifs"] if "codon" in motif])
    if type(info["FEATURE:all cds"]) == str: 
        info, fails = get_motif(info, fails, params, info["FEATURE:all cds"], "cds", framestart=0, selected_motifs=[motif for motif in params["motifs"] if "codon" in motif],
                                relative_ptc_index=ptc_cds_index, variant_id=variant_id)
    
    # search for codon patterns in the ptc cds
    if type(info["FEATURE:ptc cds"]) == str and upstream_cds_size != None:
        framestart = (ptc_cds_index-upstream_cds_size) - int((ptc_cds_index-upstream_cds_size)/3)*3
        if framestart != 0: framestart = 3-framestart
        info, fails = get_motif(info, fails, params, info["FEATURE:ptc cds"], "ptc cds", framestart=framestart, relative_ptc_index=ptc_cds_index-upstream_cds_size,
                                selected_motifs=[motif for motif in params["motifs"] if "codon" in motif], variant_id=variant_id)
        
    # search for patterns in the exons
    if type(info["FEATURE:all exons"]) == str: 
        info, fails = get_motif(info, fails, params, info["FEATURE:all exons"], "exons", framestart=0, relative_ptc_index=ptc_exon_index,
                                selected_motifs=["GC"], variant_id=variant_id)
    
    # search for patterns in the ptc exon
    if type(info["FEATURE:ptc"]) == str and upstream_exon_size != None:
        info, fails = get_motif(info, fails, params, info["FEATURE:ptc"], "ptc exon", framestart=0, relative_ptc_index=ptc_exon_index-upstream_exon_size,
                                selected_motifs=["GC"], variant_id=variant_id)
 

    if ptc_mode == "absolute": info["FEATURE:abs. ptc index"] = ptc_index

    if ptc_mode == "relative" and transcript.loc["strand"] == "-":
        wt_base = invert(wt_base)
        if mutated_base != None: mutated_base = invert(mutated_base) # <- added on 260109 to calculate read-through stop codons
       
    # <- added on 260109 to include read-through info
    info, fails = get_read_through_motif(info, fails, mutation_index, mutated_base, ptc_mode, variant_id)

    # added on 241126
    if len(wt_base) > 1: wt_base = wt_base[0]

    if (ptc_mode == "absolute" and (mutation_index >= len(genome[chr]) or
        (mutation_index < len(genome[chr]) and wt_base != "-" and genome[chr][mutation_index] != wt_base))):
        fails.append("mutated_wt_base_not_found")

    if (ptc_mode == "relative" and (mutation_index >= len(info["FEATURE:all cds"])
        or (mutation_index < len(info["FEATURE:all cds"]) and wt_base != "-" and info["FEATURE:all cds"][mutation_index] != wt_base))):
        fails.append("mutated_wt_base_not_found")

    if cds_size != None and cds_size != len(info["FEATURE:all cds"])-3:
        fails.append("inconsistent_cds_size")
        print(variant_id, cds_size, "/", len(info["FEATURE:all cds"])-3)
        print(fails)

    if info["FEATURE:ptc cds position"] == len(info["FEATURE:all cds"])-3:
        fails.append("ptc_equals_stop_codon")

    if show == True:
        print(info); print(fails)

    return info, fails


def get_motif(info, fails, params, seq, seq_name, framestart=0, relative_ptc_index=None, selected_motifs=None, variant_id=None):
    if selected_motifs == None:
        motifs         = params["motifs"]
        motif_distance = params["motif_distance"]
        motif_inframe  = params["motif_inframe"]
        motif_relative = params["motif_relative"]

    else:
        motifs         = {motif: params["motifs"][motif] for motif in params["motifs"] if motif in selected_motifs}
        motif_distance = {motif: params["motif_distance"][motif] for motif in params["motifs"] if motif in selected_motifs}
        motif_inframe  = {motif: params["motif_inframe"][motif] for motif in params["motifs"] if motif in selected_motifs}
        motif_relative = {motif: params["motif_relative"][motif] for motif in params["motifs"] if motif in selected_motifs}

    keycount = 0
    for key in motifs:
        #print(key, seq_name, len(seq), relative_ptc_index)
        motifcount = 0
        upstream_entry_made = False; downstream_entry_made = False
        
        for motif in motifs[key]:
            if motif_inframe[key] == True:    step = 3; current_framestart = int(framestart)
            elif motif_inframe[key] == False: step = 1; current_framestart = 0
            
            # count total occurrences of patterns (might deviate from upstream + downstream because here, ptc region is taken into account)
            index = [i for i in range(current_framestart, len(seq), step) if seq.startswith(motif, i)]
            if motif_relative[key] == False: info["FEATURE:total " + key + " count " + seq_name] += len(index) #print(seq_name, "total", index, "step", step, "current_framestart", current_framestart)
            else:                            info["FEATURE:total " + key + " count " + seq_name] += len(index) / len(seq)

            if relative_ptc_index != None:
                # count upstream instances of motif
                index = [i for i in range(current_framestart, int(relative_ptc_index)-len(motif)+1, step) if seq[i:i+len(motif)] == motif]
                if motif_relative[key] == False:
                    info["FEATURE:upstream " + key + " ptc count " + seq_name] += len(index)
                    #print(seq_name, "upstream", index, "current_framestart", current_framestart)

                elif int(relative_ptc_index)-len(motif) > 0:
                    if ((int(relative_ptc_index)-len(motif)+1)-current_framestart) >= params["pattern_cutoff"]: info["FEATURE:upstream " + key + " ptc count " + seq_name] += len(index) / ((int(relative_ptc_index)-len(motif)+1)-current_framestart)#len(seq)
                    else:                                                                                       info["FEATURE:upstream " + key + " ptc count " + seq_name]  = 0.5
                    #print(seq_name, "upstream", len(index), ((int(relative_ptc_index)-len(motif))+1-current_framestart))
                
                if variant_id == "":
                    print("motif", motif, motifs[key], "stats", "FEATURE:upstream " + key + " ptc count " + seq_name,
                        info["FEATURE:upstream " + key + " ptc count " + seq_name], "seq", len(seq), "index", relative_ptc_index, int(relative_ptc_index)-len(motif))

                # calculate upstream distance from ptc to the first instance of the motif
                if motif_distance[key] == True:
                    if len(index) > 0 and (motifcount == 0 or (motifcount > 0 and int(info["FEATURE:first upstream " + key + " ptc distance " + seq_name]) > index[len(index)-1])):
                        #info["FEATURE:first upstream " + key + " ptc distance " + seq_name] = str(int(relative_ptc_index)-index[len(index)-1]) #convert to string for later easy replacement by "-" skip tag
                        info["FEATURE:first upstream " + key + " ptc distance " + seq_name] = int(relative_ptc_index)-index[len(index)-1]
                        upstream_entry_made = True
                        if seq_name == "cds" and key == "start" and index[0] != 0: fails.append("motif_start_not_found")

                # count downstream instances of motif
                if motif_inframe[key] == True:    start_index = (int(relative_ptc_index/step)+1)*step + current_framestart
                elif motif_inframe[key] == False: start_index = int(relative_ptc_index)+1

                index = [i for i in range(start_index, len(seq)-len(motif)+1, step) if seq[i:i+len(motif)] == motif]
                if motif_relative[key] == False:
                    info["FEATURE:downstream " + key + " ptc count " + seq_name] += len(index)

                else:
                    if ((len(seq)-len(motif)+1)-start_index) >= params["pattern_cutoff"]: info["FEATURE:downstream " + key + " ptc count " + seq_name] += len(index) / ((len(seq)-len(motif)+1)-start_index)
                    else:                                                                 info["FEATURE:downstream " + key + " ptc count " + seq_name]  = 0.5 
                
                if motif_distance[key] == True:
                    # calculate downstream distance from ptc to the first instance of the motif
                    if len(index) > 0 and (motifcount == 0 or (motifcount > 0 and int(info["FEATURE:first downstream " + key + " ptc distance " + seq_name]) > index[0])):
                        #info["FEATURE:first downstream " + key + " ptc distance " + seq_name] = str(index[0]-int(relative_ptc_index)) #convert to string for later easy replacement by "-" skip tag
                        info["FEATURE:first downstream " + key + " ptc distance " + seq_name] = index[0]-int(relative_ptc_index)
                        downstream_entry_made = True

            motifcount += 1
        
        if motif_distance[key] == True:
            if upstream_entry_made == False:   info["FEATURE:first upstream " + key + " ptc distance " + seq_name]   = None
            if downstream_entry_made == False: info["FEATURE:first downstream " + key + " ptc distance " + seq_name] = None

        keycount += 1

    return info, fails


# gets first position from mutation info on cDNA level, e.g. c.2333_2356delAGCCATGGGGCGAGATGAAGGAGC -> 2333
def get_position_from_mutation(mutation_info):
    position = ""
    i = 2
    while mutation_info[i].isnumeric() == True and i < len(mutation_info):
        position += mutation_info[i]
        i += 1

    if len(position) == 0: return None
    else:                  return int(position)


# <- added on 260109 to include read-through info
def get_read_through_motif(info, fails, mutation_index, mutated_base, ptc_mode, variant_id):
    if ptc_mode == "absolute" and mutated_base == None:
        info["FEATURE:read-through+1"] = None
        info["FEATURE:read-through+2"] = None
        
    else:
        # replace mutated base
        cds = "".join(base if i != mutation_index else mutated_base for i, base in enumerate(info["FEATURE:all cds"]))

        if cds[info["FEATURE:ptc cds position"]-2:info["FEATURE:ptc cds position"]+1] not in stop_codons:
            print("< error @get_read_through_motif. stop codon not detected for variant id:", variant_id)
            print(" context:", cds[info["FEATURE:ptc cds position"]-3:info["FEATURE:ptc cds position"]+3], "/",
                               info["FEATURE:all cds"][info["FEATURE:ptc cds position"]-3:info["FEATURE:ptc cds position"]+3])
            info["FEATURE:read-through+1"] = None
            info["FEATURE:read-through+2"] = None
            fails.append("mutated_stop_codon_not_found")
            
        else:
            info["FEATURE:read-through+1"] = read_through_1.index(cds[info["FEATURE:ptc cds position"]-2:info["FEATURE:ptc cds position"]+2])
            info["FEATURE:read-through+2"] = read_through_2.index(cds[info["FEATURE:ptc cds position"]-2:info["FEATURE:ptc cds position"]+3])
    
    return info, fails


# <- added on 251027 averaging of scores not supported
def get_control_selection(df, params):
    selection_df = {"ID:sample id": [], "ID:strand": [], "ID:chromosome": [], "ID:gene symbol": [],"ID:transcript id": [],
                    "ID:variant id": [], "ID:HGVSc": [], "ID:HGVSp": [], "ID:variant classification": [],
                    "ID:abs. mutation index": [], "ID:cds size": [], "ID:mutation index": [], "ID:wt base": [], "ID:ptc index": [],
                    "ID:wt expression": [], "ID:wt expression, direct avg": [], "ID:allelic fraction": [], "ID:base mean": [],
                    "LABEL:NMD score": [], "LABEL:NMD score padj": []}
    
    bar = IncrementalBar(set_bar("selecting control data"), max=df.shape[0])
    for i in range(df.shape[0]):
        selection_df["ID:sample id"].append(df.iloc[i].loc["SAMPLE_ID"])
        selection_df["ID:chromosome"].append(df.iloc[i].loc["Chromosome"])
        selection_df["ID:strand"].append(df.iloc[i].loc["Strand"])
        selection_df["ID:gene symbol"].append(df.iloc[i].loc["Gene_Hugo"])
        selection_df["ID:transcript id"].append(df.iloc[i].loc["Ensembl_transcript_id"].split(".")[0])
        selection_df["ID:variant id"].append(df.iloc[i].loc["Ensembl_transcript_id"].split(".")[0]+"_"+str(df.iloc[i].loc["Start"])+"_"+str(df.iloc[i].loc["End"]))
        selection_df["ID:variant classification"].append(df.iloc[i].loc["Type_1"]) # <- added on 251027

        selection_df["ID:abs. mutation index"].append(df.iloc[i].loc["Start"]-1)

        selection_df["ID:HGVSc"].append(df.iloc[i].loc["Change_cDNA"])
        selection_df["ID:HGVSp"].append(df.iloc[i].loc["Change_Protein"])

        selection_df["ID:mutation index"].append(get_numbers(df.iloc[i].loc["Change_cDNA"])[0]-1)
        selection_df["ID:wt base"].append(df.iloc[i].loc["Wild_Type"])

        # required as placeholder
        selection_df["ID:cds size"].append(None)
        selection_df["ID:ptc index"].append(df.iloc[i].loc["ptc index"])
        selection_df["ID:base mean"].append(df.iloc[i].loc["baseMean"])
        selection_df["ID:wt expression"].append(df.iloc[i].loc["wt expression"])
        selection_df["ID:wt expression, direct avg"].append(df.iloc[i].loc["wt expression, direct avg"])
        selection_df["ID:allelic fraction"].append(df.iloc[i].loc["Allelic_Fraction_Tumor"])
        selection_df["LABEL:NMD score"].append(df.iloc[i].loc["log2FoldChange"])
        selection_df["LABEL:NMD score padj"].append(df.iloc[i].loc["padj"])
        bar.next()
    bar.finish()
    
    selection  = pd.DataFrame.from_dict(selection_df)
    init_shape = selection.shape[0]
    selection  = selection[~selection["ID:ptc index"].isna()]
    print("< removal of missing ptc index reduced values from", init_shape, "to", selection.shape[0])
    # index needs to be changed to int, else causing downstream error @get_info
    selection = selection.astype({"ID:mutation index": "int32", "ID:ptc index": "int32"})
    return selection


# averaging of scores not supported
def get_mskcc_selection(df, params):
    bar = IncrementalBar(set_bar("selecting MSKCC data"), max=df.shape[0])
    selection_df = {"ID:entrez gene id": [], "ID:gene symbol": [], "ID:patient id": [], "ID:sample id": [], "ID:transcript id": [],
                    "ID:variant id": [], "ID:chromosome": [], "ID:HGVSc": [], "ID:HGVSp": [],
                    "ID:variant classification": [], # <- added on 251027
                    "ID:abs. mutation index": [], "ID:cds size": [], "ID:mutation index": [], "ID:wt base": [], "ID:ptc index": [],
                    "ID:cds size": []}
     
    removed_ptcs = []

    # add row (last row is not considered)
    df.loc[len(df.index)] = [None for _ in df.columns]

    for i in range(df.shape[0]):
        if pd.isna(df.iloc[i].loc["HGVSp"]) == False:
            if get_position_from_mutation(df.iloc[i].loc["HGVSc"].split(":")[1]) != None:
                selection_df["ID:entrez gene id"].append(df.iloc[i].loc["Entrez_Gene_Id"])
                selection_df["ID:gene symbol"].append(df.iloc[i].loc["Hugo_Symbol"])
                selection_df["ID:patient id"].append(df.iloc[i].loc["Tumor_Sample_Barcode"].split("-")[0]+"-"+df.iloc[i].loc["Tumor_Sample_Barcode"].split("-")[1])
                selection_df["ID:sample id"].append(df.iloc[i].loc["Tumor_Sample_Barcode"])
                selection_df["ID:transcript id"].append(df.iloc[i].loc["Transcript_ID"])
                selection_df["ID:variant id"].append(df.iloc[i].loc["Transcript_ID"]+"_"+str(df.iloc[i].loc["Start_Position"])+"_"+str(df.iloc[i].loc["End_Position"]))
                selection_df["ID:variant classification"].append(df.iloc[i].loc["Variant_Classification"]) # <- added on 251027

                selection_df["ID:abs. mutation index"].append(df.iloc[i].loc["Start_Position"]-1)
                selection_df["ID:chromosome"].append("chr"+df.iloc[i].loc["Chromosome"])

                selection_df["ID:HGVSc"].append(df.iloc[i].loc["HGVSc"])
                selection_df["ID:HGVSp"].append(df.iloc[i].loc["HGVSp"])

                selection_df["ID:mutation index"].append(get_numbers(df.iloc[i].loc["HGVSc"].split(":")[1])[0]-1)
                selection_df["ID:wt base"].append(df.iloc[i].loc["Reference_Allele"])

                # required as placeholder
                selection_df["ID:cds size"].append(None)

                # changed on 241112
                ptc = get_ptc(df.iloc[i].loc["HGVSp"])

                # <- added / removed on 251027 to exclude HGVSp with '?'
                # if ptc != None: # <- removed
                if ptc != None and get_frameshift(df.iloc[i].loc["HGVSp"]) != None: # <- added
                    # marked (<-) added / replaced on 250425 to calculate with shifted relative PTC index (+2) to avoid assignment to wrong exons of full stop codon
                    # selection_df["ID:ptc index"].append(3*(ptc-1)) # <- replaced
                    selection_df["ID:ptc index"].append(3*(ptc-1)+2) # <- added

                else:
                    selection_df["ID:ptc index"].append(None)    
                    removed_ptcs.append(df.iloc[i].loc["HGVSp"])
                    print("removed", df.iloc[i].loc["HGVSp"])


        bar.next()
    bar.finish()
    
    selection  = pd.DataFrame.from_dict(selection_df)
    init_shape = selection.shape[0]
    selection  = selection[~selection["ID:ptc index"].isna()]

    print("<", init_shape-selection.shape[0], "unrecognized ptc indices were removed.")
    for removed_ptc in removed_ptcs:
        print(removed_ptc)

    # index needs to be changed to int, else causing downstream error @get_info
    selection = selection.astype({"ID:mutation index": "int32", "ID:ptc index": "int32"})
    return selection


# tcga dataset still contains non-stop mutations that are filtered out at get_info
def get_tcga_selection(df, params):
    # marked (<-) added on 250424 (shifted from prepare_data for consistency with other data types)
    df = df.sort_values(by=["variant_id"]) # <- added
    
    bar = IncrementalBar(set_bar("selecting TCGA data"), max=df.shape[0])
    selection_df = {"ID:case id": [],  "ID:gene id": [], "ID:gene symbol": [], "ID:project": [], "ID:sample id": [], "ID:transcript id": [], "ID:variant id": [], #"ID:alt. transcript id": [], 
                    "ID:chromosome": [], "ID:ptc reads": [], "ID:noptc reads": [], "ID:HGVSc": [], "ID:HGVSp": [], "ID:variant classification": [], "ID:cnv total": [], "ID:cnv minor": [], "ID:cnv avg": [],
                    "ID:t_ref_count": [], "ID:t_alt_count": [], "ID:VAF": [],
                    "ID:1000G_AF": [], "ID:PolyPhen": [], "ID:SIFT": [], 
                    "ID:noptc counts": [], "ID:noptc variance": [],
                    "ID:abs. mutation index": [], "ID:cds size": [], "ID:mutation index": [], "ID:wt base": [], "ID:ptc index": [],
                    "ID:cds size": [], "ID:aggregated variants": [], "LABEL:NMD score": [], "LABEL:alt. NMD score": []}
     
    removed_ptcs = []
    scores     = []
    last_variant_id = df.iloc[0].loc["variant_id"]

    # add row (last row is not considered)
    df.loc[len(df.index)] = [None for _ in df.columns]

    test = 0

    for i in range(df.shape[0]):
        if ((params["full_output"] == False and df.iloc[i].loc["variant_id"] != last_variant_id or i == df.shape[0]-1)
            or params["full_output"] == True):
            if len(scores) > 0:
                avg_score = np.average(scores)
                if params["full_output"] == True and len(scores) > 1: print("< warning. labels are averaged.")
                index = i-1
                #print("i", i, df.shape[0], df.iloc[index].loc["variant_id"], avg_score, len(scores), scores)

                if params["full_output"] == True or params["weighed_output"] == False:    selections = 1
                elif params["full_output"] == False and params["weighed_output"] == True: selections = len(scores)

                for _ in range(selections):
                    if pd.isna(df.iloc[index].loc["HGVSp"]) == False:
                        if get_position_from_mutation(df.iloc[index].loc["HGVSc"]) != None:
                            selection_df["ID:case id"].append(df.iloc[index].loc["case_id"])
                            selection_df["ID:gene id"].append(df.iloc[index].loc["Gene"])
                            selection_df["ID:gene symbol"].append(df.iloc[index].loc["SYMBOL"])
                            selection_df["ID:project"].append(df.iloc[index].loc["project"])
                            selection_df["ID:sample id"].append(df.iloc[index].loc["Tumor_Sample_UUID"])
                            selection_df["ID:transcript id"].append(df.iloc[index].loc["Transcript_ID"])
                            selection_df["ID:variant id"].append(df.iloc[index].loc["variant_id"])
                            selection_df["ID:variant classification"].append(df.iloc[index].loc["Variant_Classification"]) # <- added on 250604

                            selection_df["ID:abs. mutation index"].append(df.iloc[index].loc["Start_Position"]-1)
                            selection_df["ID:cds size"].append(3*int(df.iloc[index].loc["Protein_position"].split("/")[1]))
                            selection_df["ID:chromosome"].append(df.iloc[index].loc["Chromosome"])
                            selection_df["ID:t_ref_count"].append(df.iloc[index].loc["t_ref_count"])
                            selection_df["ID:t_alt_count"].append(df.iloc[index].loc["t_alt_count"])
                            selection_df["ID:VAF"].append(df.iloc[index].loc["t_alt_count"]/(df.iloc[index].loc["t_ref_count"]+df.iloc[index].loc["t_alt_count"]))

                            selection_df["ID:HGVSc"].append(df.iloc[index].loc["HGVSc"])
                            selection_df["ID:HGVSp"].append(df.iloc[index].loc["HGVSp"])
                            selection_df["ID:mutation index"].append(get_numbers(df.iloc[index].loc["CDS_position"])[0]-1)

                            # changed on 241112
                            ptc = get_ptc(df.iloc[index].loc["HGVSp"])
                            if ptc != None:
                                # marked (<-) added / replaced on 250425 to calculate with shifted relative PTC index (+2) to avoid assignment to wrong exons of full stop codon
                                # selection_df["ID:ptc index"].append(3*(ptc-1)) # <- replaced
                                selection_df["ID:ptc index"].append(3*(ptc-1)+2) # <- added
                                #print("i", i, df.iloc[index].loc["variant_id"], df.iloc[index].loc["HGVSc"], 3*(ptc-1)+2, 3*(ptc-1))
                            
                            else:
                                selection_df["ID:ptc index"].append(None)    
                                removed_ptcs.append(df.iloc[index].loc["HGVSp"])

                            # marked (<-) added / replaced on 250425 to calculate with shifted relative PTC index (+2) to avoid assignment to wrong exons of full stop codon
                            # if selection_df["ID:ptc index"][-1] >= int(df.iloc[index].loc["CDS_position"].split("/")[1].replace("\'", "").replace("]", "")): # <- replaced
                            if selection_df["ID:ptc index"][-1] >= int(df.iloc[index].loc["CDS_position"].split("/")[1].replace("\'", "").replace("]", ""))+2: # <- added
                                test += 1
                                print(df.iloc[index].loc["variant_id"], selection_df["ID:ptc index"][-1], df.iloc[index].loc["Protein_position"].split("/")[1], "/", df.iloc[index].loc["CDS_position"],
                                      df.iloc[index].loc["HGVSc"], "/", df.iloc[index].loc["HGVSp"],
                                      df.iloc[index].loc["Reference_Allele"], selection_df["ID:mutation index"][-1], test)
                                
                             # if "analytical" is chosen, basal expression is inferred from project-specific average (TCGA)
                            if params["extract_expression"] == True:
                                if params["output_type"] == "analytical" and pd.isna(df.iloc[index].loc["RNASEQ_noptc_fpkm_unstranded"]) == True:
                                    expressions = df[df["project"] == df.iloc[index].loc["project"]]["RNASEQ_noptc_fpkm_unstranded"].tolist()
                                    selection_df["ID:noptc reads"].append(np.mean([expression for expression in expressions if pd.isna(expression) == False]))

                                else:
                                    selection_df["ID:noptc reads"].append(df.iloc[index].loc["RNASEQ_noptc_fpkm_unstranded"])

                                selection_df["ID:noptc counts"].append(df.iloc[index].loc["RNASEQ_noptc_fpkm_unstranded_counts"])    
                                selection_df["ID:noptc variance"].append(df.iloc[index].loc["RNASEQ_noptc_fpkm_unstranded_stats"])
                                selection_df["ID:ptc reads"].append(df.iloc[index].loc["RNASEQ_ptc_fpkm_unstranded"])
                                selection_df["ID:cnv total"].append(df.iloc[index].loc["RNASEQ_ptc_cnv_total"])
                                selection_df["ID:cnv minor"].append(df.iloc[index].loc["RNASEQ_ptc_cnv_minor"])
                                selection_df["ID:cnv avg"].append(df.iloc[index].loc["RNASEQ_ptc_cnv_avg"])

                            else:
                                selection_df["ID:noptc reads"].append(None)
                                selection_df["ID:noptc counts"].append(None)    
                                selection_df["ID:noptc variance"].append(None)
                                selection_df["ID:ptc reads"].append(None)
                                selection_df["ID:cnv total"].append(None)
                                selection_df["ID:cnv minor"].append(None)
                                selection_df["ID:cnv avg"].append(None)
                            
                            selection_df["ID:1000G_AF"].append(df.iloc[index].loc["1000G_AF"])
                            selection_df["ID:PolyPhen"].append(df.iloc[index].loc["PolyPhen"])
                            selection_df["ID:SIFT"].append(df.iloc[index].loc["SIFT"])

                            selection_df["ID:aggregated variants"].append(len(scores))
                            # changed on 241126
                            #selection_df["ID:wt base"].append(df.iloc[index].loc["Reference_Allele"][0])
                            selection_df["ID:wt base"].append(df.iloc[index].loc["Reference_Allele"])

                            if params["output_type"] == "analytical" and pd.isna(avg_score) == True:
                                selection_df["LABEL:NMD score"].append(None)
                                selection_df["LABEL:alt. NMD score"].append(None)

                            else:
                                selection_df["LABEL:NMD score"].append(avg_score)
                                selection_df["LABEL:alt. NMD score"].append(df.iloc[index].loc["alt. NMD score"])

                scores = []

            last_variant_id = df.iloc[i].loc["variant_id"]

        if df.iloc[i].loc["variant_id"] == last_variant_id:
            scores.append(df.iloc[i].loc["NMD score"])

        bar.next()
    bar.finish()
    selection  = pd.DataFrame.from_dict(selection_df)
    init_shape = selection.shape[0]
    selection  = selection[~selection["ID:ptc index"].isna()]
    print("<", init_shape-selection.shape[0], "unrecognized ptc indices were removed.")
    for removed_ptc in removed_ptcs:
        print(removed_ptc)
    print("test", test)
    return selection


def get_teran_selection(df, params):
    df = df.sort_values(by=["VARIANT_ID"])
    
    selection_df = {"ID:variant id": [], "ID:alt. variant id": [], "ID:gene id": [], "ID:gene symbol": [], "ID:sample id": [], "ID:subject id": [], "ID:tissue id": [], "ID:transcript id": [],
                    "ID:chromosome": [], "ID:abs. mutation index": [], "ID:cds size": [], "ID:mutation index": [], "ID:ptc index": [], "ID:wt base": [],
                    "ID:mutated base": [],
                    "ID:aggregated variants": [], "LABEL:NMD score": []}
        
    alt_counts     = []
    ref_counts     = []
    total_counts   = []

    # add row (last row is not considered)
    bar = IncrementalBar(set_bar("creating teran selection"), max=df.shape[0])
    df.loc[len(df.index)] = [None for _ in df.columns]

    last_variant_id = df.iloc[0].loc["VARIANT_ID"]

    for i in range(df.shape[0]):
        if ((params["full_output"] == False and df.iloc[i].loc["VARIANT_ID"] != last_variant_id or i == df.shape[0]-1)
            or params["full_output"] == True):
            if len(ref_counts) > 0:
                label = np.average(ref_counts)

                if params["full_output"] == True and len(ref_counts) > 1:
                    print("< warning. labels are averaged.")

                index = i-1

                if params["full_output"] == True or params["weighed_output"] == False:    selections = 1
                elif params["full_output"] == False and params["weighed_output"] == True: selections = len(alt_counts)

                for _ in range(selections): 
                    selection_df["ID:variant id"].append(df.iloc[index].loc["VARIANT_ID"])
                    selection_df["ID:alt. variant id"].append(df.iloc[index].loc["Gene"]+"_"+str(df.iloc[index].loc["Protein_position_only"]))
                    selection_df["ID:gene id"].append(df.iloc[index].loc["Gene"])
                    selection_df["ID:gene symbol"].append(df.iloc[index].loc["SYMBOL"])
                    selection_df["ID:sample id"].append(df.iloc[index].loc["SAMPLE_ID"])
                    selection_df["ID:subject id"].append(df.iloc[index].loc["SUBJECT_ID"])
                    selection_df["ID:tissue id"].append(df.iloc[index].loc["TISSUE_ID"])
                    selection_df["ID:transcript id"].append(df.iloc[index].loc["Feature"])

                    selection_df["ID:aggregated variants"].append(len(ref_counts))
                    selection_df["ID:chromosome"].append(df.iloc[index].loc["CHR"])
                    selection_df["ID:cds size"].append(3*df.iloc[index].loc["Protein_length"])
                    selection_df["ID:abs. mutation index"].append(df.iloc[index].loc["POS"]-1)
                    selection_df["ID:mutation index"].append(df.iloc[index].loc["CDS_position_only"]-1)
                    selection_df["ID:ptc index"].append(3*int((df.iloc[index].loc["CDS_position_only"]-1)/3)+2) # <- added
                    selection_df["ID:wt base"].append(df.iloc[index].loc["REF_ALLELE"])
                    selection_df["ID:mutated base"].append(df.iloc[index].loc["ALT_ALLELE"]) # <-
                    selection_df["LABEL:NMD score"].append(label)
                                             		
                alt_counts   = []
                ref_counts   = []
                total_counts = []

            last_variant_id = df.iloc[i].loc["VARIANT_ID"]

        if df.iloc[i].loc["VARIANT_ID"] == last_variant_id:
            alt_counts.append(df.iloc[i].loc["ALT_COUNT"])
            ref_counts.append(df.iloc[i].loc["REF_RATIO"])
            total_counts.append(df.iloc[i].loc["TOTAL_COUNT"])

        bar.next()
    bar.finish()

    selection = pd.DataFrame.from_dict(selection_df)
    return selection


def get_ucids_by_position(selection, knowngene):
    knowngene_dict = create_knowngene_dict_by_position(knowngene, stepsize=10000000)
    uc_ids         = []
    misses         = 0
    bar = IncrementalBar(set_bar("getting ucids by positions"), max=selection.shape[0])

    for i in range(selection.shape[0]):
        current_uc_ids = search_knowngene_dict_by_position(knowngene_dict, knowngene, chrs[selection.iloc[i].loc["ID:chromosome"][3::]],
                                                           selection.iloc[i].loc["ID:abs. mutation index"], stepsize=10000000)
        
        uc_ids.append(current_uc_ids)
        bar.next()
    bar.finish()

    selection["uc id"] = uc_ids
    print("<", misses, "missing values")
    return selection


def invert_sequence(sequence):
    bases = {"A": "T", "C": "G", "G": "C", "T": "A"}
    inverted_sequence = ""
    for i in range(len(sequence)):
        inverted_sequence += bases[sequence[len(sequence)-1-i]]

    return inverted_sequence


def _fill_df(selection, thread_index, genome, knowngene, knowngene_dict, params, hg_build, bar):
    if params["ptc_mode"][hg_build] == "absolute":
        mutation_identifier = "ID:abs. mutation index"
        ptc_identifier      = "ID:abs. ptc index"

    if params["ptc_mode"][hg_build] == "relative":
        mutation_identifier = "ID:mutation index"
        ptc_identifier      = "ID:ptc index"
    
    # initialize error report
    error_report = {key: {error: [] for error in [*params["error_handling"], "no_error"]} 
                    for key in ["+del1", "+delins1", "+dup1", "+ins1", "+nonsense", "-del1", "-delins1", "-dup1", "-ins1", "-nonsense",
                                "+del>1", "+delins>1", "+dup>1", "+ins>1", "-del>1", "-delins>1", "-dup>1", "-ins>1", "total", "+total", "-total"]}
    
    for i in range(thread_index[0], thread_index[1]+1, 1):
        if selection.iloc[i].loc["success"] == False and pd.isna(selection.iloc[i].loc[ptc_identifier]) == False and pd.isna(selection.iloc[i].loc[mutation_identifier]) == False:
            total_fail         = True
            appris_annotations = []
            fail_info          = []
            infos              = []

            for mapping_target in selection.iloc[i].loc["uc id"]:
                if "." in mapping_target:
                    mapping_target = mapping_target.split(".")[0]

                if mapping_target in knowngene_dict.keys():
                    j = knowngene_dict[mapping_target]

                    if "ID:mutated base" in selection.columns: mutated_base = selection.iloc[i].loc["ID:mutated base"]
                    else:                                      mutated_base = None

                    info, fails = get_info(knowngene, selection.iloc[i].loc[params["block_identifier"]], j, selection.iloc[i].loc[ptc_identifier],
                                           selection.iloc[i].loc[mutation_identifier], selection.iloc[i].loc["ID:wt base"], mutated_base,
                                           selection.iloc[i].loc["ID:cds size"], genome, params, ptc_mode=params["ptc_mode"][hg_build])
                        
                    appris_annotations.append(knowngene.iloc[j].loc["appris annotation"])
                    fail_info.append(str(i) + ", variant id: " + selection.iloc[i].loc["ID:variant id"] + ", transcript id: " + selection.iloc[i].loc["ID:transcript id"]
                                     + ", mapping_target: " + mapping_target + ", fails: " + json.dumps(fails))

                    if check_fails(fails, params) == False:
                        infos.append(info)
                        total_fail = False

                    # fill error report
                    error_report = fill_error_report(error_report, selection, fails, knowngene.iloc[j].loc["strand"], i)

            if total_fail == True:
                if params["verbose"] == True:
                    print("gene symbol", selection.iloc[i].loc["ID:gene symbol"], "pos", selection.iloc[i].loc[ptc_identifier])
                    for j in range(len(fail_info)):
                        print(fail_info[j])

            # the following section works only if cds size is given (not for custom WXS data)
            elif pd.isna(selection.iloc[i].loc["ID:cds size"]) == False:
                # find best match among possible solutions with priority order: identical cds length - PRINCIPAL isoform
                match_index     = None
                principal_index = None
                j = 0
                error = []

                while j < len(infos) and match_index == None:
                    error.append("j " + str(j) + " " + str(3*selection.iloc[i].loc["ID:cds size"]) + " " + str(infos[j]["FEATURE:total cds size"]) + " " + appris_annotations[j])
                    if selection.iloc[i].loc["ID:cds size"] == infos[j]["FEATURE:total cds size"]-3:
                        match_index = j # priority 1: indicated lengths match dictionary length

                    if appris_annotations[j] == "PRINCIPAL:0" or appris_annotations[j] == "PRINCIPAL:1":
                        principal_index = j # priority 2: component is PRICIPAL:0 (no other isoforms) or PRINCIPAL:1 (main isoform)

                    j += 1

                if match_index == None and principal_index != None:
                    match_index = principal_index

                # priority 3: priorities 1 and 2 do not apply but selection is unambiguous
                if match_index == None and len(infos) == 1:
                    match_index = 0

                if match_index != None:
                    selection.at[selection.index[i], "success"] = True

                    for key, val in infos[match_index].items():
                        selection.at[selection.index[i], key] = val
            
            # if cds size is not known, first info is taken arbitrarily
            else:
                selection.at[selection.index[i], "success"] = True

                for key, val in infos[0].items():
                    selection.at[selection.index[i], key] = val

        if params["verbose"] == False: bar.next()

    bar.finish()
    print("< error report")
    for key1 in error_report:
        for key2 in error_report[key1]:
            print(key1, key2, len(error_report[key1][key2]))

    return selection


def fill_df(selection, genome, knowngene, knowngene_dict, hg_build, params):
    thread_index = split_index(selection.shape[0], params["threads"])
    
    threads      = []
    bar          = IncrementalBar(set_bar("extracting sequence info"), max=selection.shape[0])
    for i in range(len(thread_index)):
        thread = threading.Thread(target=_fill_df, args=(selection, thread_index[i], genome, knowngene, knowngene_dict, params, hg_build, bar))
        threads.append(thread)
        thread.start()

    for i, thread in enumerate(threads):
        thread.join()

    bar.finish()
    return selection