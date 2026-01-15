import os
import pandas as pd
import platform

from prepare_data_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


# parameters (default setting: input from Teran et al. 2021 (mmc4.txt))
params = {
         "block_identifier":  "ID:variant id",
         "cohorts":           5,
         "datatype":          "teran", # "control" "mskcc" "teran" "tcga"
         "delimiters":        None,
         "dictionary_file":   "gene_symbol_ucid_dictionary_with_synonyms_manual_edit.txt",
         "data_dir":          parent_dir+r"\data",
         "error_handling":    { # defines handling of errors during sequence information extraction (True means excluding errors)
                              "chromosome_incompatible"         : True, # in absolute mode, alternative chromosomes cannot be used
                              "chromosome_not_found"            : True,
                              # this error is taken out because some entries are outdated
                              "inconsistent_cds_size"           : True, # only applied if value is provided
                              "motif_start_not_found"           : True,
                              "mutated_wt_base_not_found"       : False,
                              # marked (<-) added on 260109 to treat error in novel get_read_through_motif function
                              "mutated_stop_codon_not_found"    : False,
                              # marked (<-) added to account for no cds (relevant for calculation of "FEATURE:cds EJC density" and "FEATURE:downstream cds EJC density")
                              "no_cds"                          : True, # <-
                              "no_exons"                        : True,
                              "no_start"                        : True,
                              "no_stop"                         : True,
                              "ptc_cds_index_not_found"         : True,
                              "ptc_equals_stop_codon"           : False,
                              "ptc_exon_index_not_found"        : True,
                              "ptc_not_in_cds"                  : True,
                              "ptc_not_in_exon"                 : True,
                              "ptc_in_utr"                      : True
                              },
         "extract_expression": True,
         "features":          ["FEATURE:total exon size", "FEATURE:total cds size", "FEATURE:5'utr-size", "FEATURE:3'utr-size", "FEATURE:abs. ptc index", 
                               "FEATURE:cds", "FEATURE:exons", "FEATURE:upstream exons", "FEATURE:downstream exons", 
                               # marked (<-) features added on 250425 to account for cds exons separately
                               "FEATURE:cds exons", "FEATURE:upstream cds exons", "FEATURE:downstream cds exons", # <- added
                               "FEATURE:EJC density", "FEATURE:downstream EJC density",
                               # marked (<-) features added on 250425 to account for cds exons separately
                               "FEATURE:cds EJC density", "FEATURE:downstream cds EJC density", # <- added
                               "FEATURE:ptc upstream distance", "FEATURE:ptc downstream distance", "FEATURE:dist. from last EJC",
                               "FEATURE:ptc cds position", "FEATURE:ptc exon position", "FEATURE:ptc-wt stop codon distance",
                               "FEATURE:last exon", "FEATURE:last cds exon", "FEATURE:50 nt to last EJC", "FEATURE:50 nt to last cds EJC", "FEATURE:ptc exon size", # <- added

                               "FEATURE:total start count cds", "FEATURE:upstream start ptc count cds", "FEATURE:downstream start ptc count cds",
                               "FEATURE:first upstream start ptc distance cds", "FEATURE:first downstream start ptc distance cds",
                               "FEATURE:total start count ptc cds", "FEATURE:upstream start ptc count ptc cds", "FEATURE:downstream start ptc count ptc cds",
                               "FEATURE:first upstream start ptc distance ptc cds", "FEATURE:first downstream start ptc distance ptc cds",

                               "FEATURE:total GC count exons", "FEATURE:upstream GC ptc count exons", "FEATURE:downstream GC ptc count exons",
                               "FEATURE:total GC count ptc exon", "FEATURE:upstream GC ptc count ptc exon", "FEATURE:downstream GC ptc count ptc exon",
                               # marked (<-) features added on 260109 to integrate read-through motifs
                               "FEATURE:read-through+1", "FEATURE:read-through+2",
                               "FEATURE:5'utr", "FEATURE:ptc cds", "FEATURE:ptc", "FEATURE:5'ptc", "FEATURE:3'ptc", "FEATURE:5'ejc", "FEATURE:3'ejc",
                               "FEATURE:3'utr", "FEATURE:all cds", "FEATURE:all exons"],
         "file_tags":         None,
         "full_output":       True, # True means no averaging over identical variants
         "hg_build":          ["hg38"], #["hg19", "hg38"], priority list for hg builds
         "mapping_target":    {"hg38": "transcript id"}, #{"hg19": "by position", "hg38": "transcript id"}, # {"hg38": "transcript id"} # transcript id, by_position
         "motifs":            {
                              "start":  ["ATG"],
                              "GC":     ["G", "C"],
                              },
         "motifs_combined":   False,
         "motif_distance":    {"start": True, "GC": False},
         "motif_inframe":     {"start": True, "GC": False},
         "motif_relative":    {"start": False, "GC": True},
         "os_sep":             "//",
         "outfname":           "teran_read_through",
         "output_type":        "analytical", # if "analytical" is chosen, missing NMD scores are not dropped, basal expression is inferred from project-specific average (TCGA)
         "pattern_cutoff":     5, # inclusive threshold for the minimum number of bases to determine a pattern (e.g. GC content); if below, default value is used instead
         "ptc_cutoff":         None,
         "ptc_mode":           {"hg38": "relative"}, #{"hg19": "absolute", "hg38": "relative"}, # {"hg38": "relative"} # absolute, relative
         "randomize":          False,
         "threads":            1,
         "transcript_column":  "alt. transcript id", # applicable to TCGA dataset (alternative transcript id, based on isoform analysis)
         "transcript_filter":  "ext. appris annotation", # applicable to TCGA dataset (alternative transcript id, based on isoform analysis)
         "verbose":            False,
         "weighed_output":     False
         }


if platform.system() == "Windows":                                         params["os_sep"] = "\\"
if params["motifs_combined"] == True:                                      params = append_features(params)
if params["full_output"] == True and params["weighed_output"] == True:     print("< warning. full_output=True and weighed_output=True are mutually exclusive.")
if params["output_type"] == "analytical" and params["datatype"] != "tcga": print("< warning. output_type is set to 'analytical' which is applied to tcga data only.")

if len([params["ptc_mode"][key] for key in params["ptc_mode"] if params["ptc_mode"][key] == "absolute"]) > 0:
    print("< absolute mode is currently not supported.")
    exit()

# default settings depending on datatypes
if params["datatype"] == "control":
    params["delimiters"] = [","]
    params["file_tags"]  = ["control_raw.txt"]

if params["datatype"] == "mskcc":
    params["delimiters"] = ["\t"]
    params["file_tags"]  = ["msk_chord_data_mutations.txt"]

if params["datatype"] == "tcga":
    params["delimiters"] = ["\t"]
    params["file_tags"]  = ["cptac3_raw.txt"]

if params["datatype"] == "teran":
    params["delimiters"] = [" ", ","]
    params["file_tags"]  = ["mmc4.txt", "mmc4_lindeboom_predictions.txt"]


def main():
    # load data
    data = []     
    for i in range(len(params["file_tags"])):
        file_tag = params["file_tags"][i]
        with open(params["data_dir"]+params["os_sep"]+file_tag, 'r') as _:
            data.append(pd.read_csv(params["data_dir"]+params["os_sep"]+file_tag, delimiter=params["delimiters"][i], index_col=False))

    if params["datatype"] == "control":
        for hg_build in params["ptc_mode"]:
            if params["ptc_mode"][hg_build] == "absolute":
                print("< absolute mode not working for control dataset. set to 'relative'.")
                params["ptc_mode"][hg_build] = "relative"

    # apply filtering steps to MSKCC dataset
    if params["datatype"] == "mskcc":
        data[0] = data[0][~data[0]["HGVSp"].isna()]
        data[0] = data[0][data[0]["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"])]

    # apply filtering steps to TCGA dataset
    if params["datatype"] == "tcga":
        for hg_build in params["ptc_mode"]:
            if params["ptc_mode"][hg_build] == "absolute":
                print("< absolute mode not working for TCGA dataset. set to 'relative'.")
                params["ptc_mode"][hg_build] = "relative"

        data[0] = data[0][~data[0]["Protein_position"].isna()]
        if params["output_type"] != "analytical": data[0].dropna(subset=["NMD score"], inplace=True)

        # if no label is passed, create placeholders
        if "NMD score" not in data[0].columns:
            data[0]["NMD score"]      = [None for _ in range(data[0].shape[0])]
            data[0]["alt. NMD score"] = [None for _ in range(data[0].shape[0])]


    # create a list of non-redundant variant ids
    # marked (<-) newly added on 251027 to integrate control data
    if params["datatype"] == "control": selection = get_control_selection(data[0], params)
    if params["datatype"] == "mskcc":   selection = get_mskcc_selection(data[0], params)
    if params["datatype"] == "tcga":    selection = get_tcga_selection(data[0], params)
    if params["datatype"] == "teran":   selection = get_teran_selection(data[0], params)
    print("< selection created.")


    # append selection to prepare extraction
    for col in params["features"]:
        if col not in selection.columns:
            selection[col] = [None for _ in range(selection.shape[0])]

    selection["success"] = [False for _ in range(selection.shape[0])]
    
    for hg_build in params["hg_build"]:
        print("<", hg_build)

        if hg_build == "hg19":
            hg_build_dir    = "hg19"
            knowngene_fname = "hg19_knownGene_appended.txt"

        if hg_build == "hg38":
            hg_build_dir    = "hg38.p14"
            knowngene_fname = "hg38_knownGene_appended.txt"


        # load genome
        genome = load_split_genome(params["data_dir"]+params["os_sep"]+hg_build_dir, os_sep=params["os_sep"])
        

        # load knowGene dictionary
        with open(params["data_dir"]+params["os_sep"]+knowngene_fname, 'r') as _:
            knowngene = pd.read_csv(params["data_dir"]+params["os_sep"]+knowngene_fname, delimiter=",", index_col=False).sort_values(by=["chr", "cdsstart"])
        
        knowngene_dict = create_knowngene_dict(knowngene, "uc id")
        print("< KnownGene dictionary created.")
        
        if params["mapping_target"][hg_build] == "by position":
            selection = get_ucids_by_position(selection, knowngene)

        if params["mapping_target"][hg_build] == "gene id":
            knowngene["transcript id"] = [knowngene.iloc[i].loc["gene id"].split(".")[0] for i in range(knowngene.shape[0])]
            ucids                      = append_df_with_mapping2([selection, knowngene], "ID:gene id", "gene id", "uc id", verbose=False)
            selection["uc id"]         = ucids
        
        if params["mapping_target"][hg_build] == "transcript id":
            knowngene["transcript id"] = [knowngene.iloc[i].loc["transcript id"].split(".")[0] for i in range(knowngene.shape[0])]
            ucids                      = append_df_with_mapping2([selection, knowngene], "ID:transcript id", "transcript id", "uc id", verbose=True)
            selection["uc id"]         = ucids

        # fill dataframe
        selection = fill_df(selection, genome, knowngene, knowngene_dict, hg_build, params)

    # apply filters
    pd.set_option('display.max_columns', None)
    total_selection_size = selection.shape[0]
    selection = selection[selection["success"] == True]

    print("< filling dataframe successful for", selection.shape[0], "of", total_selection_size, "variants")

    # remove temporary index, rename additional features
    selection = selection.drop(["success", "uc id"], axis=1)
    selection.to_csv(path_or_buf=params["data_dir"]+params["os_sep"]+params["outfname"] + ".txt", sep=",", index=False)
    print("< file printed.")


    #create cohort file if cohort list was loaded
    if len(data) == 1 and params["randomize"] == True:
        selection = get_random_cohorts(selection, params["cohorts"])
        selection.to_csv(path_or_buf=params["data_dir"]+params["os_sep"]+params["outfname"]+"_N" + str(params["cohorts"]) + ".txt",
                         sep=",", index=False)

    if len(data) == 2:
        selection = append_lindeboom_predictions(selection, data[1])
        if params["randomize"] == True: selection = get_random_cohorts(selection, params["cohorts"])
        selection.to_csv(path_or_buf=params["data_dir"]+params["os_sep"]+params["outfname"]+"_lindeboom_N" + str(params["cohorts"]) + ".txt",
                         sep=",", index=False)

    # store parameters of current session
    params_fname = params["outfname"] + "_params.json"
    with open(params["data_dir"]+params["os_sep"]+params_fname, "w") as f:
        f.write(json.dumps(params, indent=4))

if __name__ == '__main__':
    main()