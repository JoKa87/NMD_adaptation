import pandas as pd
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from apply_lindeboom_utils import *
from shared_utils import *


params = {
         "block_identifier"      : "VARIANT_ID", # "ID:variant id", "VARIANT_ID"
         "blockwise"             : False,
         "data_dir"              : parent_dir+r"\data",
         "datatype"              : "cuomo19", # "cuomo19", "cuomo38", "custom", "tcga", "teran" # datatype 'kim' added 250410
         "dictionary_file"       : "",
         "decision_threshold"    : 0.65,
         "os_sep"                : "\\",
         "lindeboom_file"        : "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf", # "hg19_NMDetectiveA_Lindeboom_et_al.v2.gtf", "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf",
         "mapping_target"        : "transcript id", # gene symbol, transcript id
         "mode"                  : "relative", # "absolute", "relative"
         "separator"             : " ",
         "target_file"           : "model_test_variants.txt" # "mmc4.txt" "model_test_variants.txt"
         }

# marked (<-) added/removed on 250429 to allow calculation of cuomo hg38 along with hg19
# if params["datatype"] == "cuomo": # <- removed
if params["datatype"] == "cuomo19": # <- added
    params["block_identifier"] = "ID:variant id"
    hg_build_dir               = "hg19"
    params["lindeboom_file"]   = "hg19_NMDetectiveA_Lindeboom_et_al.v2.gtf"
    params["separator"]        = ","
    knowngene_fname            = "hg19_knownGene.txt"

# marked (<-) added/removed on 250429 to allow calculation of cuomo hg38 along with hg19
if params["datatype"] == "cuomo38": # <- added from here
    params["block_identifier"] = "ID:variant id"
    hg_build_dir               = "hg38.p14"
    params["lindeboom_file"]   = "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"
    params["separator"]        = ","
    knowngene_fname            = "hg38_knownGene.txt"
    # <- until here

if params["datatype"] == "custom":
    params["block_identifier"] = "ID:variant id"
    hg_build_dir               = "hg38.p14"
    params["lindeboom_file"]   = "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"
    params["separator"]        = ","
    knowngene_fname            = "hg38_knownGene.txt"

# marked (<-) added/removed on 250410 to integrate data from Kim et al.
# if params["datatype"] == "tcga": # <- removed
if params["datatype"] == "kim" or params["datatype"] == "tcga": # <- added
    params["block_identifier"] = "ID:variant id"
    hg_build_dir               = "hg38.p14"
    params["lindeboom_file"]   = "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"
    params["separator"]        = ","
    knowngene_fname            = "hg38_knownGene.txt"

if params["datatype"] == "teran":
    params["block_identifier"] = "VARIANT_ID"
    hg_build_dir               = "hg38.p14"
    params["lindeboom_file"]   = "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"
    params["separator"]        = " "
    knowngene_fname            = "hg38_knownGene.txt"


# load genome data
genome = load_split_genome(params["data_dir"]+params["os_sep"]+hg_build_dir)
print("<", hg_build_dir, "loaded.")

# load lindeboom file
with open(params["data_dir"]+params["os_sep"]+params["lindeboom_file"], "r") as f:
    lines = f.readlines()

# create gene dictionary
lindeboom_dictionary = create_lindeboom_dictionary(lines)
print("< gene dictionary created.")

# load target file
targets = pd.read_csv(params["data_dir"]+params["os_sep"]+params["target_file"], delimiter=params["separator"])
print("< targets loaded.")

# marked (<-) added/removed on 250429 to allow calculation of cuomo hg38 along with hg19
if params["datatype"] == "cuomo38":
    targets = make_liftover(targets) # <- added
    targets = targets[~targets["ID:abs. mutation index"].isna()] # <- added

# get ucids for targets
# marked (<-) added/removed on 250429 to allow calculation of cuomo hg38 along with hg19
# if params["datatype"] == "cuomo": # <- removed
if "cuomo" in params["datatype"]: # <- added
    targets["uc id"] = [[targets.iloc[i].loc["ID:uc id"]] for i in range(targets.shape[0])]

else:
    if params["mapping_target"] == "transcript id":
        # load knowGene dictionary
        with open(params["data_dir"]+params["os_sep"]+knowngene_fname, 'r') as _:
            knowngene = pd.read_csv(params["data_dir"]+params["os_sep"]+knowngene_fname, delimiter=",", index_col=False).sort_values(by=["chr", "cdsstart"])
        
        knowngene["transcript id"] = [knowngene.iloc[i].loc["transcript id"].split(".")[0] for i in range(knowngene.shape[0])]
        if params["datatype"] == "teran": ucids = append_df_with_mapping2([targets, knowngene], "Feature", "transcript id", "uc id", verbose=True)
        else:                             ucids = append_df_with_mapping2([targets, knowngene], "ID:transcript id", "transcript id", "uc id", verbose=True)
        targets["uc id"] = ucids


# get lindeboom predictions
targets, report = get_lindeboom_predictions(targets, lindeboom_dictionary, genome, params)
print("< predictions calculated.", report)
plt.hist(targets["LABEL:NMD score"], bins=40, histtype="step")
plt.show()

# evaluate
if "LABEL:NMD score" in targets.columns: evaluate(targets, params)

# print predictions
print_results(targets, params)