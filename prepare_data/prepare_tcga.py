import pandas as pd
import sys

from prepare_mskcc_utils import *
sys.path.insert(0, r"C:\Programming\Translational_genomics\NMD_analysis\NMD_activity")
from tcga_tools_utils import *


tcga_profile = {
                "data_dir"          : r"C:\Programming\Translational_genomics\NMD_analysis\data",
                "fnames"            : {
                                        "patient"             : "TCGA-CDR-SupplementalTableS1_extracted_appended.txt",
                                        "ptcs"                : "tcga.txt",
                                        },
                "gender_target"     : "GENDER",
                "label_ids"         : ["LABEL:OS", "LABEL:OS.time", "LABEL:DSS", "LABEL:DSS.time", "LABEL:DFI", "LABEL:DFI.time", "LABEL:PFI", "LABEL:PFI.time"],
                "outfname"          : "tcga_zeros_included_filtered.txt",
               }


profile = tcga_profile

# parameters (default setting: input from Teran et al. 2021 (mmc4.txt))
params = {
         "cohorts"              : 5,
         "data_dir"             : profile["data_dir"],
         "expression_target"    : "median_fpkm_unstranded",
         "expression_path"      : r"C:\Programming\Translational_genomics\NMD_analysis\data\tcga_expressions_info.txt",
         "feature_targets"      : [],
         "filter_by_msk"        : True,
         "fnames"               : profile["fnames"],
         "gender_target"        : profile["gender_target"],
         "label_ids"            : profile["label_ids"],
         "msk_path"             : r"C:\Programming\Translational_genomics\NMD_analysis\data\msk_chord_appended.txt",
         "os_sep"               : "\\",
         "outfname"             : profile["outfname"],
         "processed_fname"      : "tcga_zeros_included.txt",
         # path of reference file to avoid redundance of patients in multiple datasets
         "reference_filter"     : False,
         "score_threshold"      : {"FEATURE:escape": (0, 0.57), "FEATURE:target": (0.64, 1)},
         "split_mode"           : "",
         "zeros_included"       : True
         }


def count_by_threshold(scores, score_threshold):
    return len([score for score in scores if score >= score_threshold[0] and score <= score_threshold[1]])


def main():
    # load patient data
    patient_data = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["patient"], sep=",")

    # load extracted ptcs
    ptcs         = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

    # load msk data
    if params["filter_by_msk"] == True:
        msk  = pd.read_csv(params["msk_path"], sep=",")
        ptcs = ptcs[ptcs["ID:transcript id"].isin(msk["ID:transcript id"])]

    # load expression data
    expression_data    = pd.read_csv(params["expression_path"], sep=",")
    expressions        = append_df_with_mapping([ptcs, expression_data], "ID:gene symbol", "gene_name", params["expression_target"],
                                                 set_bar("mapping expression data"), non_redundant=True, reverse=False)
    median_expression  = expression_data[params["expression_target"]].median()

    ptcs["expression"] = [float(expression) if expression != "-" else median_expression for expression in expressions] # <- added on 251028

    # prepare analysis
    updated_cols = {}
    for col in patient_data.columns:
        if col in params["label_ids"]:         updated_cols[col] = "LABEL:"+col
        elif col in params["feature_targets"]: updated_cols[col] = "FEATURE:"+col
        else:                                  updated_cols[col] = "ID:"+col

    patient_data = patient_data.rename(columns=updated_cols)

    # add categories
    ptcs["index"] = [i for i in range(ptcs.shape[0])]
    mapped_index  = append_df_with_mapping([patient_data, ptcs], "ID:case id", "ID:case id", "index",
                                            set_bar("mapping TCGA data"), non_redundant=False, reverse=True)
    
    extracted_mutations = {**{"FEATURE:expression":             [None for _ in mapped_index]},
                           **{"FEATURE:frameshift":             [None for _ in mapped_index]},
                           **{"FEATURE:frameshift_mutations":   [None for _ in mapped_index]},
                           **{"FEATURE:escape":                 [None for _ in mapped_index]},
                           **{"FEATURE:prediction":             [None for _ in mapped_index]},
                           **{"FEATURE:target":                 [None for _ in mapped_index]},
                           **{"FEATURE:ptc_mutations":          [None for _ in mapped_index]}}

    bar = IncrementalBar(set_bar("processing mutation data"), max=len(mapped_index))
    for i in range(len(mapped_index)):
        selected_ptcs = ptcs.iloc[[int(j) for j in mapped_index[i]]]
        selected_ptcs = selected_ptcs[selected_ptcs["ID:variant classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"])]

        if selected_ptcs.shape[0] > 0:
            extracted_mutations["FEATURE:expression"][i]           = selected_ptcs["expression"].mean()
            extracted_mutations["FEATURE:frameshift"][i]           = selected_ptcs["FEATURE:frameshift"].mean()
            extracted_mutations["FEATURE:frameshift_mutations"][i] = selected_ptcs["FEATURE:frameshift_mutations"].mean()
            
            # predictions
            scores = [selected_ptcs.iloc[j].loc["FEATURE:prediction"] for j in range(selected_ptcs.shape[0])
                        if pd.isna(selected_ptcs.iloc[j].loc["FEATURE:prediction"]) == False]

            if len(scores) > 0:
                extracted_mutations["FEATURE:escape"][i]     = count_by_threshold(scores, params["score_threshold"]["FEATURE:escape"]) / len(scores)
                extracted_mutations["FEATURE:prediction"][i] = np.mean(scores)
                extracted_mutations["FEATURE:target"][i]     = count_by_threshold(scores, params["score_threshold"]["FEATURE:target"]) / len(scores)

            elif params["zeros_included"] == True:
                extracted_mutations["FEATURE:escape"][i]     = 0
                extracted_mutations["FEATURE:target"][i]     = 0

            
        elif params["zeros_included"] == False:
            extracted_mutations["FEATURE:frameshift_mutations"][i] = 0

        elif params["zeros_included"] == True:
            extracted_mutations["FEATURE:escape"][i]               = 0
            extracted_mutations["FEATURE:expression"][i]           = 0
            extracted_mutations["FEATURE:frameshift"][i]           = 0
            extracted_mutations["FEATURE:frameshift_mutations"][i] = 0
            extracted_mutations["FEATURE:target"][i]               = 0

        extracted_mutations["FEATURE:ptc_mutations"][i] = selected_ptcs.shape[0]
        
        bar.next()
    bar.finish()

    for feature in extracted_mutations:
        if None not in extracted_mutations[feature]: print(feature, np.mean(extracted_mutations[feature]))
        patient_data.insert(patient_data.shape[1], feature, extracted_mutations[feature])

    print("< total values", patient_data.shape[0],
            ", existing frameshifts", patient_data[patient_data["FEATURE:frameshift"] > 0].shape[0],
            ", existing prediction values", patient_data[~patient_data["FEATURE:prediction"].isna()].shape[0])
    
    patient_data = get_random_cohorts(patient_data, params["cohorts"], targets=["ID:case id"])
    patient_data = patient_data.replace({"#NV": None})
    patient_data.to_csv(path_or_buf=params["data_dir"]+params["os_sep"]+params["outfname"], sep=",", index=False)
    

if __name__ == '__main__':
    main()