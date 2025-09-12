import os
import pandas as pd

from prepare_mskcc_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


mskcc_chord_profile = {
                       "data_dir"          : parent_dir+r"\data",
                       "fnames"            : {
                                              # <- "markers" added on 250613
                                              "markers"             : {
                                                                       "ca_15-3"               : {
                                                                                                  "fname"           : "data_timeline_ca_15-3_labs.txt",
                                                                                                  "target_col"      : "RESULT"
                                                                                                 },
                                                                       "ca_19-9"               : {
                                                                                                  "fname"           : "data_timeline_ca_19-9_labs.txt",
                                                                                                  "target_col"      : "RESULT"
                                                                                                 },
                                                                       "cea"                   : {
                                                                                                  "fname"           : "data_timeline_cea_labs.txt",
                                                                                                  "target_col"      : "RESULT"
                                                                                                 },
                                                                       "gleason"               : {
                                                                                                  "fname"           : "data_timeline_gleason.txt",
                                                                                                  "target_col"      : "GLEASON_SCORE"
                                                                                                 },
                                                                       "mmr"                   : {
                                                                                                  "fname"           : "data_timeline_mmr.txt",
                                                                                                  "target_col"      : "MMR_ABSENT"
                                                                                                 },
                                                                       "pdl1"                  : {
                                                                                                  "fname"           : "data_timeline_pdl1.txt",
                                                                                                  "target_col"      : "PDL1_POSITIVE"
                                                                                                 },
                                                                       "psa"                   : {
                                                                                                  "fname"           : "data_timeline_psa_labs.txt",
                                                                                                  "target_col"      : "RESULT"
                                                                                                 }
                                                                      },
                                              "mutations"           : "data_mutations.txt",
                                              "patient"             : "data_clinical_patient_edited.txt",
                                              "ptcs"                : "msk_chord_variants.txt",
                                              "sample"              : "data_clinical_sample_edited.txt",
                                              "specimen"            : "data_timeline_specimen.txt",
                                              "treatment"           : "data_timeline_treatment.txt",
                                             },
                       "gender_target"     : "GENDER",
                       "label_ids"         : ["OS_MONTHS", "OS_STATUS"],
                       "outfname"          : "msk_chord_survival_analysis.txt",
                       "reference_path"    : r"C:\Programming\Translational_genomics\NMD_analysis\data\tmb_mskcc_2018\data_clinical_sample_edited.txt"
                      }

profile = mskcc_chord_profile

# parameters (default setting: input from Teran et al. 2021 (mmc4.txt))
params = {
         "adjustment_params"    : { # <- added on 250527 to compare msk subsets based on ptc mutation occurrence,
                                   "min_expression" : 0,
                                   "mode"           : "ID:psa_diff", # "immuno" "metastasis" 
                                   "repeats"        : 5
                                  },
         "cancer_wise"          : False,
         "cohorts"              : 5,
         "data_dir"             : profile["data_dir"],
         "expression_target"    : "median_fpkm_unstranded",
         "expression_path"      : r"C:\Programming\Translational_genomics\NMD_analysis\data\tcga_avg_expressions.txt",
         "feature_targets"      : [],
         "fnames"               : profile["fnames"],
         "gender_target"        : profile["gender_target"],
         "label_ids"            : profile["label_ids"],
         # "adjust_to_ptcs" "append_patient_data" "append_treatment" "append_treatment_by_time" "append_treatment_by_time2" "filter_mutations"
         # "map_expressions" # <- added on 250527 to compare msk subsets based on ptc mutation occurrence
         # "prepare_analysis" "split_by_types" "test_patient_data" "test_ptcs"
         "mode"                 : "prepare_analysis",
         "os_sep"               : "\\",
         "outfname"             : profile["outfname"],
         "processed_fname"      : "msk_chord_survival_analysis.txt",
         # path of reference file to avoid redundance of patients in multiple datasets
         "reference_filter"     : False,
         "reference_path"       : profile["reference_path"],
         "score_threshold"      : {"FEATURE:escape": (0, 0.57), "FEATURE:target": (0.64, 1)},
         "split_mode"           : "",
         "zeros_included"       : True
         }


def main():
    pmu = Prepare_mskcc_utils(params)

    # marked (<-) function added on 250527
    if params["mode"] == "adjust_to_ptcs":
        # load mutation data
        mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        # create msk datasets
        pmu.adjust_to_ptcs(mutations)


    if params["mode"] == "append_patient_data":
        # load mutation data
        mutations    = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        # load patient data
        patient_data = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["patient"], sep="\t")

        # load sample data
        sample_data  = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["sample"], sep="\t")

        # append patient data
        pmu.append_patient_data(mutations, patient_data, sample_data)


    if params["mode"] == "append_treatment":
        # load mutation data
        mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+params["processed_fname"], sep=",")

        # load treatment data
        treatment = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["treatment"], sep="\t")

        mutations = pmu.append_treatment(mutations, treatment)
        mutations.to_csv(params["data_dir"]+params["os_sep"]+params["processed_fname"].split(".")[0]+"_appended.txt", sep=",", index=False)


    # marked (<-) added on 250521 to allow appending features to processed data as input for survival analysis (can be merged later with append_treatment_by_time)
    # if params["mode"] == "append_treatment_by_time": # <- removed
    if params["mode"] == "append_treatment_by_time" or params["mode"] == "append_treatment_by_time2": # <- added from here
        if params["mode"] == "append_treatment_by_time":
            mutations_fname = params["fnames"]["ptcs"]
            patient_col     = "ID:patient id"
            sample_col      = "ID:sample id"

        if params["mode"] == "append_treatment_by_time2":
            mutations_fname = params["processed_fname"]
            patient_col     = "ID:PATIENT_ID"
            sample_col      = "ID:SAMPLE_ID"

        # <- until here

        # load mutation data
        # marked (<-) added on 250521 to allow appending features to processed data as input for survival analysis (can be merged later with append_treatment_by_time)
        # mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",") # <- removed
        mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+mutations_fname, sep=",") # <- added

        # load sequence timeline to exclude any samples with time of sequencing other than 0
        specimen  = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["specimen"], sep="\t")
        specimen  = specimen[specimen["START_DATE"] == 0]

        # test whether different sample ids exist for identical patients at specimen time 0
        # important to exclude the possibility of redundant samples fed into analyze_predictions
        patients = specimen.drop_duplicates(subset=["PATIENT_ID"])["PATIENT_ID"].tolist()
        test     = [specimen[specimen["PATIENT_ID"] == patient].drop_duplicates(subset=["SAMPLE_ID"]).shape[0] for patient in patients]
        if np.unique(test).shape[0] != 1 or np.unique(test)[0] != 1:
            print("< warning @append_treatment_by_time. test shows different samples for same patient at specimen time 0.")
            for patient in patients:
                if specimen[specimen["PATIENT_ID"] == patient].drop_duplicates(subset=["SAMPLE_ID"]).shape[0] > 1:
                    print(patient, specimen[specimen["PATIENT_ID"] == patient].drop_duplicates(subset=["SAMPLE_ID"])["SAMPLE_ID"].tolist())

        specimen_times = append_df_with_mapping([mutations, specimen], sample_col, "SAMPLE_ID", "START_DATE",
                                                set_bar("mapping start times"), non_redundant=True, reverse=False)
        specimen_times = [int(specimen_time) if specimen_time != "-" else None for specimen_time in specimen_times]

        if "ID:specimen time" not in mutations.columns: mutations.insert(mutations.shape[1], "ID:specimen time", specimen_times)
        else:                                           mutations["ID:specimen time"] = specimen_times

        # check whether any values other than 0 exist (should fail)
        if mutations.drop_duplicates(subset=["ID:specimen time"]).shape[0] > 2:
            print("< inconsistent mapping @append_treatment_by_time")
            exit()

        # load treatment timeline to map against mutations (only treatments before or during sequencing)
        treatment        = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["treatment"], sep="\t")
        immune_treatment = treatment[treatment["START_DATE"] < 0]
        immune_treatment = immune_treatment[immune_treatment["SUBTYPE"] == "Immuno"]

        subtypes         = append_df_with_mapping([mutations, immune_treatment], patient_col, "PATIENT_ID", "SUBTYPE",
                                                  set_bar("mapping subtypes"), non_redundant=True, reverse=False)

        # marked (<-) added on 250527 to include treatment times
        treatment_starts = append_df_with_mapping([mutations, immune_treatment], patient_col, "PATIENT_ID", "START_DATE",
                                                  set_bar("mapping treatment starts"), non_redundant=True, reverse=False) # <- added
        
        # set 'None' if treatment data are not available
        patients = treatment["PATIENT_ID"].tolist()
        subtypes = [subtype if mutations.iloc[i].loc[patient_col] in patients else None for i, subtype in enumerate(subtypes)]
        if "ID:subtype" not in mutations.columns: mutations.insert(mutations.shape[1], "ID:subtype", subtypes)
        else:                                     mutations["ID:subtype"] = subtypes

        # marked (<-) added on 250527 to include treatment times
        treatment_starts = [int(treatment_start) if treatment_start != "-" else None for treatment_start in treatment_starts] # <- added
        if "ID:subtype" not in mutations.columns: mutations.insert(mutations.shape[1], "ID:treatment start", treatment_starts) # <- added
        else:                                     mutations["ID:treatment start"] = treatment_starts # <- added

        # check whether any values other than Immuno, placeholders "-" and "None" exist (should fail)
        if mutations.drop_duplicates(subset=["ID:subtype"]).shape[0] > 3:
            print("< inconsistent mapping @append_treatment_by_time")
            exit()

        # marked (<-) added on 250613
        mutations = pmu.append_markers(mutations, patient_col) # <- added
        
        # marked (<-) added on 250521 to allow appending features to processed data as input for survival analysis (can be merged later with append_treatment_by_time)
        # mutations.to_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"].split(".")[0]+"_appended.txt", sep=",", index=False) # <- removed
        mutations.to_csv(params["data_dir"]+params["os_sep"]+mutations_fname.split(".")[0]+"_appended.txt", sep=",", index=False) # <- added
        # pmu.test_ptcs(mutations, specimen, treatment) # <- removed
        pmu.test_ptcs(mutations, specimen, treatment, patient_col, sample_col) # <- added


    # <- if-statement added on 250616
    if params["mode"] == "filter_mutations":
        # load mutation data
        mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        # filter mutations to excluded redundant variants
        mutations = pmu.filter_mutations(mutations)

        # print to file
        mutations.to_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"].split(".")[0]+"_filtered.txt", sep=",", index=False)


    if params["mode"] == "map_expressions":
        # load mutation data
        mutations                       = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        # load expression data
        expression_data                 = pd.read_csv(params["expression_path"], sep=",")

        # assign basal expression to mutation data
        expressions                     = append_df_with_mapping([mutations, expression_data], "ID:gene symbol", "gene_name", params["expression_target"],
                                                                  set_bar("mapping expression data"), non_redundant=True, reverse=False)
        median_expression               = expression_data[params["expression_target"]].median()
        mutations["FEATURE:expression"] = [float(expression) if expression != "-" else median_expression for expression in expressions]
        mutations                       = mutations.to_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"].split(".")[0]+"_appended.txt", sep=",")


    if params["mode"] == "prepare_analysis" or params["mode"] == "test_patient_data":
        # load sample data
        sample_data = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["sample"], sep="\t")

        # load reference data (if reference_filter=True)
        if params["reference_filter"] == True:
            reference_sample_data = pd.read_csv(params["reference_path"], sep="\t")
            sample_data           = pmu.filter_samples(sample_data, reference_sample_data)

        else:
            # if sample ids with same patient id are present, only primary tumors are considered
            # if reference path is not None, only sample data not present in reference path are considered
            sample_data = pmu.filter_samples(sample_data)

        # load patient data
        patient_data = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["patient"], sep="\t")

        # map patient data to sample data
        sample_data = pmu.map_patient_data(sample_data, patient_data)

        # load raw mutation data
        mutations   = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["mutations"], sep="\t", low_memory=False)

        init_shape  = mutations.shape[0]
        mutations   = mutations[~mutations["HGVSp"].isna()]
        #mutations   = mutations[mutations["Tumor_Sample_Barcode"].isin(["P-0049359-T01-IM6", "P-0029588-T01-IM6", "P-0052105-T01-IM6"])]
        print("< removal of missing HGVSp value reduced data from", init_shape, "to", mutations.shape[0])

        # load ptc mutation data comtaining NMD score predictions
        mutations.insert(mutations.shape[1], "ID:variant id",
                         [mutations.iloc[i].loc["Transcript_ID"]+"_"+str(mutations.iloc[i].loc["Start_Position"])+"_"+str(mutations.iloc[i].loc["End_Position"])
                          for i in range(mutations.shape[0])])

        ptcs        = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        print("< format example of variant ids:", mutations["ID:variant id"].tolist()[0], "/", ptcs["ID:variant id"].tolist()[0])

        preds       = append_df_with_mapping([mutations, ptcs], "ID:variant id", "ID:variant id", "FEATURE:prediction",
                                              set_bar("mapping PTC data"), non_redundant=True, reverse=False)
        mutations.insert(mutations.shape[1], "FEATURE:prediction", [float(pred) if pred != "-" else None for pred in preds])

        # load expression data
        expression_data = pd.read_csv(params["expression_path"], sep=",")

        # assign basal expression to mutation data
        expressions             = append_df_with_mapping([mutations, expression_data], "Hugo_Symbol", "gene_name", params["expression_target"],
                                                          set_bar("mapping expression data"), non_redundant=True, reverse=False)
        median_expression       = expression_data[params["expression_target"]].median()
        mutations["expression"] = [float(expression) if expression != "-" else median_expression for expression in expressions]

        # prepare analysis
        if params["mode"] == "prepare_analysis":
            processed_patient_data = pmu.prepare_analysis(sample_data, mutations)

        if params["mode"] == "test_patient_data":
            processed_patient_data = pd.read_csv(params["data_dir"]+params["os_sep"]+params["outfname"], sep=",")

            # load treatment timeline to map against mutations (only treatments before or during sequencing)
            treatment = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["treatment"], sep="\t")
            #processed_patient_data = processed_patient_data[processed_patient_data["ID:SAMPLE_ID"].isin(["P-0049359-T01-IM6", "P-0029588-T01-IM6", "P-0052105-T01-IM6"])]

            pmu.test_patient_data(processed_patient_data, sample_data, mutations, treatment)


    if params["mode"] == "test_ptcs":
        # load mutation data
        mutations = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["ptcs"], sep=",")

        # load sequence timeline to exclude any samples with time of sequencing other than 0
        specimen  = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["specimen"], sep="\t")

        # load treatment timeline to map against mutations (only treatments before or during sequencing)
        treatment = pd.read_csv(params["data_dir"]+params["os_sep"]+params["fnames"]["treatment"], sep="\t")

        pmu.test_ptcs(mutations, specimen, treatment)


if __name__ == '__main__':
    main()