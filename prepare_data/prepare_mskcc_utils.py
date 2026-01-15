import os
import pandas as pd
import platform
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


class Prepare_mskcc_utils():
    def __init__(self, params):
        self.params      = params

        with open(self.params["data_dir"]+self.params["os_sep"]+self.params["outfname"].split(".")[0]+"_params.json", "w") as f:
            f.write(json.dumps(self.params, indent=4))


    # marked (<-) function added on 250527
    def _adjust_to_ptcs(self, adjusted_data, subset1, subset2, cancer_type, it):
        subset1_ptcs = self.create_ptcs(subset1)
        subset2_ptcs = self.create_ptcs(subset2)
        # remove mutations not present in either dataset
        # only required for subset1 as dataset is built later based on subset1 only (subset2-specific counts are not considered)
        for ptc_count in subset1_ptcs:
            if ptc_count not in subset2_ptcs: 
                subset1_ptcs[ptc_count] = [] # entries deleted so that they are ignored during reconstruction of the dataset

        adjusted_subset1 = []; adjusted_subset2 = []
        for ptc_count in subset1_ptcs:
            if len(subset1_ptcs[ptc_count]) > 0 and ptc_count not in subset2_ptcs:
                print("< unknown ptc_count", ptc_count, "detected @_adjust_to_ptcs")
                exit()

            elif len(subset1_ptcs[ptc_count]) > 0:
                #print(ptc_count, len(subset1_ptcs[ptc_count]), "/", len(subset2_ptcs[ptc_count]))
                if len(subset1_ptcs[ptc_count]) <= len(subset2_ptcs[ptc_count]):
                    adjusted_subset1, adjusted_subset2 = self.create_adjusted_subsets(adjusted_subset1, adjusted_subset2, subset1, subset2,
                                                                                      subset1_ptcs[ptc_count], subset2_ptcs[ptc_count])
                    
                elif len(subset1_ptcs[ptc_count]) > len(subset2_ptcs[ptc_count]):
                    adjusted_subset2, adjusted_subset1 = self.create_adjusted_subsets(adjusted_subset2, adjusted_subset1, subset2, subset1,
                                                                                      subset2_ptcs[ptc_count], subset1_ptcs[ptc_count])
        if len(adjusted_subset1) > 0 and len(adjusted_subset2) > 0:
            adjusted_subset1 = pd.concat(adjusted_subset1)
            adjusted_subset2 = pd.concat(adjusted_subset2)

            if adjusted_subset1.shape[0] != adjusted_subset2.shape[0]:
                print("< inconsistent sizes1 @_adjust_to_ptcs ("+str(adjusted_subset1.shape[0])+"/"+str(adjusted_subset2.shape[0])+")")
                exit()

            if adjusted_subset1.shape[0] > subset1.shape[0] and adjusted_subset2.shape[0] > subset2.shape[0]:
                print("< inconsistent sizes2 @_adjust_to_ptcs ("+str(adjusted_subset1.shape[0])+"/"+str(subset1.shape[0])+", "+str(adjusted_subset2.shape[0])+"/"+str(subset2.shape[0])+")")
                exit()

            show_plots = False
            if it == 0: show_plots = True
            adjusted_data = self.evaluate_adjusted_data(adjusted_data, adjusted_subset1, adjusted_subset2, cancer_type, it, show_plots)
        
        return adjusted_data


    # marked (<-) function added on 250527
    def adjust_to_ptcs(self, mutations):
        cancer_types  = np.unique(mutations["ID:cancer type"])
        adjusted_data = pd.DataFrame({**{"mw-statistic"+str(i): [None for _ in range(len(cancer_types))] for i in range(self.params["adjustment_params"]["repeats"])},
                                      **{"mw-p"+str(i):         [None for _ in range(len(cancer_types))] for i in range(self.params["adjustment_params"]["repeats"])},
                                      **{"mw-p_comb":           [None for _ in range(len(cancer_types))]},
                                      **{"ks-statistic"+str(i): [None for _ in range(len(cancer_types))] for i in range(self.params["adjustment_params"]["repeats"])},
                                      **{"ks-p"+str(i):         [None for _ in range(len(cancer_types))] for i in range(self.params["adjustment_params"]["repeats"])},
                                      **{"ks-p_comb":           [None for _ in range(len(cancer_types))]},
                                      **{"size1":               [None for _ in range(len(cancer_types))]},
                                      **{"size2":               [None for _ in range(len(cancer_types))]},
                                      **{"patients1":           [None for _ in range(len(cancer_types))]},
                                      **{"patients2":           [None for _ in range(len(cancer_types))]}},
                                      index=cancer_types)

        mutations     = mutations[~mutations["FEATURE:prediction"].isna()]
        #mutations     = mutations[mutations["FEATURE:expression"] >= self.params["adjustment_params"]["min_expression"]]

        bar = IncrementalBar(set_bar("adjusting to PTCs"), max=len(cancer_types)*self.params["adjustment_params"]["repeats"])
        for cancer_type in cancer_types:
            for i in range(self.params["adjustment_params"]["repeats"]):
                subset = mutations[mutations["ID:cancer type"] == cancer_type]
    
                if self.params["adjustment_params"]["mode"] == "immuno":
                    subset1 = subset[subset["ID:subtype"] == "Immuno"]
                    subset2 = subset[subset["ID:subtype"] == "-"]

                elif self.params["adjustment_params"]["mode"] == "metastasis":
                    subset1 = subset[subset["ID:sample type"] == "Metastasis"]
                    subset2 = subset[subset["ID:sample type"] == "Primary"]

                elif self.params["adjustment_params"]["mode"] in subset.columns:
                    subset1 = subset[subset[self.params["adjustment_params"]["mode"]] < 0] #subset[self.params["adjustment_params"]["mode"]].median()]
                    subset2 = subset[subset[self.params["adjustment_params"]["mode"]] > 0] # subset[self.params["adjustment_params"]["mode"]].median()]

                else:
                    print("< column not found")
                    exit()
                
                if subset1.shape[0] > 0 and subset2.shape[0] > 0:
                    adjusted_data = self._adjust_to_ptcs(adjusted_data, subset1, subset2, cancer_type, i)

                bar.next()

            pvalues = [adjusted_data.loc[cancer_type, "mw-p"+str(i)] for i in range(self.params["adjustment_params"]["repeats"])
                       if pd.isna(adjusted_data.loc[cancer_type, "mw-p"+str(i)]) == False]
            
            if len(pvalues) > 0: adjusted_data.at[cancer_type, "mw-p_comb"] = combine_pvalues(pvalues)

            pvalues = [adjusted_data.loc[cancer_type, "ks-p"+str(i)] for i in range(self.params["adjustment_params"]["repeats"])
                       if pd.isna(adjusted_data.loc[cancer_type, "ks-p"+str(i)]) == False]

            if len(pvalues) > 0: adjusted_data.at[cancer_type, "ks-p_comb"] = combine_pvalues(pvalues)

        bar.finish()
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        print(adjusted_data[[col for col in adjusted_data.columns if "ks-statistic" not in col]])


    # <- function added on 250613
    def append_markers(self, mutations, patient_col):
        mutations             = mutations.sort_values(by=patient_col).reset_index()
        extracted_marker_data = pd.DataFrame({**{"ID:"+marker+"_pre-seq" : [None for _ in mutations.index] for marker in self.params["fnames"]["markers"]},
                                              **{"ID:"+marker+"_post-seq": [None for _ in mutations.index] for marker in self.params["fnames"]["markers"]},
                                              **{"ID:"+marker+"_diff"    : [None for _ in mutations.index] for marker in self.params["fnames"]["markers"]}})
        block_index           = create_blocks(mutations, [patient_col])

        bar = IncrementalBar(set_bar("appending markers"), max=len(block_index)*len(self.params["fnames"]["markers"]))
        for marker in self.params["fnames"]["markers"]:
            # load marker file
            marker_data = pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["fnames"]["markers"][marker]["fname"], sep="\t")

            if marker_data[self.params["fnames"]["markers"][marker]["target_col"]].dtype == object:
                marker_data[self.params["fnames"]["markers"][marker]["target_col"]] = marker_data[self.params["fnames"]["markers"][marker]["target_col"]].replace({"No": 0, "Yes": 1})

            marker_data[self.params["fnames"]["markers"][marker]["target_col"]].astype("float32")

            testcount = 0
            # mapping marker file
            for i in range(len(block_index)):
                selected_marker_data = marker_data[marker_data["PATIENT_ID"] == mutations.iloc[block_index[i]["index"][0]].loc[patient_col]]
                pre_seq              = selected_marker_data[selected_marker_data["START_DATE"] <= 0][self.params["fnames"]["markers"][marker]["target_col"]].mean()
                post_seq             = selected_marker_data[selected_marker_data["START_DATE"] > 0][self.params["fnames"]["markers"][marker]["target_col"]].mean()

                if selected_marker_data.shape[0] > 0:
                    testcount += 1

                #if i < 10: print(marker, i, selected_marker_data[selected_marker_data["START_DATE"] <= 0].shape[0], selected_marker_data[selected_marker_data["START_DATE"] > 0].shape[0],
                #                 pre_seq, post_seq, len(block_index), mutations.iloc[block_index[i]["index"][0]].loc[patient_col])

                for j in block_index[i]["index"]:
                    extracted_marker_data.at[extracted_marker_data.index[j], "ID:"+marker+"_pre-seq"]  = pre_seq
                    extracted_marker_data.at[extracted_marker_data.index[j], "ID:"+marker+"_post-seq"] = post_seq
                    extracted_marker_data.at[extracted_marker_data.index[j], "ID:"+marker+"_diff"]     = post_seq-pre_seq

                bar.next()

            if patient_col == "PATIENT_ID" and marker_data.drop_duplicates(subset="PATIENT_ID").shape[0] != testcount:
                print("< dimension error occurred @append_markers ("+str(marker_data.drop_duplicates(subset="PATIENT_ID").shape[0])+"/"+str(testcount)+")")

        bar.finish()
        return pd.concat([mutations, extracted_marker_data], axis=1)


    def append_patient_data(self, mutations, patient_data, sample_data):
        if "ID:gender" not in mutations:
            genders = append_df_with_mapping([mutations, patient_data], "ID:patient id", "PATIENT_ID", self.params["gender_target"],
                                              set_bar("mapping patient data"), non_redundant=True, reverse=True)
            
            mutations.insert(mutations.shape[1], "ID:gender", [gender.upper() for gender in genders])
        
        sample_data["index"] = [i for i in range(sample_data.shape[0])]
        index                = append_df_with_mapping([mutations, sample_data], "ID:sample id", "SAMPLE_ID", "index",
                                                       set_bar("mapping sample data"), non_redundant=True, reverse=True)
        
        if "ID:cancer type" not in mutations: mutations.insert(mutations.shape[1], "ID:cancer type", [sample_data.iloc[int(i)].loc["CANCER_TYPE"] for i in index])
        else:                                 mutations["ID:cancer type"] = [sample_data.iloc[int(i)].loc["CANCER_TYPE"] for i in index]
        if "ID:sample type" not in mutations: mutations.insert(mutations.shape[1], "ID:sample type", [sample_data.iloc[int(i)].loc["SAMPLE_TYPE"] for i in index])
        else:                                 mutations["ID:sample type"] = [sample_data.iloc[int(i)].loc["SAMPLE_TYPE"] for i in index]
        if "ID:msi score" not in mutations:   mutations.insert(mutations.shape[1], "ID:msi score", [sample_data.iloc[int(i)].loc["MSI_SCORE"] for i in index])
        else:                                 mutations["ID:msi score"] = [sample_data.iloc[int(i)].loc["MSI_SCORE"] for i in index] # <- added on 250704
        if "ID:msi type" not in mutations:    mutations.insert(mutations.shape[1], "ID:msi type", [sample_data.iloc[int(i)].loc["MSI_TYPE"] for i in index])
        else:                                 mutations["ID:msi type"] = [sample_data.iloc[int(i)].loc["MSI_TYPE"] for i in index]  # <- added on 250704
        mutations.to_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["fnames"]["ptcs"].split(".")[0]+"_appended.txt", sep=",", index=False)


    # checks whether patient was subjected to any specified treatment, regardless of treatment time
    def append_treatment(self, mutations, treatment, target_col="ID:PATIENT_ID"):
        subtypes = treatment.drop_duplicates(subset="SUBTYPE")["SUBTYPE"].tolist()
        subtypes = {"ID:"+subtype: [] for subtype in subtypes}
        print(subtypes)

        bar = IncrementalBar(set_bar("appending treatment data"), max=mutations.shape[0])
        for i in range(mutations.shape[0]):
            selected_treatment = treatment[treatment["PATIENT_ID"] == mutations.iloc[i].loc[target_col]]

            for key in subtypes:
                if key.replace("ID:", "") in selected_treatment["SUBTYPE"].tolist(): subtypes[key].append(1)
                else:                                                                subtypes[key].append(0)

            bar.next()
        bar.finish()

        for key in subtypes:
            mutations.insert(mutations.shape[1], key, subtypes[key])

        return mutations


    def count_by_threshold(self, scores, score_threshold):
        return len([score for score in scores if score >= score_threshold[0] and score <= score_threshold[1]])


    # marked (<-) function added on 250527
    def create_ptcs(self, mutations):
        ptcs = {}
        patients = np.unique(mutations["ID:patient id"])
        for patient in patients:
            ptc_count = mutations[mutations["ID:patient id"] == patient].shape[0]
            if ptc_count in ptcs: ptcs[ptc_count].append(patient)
            else:                 ptcs[ptc_count] = [patient]
        
        return ptcs


    def create_adjusted_subsets(self, adjusted_subset1, adjusted_subset2, subset1, subset2, ptc_counts1, ptc_counts2):
        rand_index = np.arange(len(ptc_counts2))
        np.random.shuffle(rand_index)
        rand_index = rand_index.tolist()

        for i in range(len(ptc_counts1)):
            adjusted_subset1.append(subset1[subset1["ID:patient id"] == ptc_counts1[i]])
            adjusted_subset2.append(subset2[subset2["ID:patient id"] == ptc_counts2[rand_index[i]]])

        return adjusted_subset1, adjusted_subset2


    # marked (<-) function added on 250527
    def evaluate_adjusted_data(self, adjusted_data, subset1, subset2, cancer_type, it, show_plots=False):
        preds1 = subset1["FEATURE:prediction"].tolist()
        preds2 = subset2["FEATURE:prediction"].tolist()

        mw_statistic = scipy.stats.mannwhitneyu(preds1, preds2).statistic
        mw_statistic = (len(preds1)*len(preds2) - mw_statistic) / mw_statistic
        mw_p         = scipy.stats.mannwhitneyu(preds1, preds2).pvalue
        ks_statistic = scipy.stats.kstest(preds1, preds2).statistic
        ks_p         = scipy.stats.kstest(preds1, preds2).pvalue

        adjusted_data.at[cancer_type, "mw-statistic"+str(it)] = mw_statistic
        adjusted_data.at[cancer_type, "mw-p"+str(it)]         = mw_p
        adjusted_data.at[cancer_type, "ks-statistic"+str(it)] = ks_statistic
        adjusted_data.at[cancer_type, "ks-p"+str(it)]         = ks_p
        adjusted_data.at[cancer_type, "size1"]                = len(preds1)
        adjusted_data.at[cancer_type, "size2"]                = len(preds2)
        adjusted_data.at[cancer_type, "patients1"]            = len(np.unique(subset1["ID:patient id"]))
        adjusted_data.at[cancer_type, "patients2"]            = len(np.unique(subset2["ID:patient id"]))

        if show_plots == True:
            plt.title(cancer_type+" ("+str(it)+")")
            plt.hist(preds1, bins=30, histtype="step", density=True,
                    label="subset1 ("+str(round(mw_statistic, 4))+" / "+str(round(mw_p, 4))+", "+str(round(ks_statistic, 4))+" / "+str(round(ks_p, 4))+")")
            plt.hist(preds2, bins=30, histtype="step", density=True, label="subset2")
            plt.legend()
            plt.show()

        return adjusted_data


    # <- function added on 250616 to filter out downstream PTC mutations of PTC variants with existing mutations
    # <- changed on 251114 from ID:patient id to ID:sample id as samples always reflect different cancer types and can thus be considered independent
    def filter_mutations(self, mutations):
        mutations = mutations.sort_values(["ID:sample id", "ID:transcript id", "FEATURE:ptc cds position"])
        blocks    = create_blocks(mutations, ["ID:sample id", "ID:transcript id"])
        #blocks    = create_blocks(mutations, ["ID:patient id", "ID:transcript id", "FEATURE:ptc cds position"])

        filtered_mutations = {col: [] for col in mutations.columns}
        bar = IncrementalBar(set_bar("filtering mutations"), max=len(blocks))
        for i in range(len(blocks)):
            if mutations.iloc[blocks[i]["index"][0]].loc["FEATURE:ptc cds position"] == mutations.iloc[blocks[i]["index"]]["FEATURE:ptc cds position"].min():
                for col in mutations.columns:
                    filtered_mutations[col].append(mutations.iloc[blocks[i]["index"][0]].loc[col])

            else:
                print("< error occurred @filter_mutations. selected PTC cds position is not most upstream")
                print(mutations.iloc[blocks[i]["index"][0]].loc["FEATURE:ptc cds position"], "/", mutations.iloc[blocks[i]["index"]].loc["FEATURE:ptc cds position"].tolist())
                exit()

            if len(blocks[i]["index"]) > 1:
                tests = []
                for j in range(len(blocks[i]["index"])):
                    tests.append((mutations.iloc[blocks[i]["index"][j]].loc["ID:sample id"]
                                 +"_"+mutations.iloc[blocks[i]["index"][j]].loc["ID:transcript id"]))
                                 #+"_"+str(mutations.iloc[blocks[i]["index"][j]].loc["FEATURE:ptc cds position"])))

                if np.unique(tests).shape[0] != 1:
                    print("< dimension error occurred @filter_mutations")
                    print(tests)
                    exit()
            
            bar.next()
        bar.finish()

        filtered_mutations = pd.DataFrame(filtered_mutations)

        if filtered_mutations.drop_duplicates(subset=["ID:sample id", "ID:transcript id"]).shape[0] != filtered_mutations.shape[0]:
            print("< error occurred.", filtered_mutations.shape[0]-filtered_mutations.drop_duplicates(subset=["ID:sample id", "ID:transcript id"]).shape[0],
                  "downstream ptcs could not be removed.")
            exit()

        else:
            print("< filtering removed mutations from", mutations.shape[0], "to", filtered_mutations.shape[0])

        return filtered_mutations


    def filter_samples(self, sample_data):
        sample_data = sample_data.sort_values(by=["PATIENT_ID"])
        sample_data = sample_data.reset_index()
        blocks      = create_blocks(sample_data, ["PATIENT_ID"])

        multiple_primaries = []; selected_index = []
        for block in blocks:
            if len(block["index"]) > 1:
                primary_data = sample_data.iloc[block["index"]][sample_data.iloc[block["index"]]["SAMPLE_TYPE"] == "Primary"]

                if primary_data.shape[0] > 1:  multiple_primaries.append(block["block id"][list(block["block id"].keys())[0]])
                if primary_data.shape[0] >= 1: selected_index.append(primary_data.index[0])
                else:                          selected_index.append(block["index"][0])

            else:
                selected_index.append(block["index"][0])

            # check validity of selected index
            if sample_data.iloc[block["index"][0]].loc["PATIENT_ID"] != block["block id"][list(block["block id"].keys())[0]]:
                print("< invalid index @filter_samples:", sample_data.iloc[block["index"][0]].loc["PATIENT_ID"], "/", block["block id"][list(block["block id"].keys())[0]])
                exit()

        if len(multiple_primaries) > 0:
            print("< multiple primary data detected for", len(multiple_primaries), "patients. only first entries were selected:")
            print(multiple_primaries)

        if sample_data.iloc[selected_index].drop_duplicates(subset="PATIENT_ID").shape[0] != len(selected_index):
            print("< wrong filtering @filter_samples:", sample_data.iloc[selected_index].drop_duplicates(subset="PATIENT_ID").shape[0], "/", len(selected_index))
            exit()

        print("< redundancy filtering reduced sample data from", sample_data.shape[0], "to", len(selected_index))
        return sample_data.iloc[selected_index]


    def map_patient_data(self, sample_data, patient_data):
        patient_data["index"] = [i for i in range(patient_data.shape[0])]
        mapped_index          = append_df_with_mapping([sample_data, patient_data], "PATIENT_ID", "PATIENT_ID", "index",
                                                        set_bar("mapping patient data"), non_redundant=True, reverse=False)
        
        for label_id in self.params["label_ids"]:
            sample_data.insert(sample_data.shape[1], label_id, [patient_data.iloc[int(i)].loc[label_id] if i != "-" else None for i in mapped_index])
        
        # convert "OS_STATUS" to 0/1
        sample_data["OS_STATUS"] = [int(sample_data.iloc[i].loc["OS_STATUS"].split(":")[0]) for i in range(sample_data.shape[0])]
        return sample_data


    # <- alternative function created on 251028 (_prepare_analysis_from_ptcs)
    def _prepare_analysis(self, sample_data, mutations):
        # rename columns
        updated_cols = {}
        for col in sample_data.columns:
            if col in self.params["label_ids"]:         updated_cols[col] = "LABEL:"+col
            elif col in self.params["feature_targets"]: updated_cols[col] = "FEATURE:"+col
            else:                                       updated_cols[col] = "ID:"+col

        sample_data = sample_data.rename(columns=updated_cols)

        # add categories
        mutations["index"]  = [i for i in range(mutations.shape[0])]
        mapped_index        = append_df_with_mapping([sample_data, mutations], "ID:SAMPLE_ID", "Tumor_Sample_Barcode", "index",
                                                      set_bar("mapping MSK data"), non_redundant=False, reverse=True)
        
        extracted_mutations = {**{"FEATURE:expression":             [None for _ in mapped_index]},
                               **{"FEATURE:frameshift":             [None for _ in mapped_index]},
                               **{"FEATURE:frameshift_mutations":   [None for _ in mapped_index]},
                               **{"FEATURE:escape":                 [None for _ in mapped_index]},
                               **{"FEATURE:prediction":             [None for _ in mapped_index]},
                               **{"FEATURE:target":                 [None for _ in mapped_index]},
                               **{"FEATURE:ptc_mutations":          [None for _ in mapped_index]}}

        bar = IncrementalBar(set_bar("processing mutation data"), max=len(mapped_index))
        for i in range(len(mapped_index)):
            for j in range(len(mapped_index[i])):
                index = int(mapped_index[i][j])

            selected_mutations = mutations.iloc[[int(j) for j in mapped_index[i]]]
            selected_mutations = selected_mutations[selected_mutations["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"])]

            if selected_mutations.shape[0] > 0:
                extracted_mutations["FEATURE:expression"][i] = selected_mutations["expression"].mean()
                # frameshift
                values = [get_frameshift(selected_mutations.iloc[j].loc["HGVSp"]) for j in range(selected_mutations.shape[0])
                          if type(selected_mutations.iloc[j].loc["HGVSp"]) == str and "?" not in selected_mutations.iloc[j].loc["HGVSp"]
                             and len(selected_mutations.iloc[j].loc["HGVSp"].split("Ter")) == 2]

                # remove None values
                values = [value for value in values if pd.isna(value) == False]
                if len(values) > 0: extracted_mutations["FEATURE:frameshift"][i] = np.mean(values)
                
                values = [value for value in values if value > 0]
                extracted_mutations["FEATURE:frameshift_mutations"][i] = len(values)
                
                # predictions
                scores                                       = [selected_mutations.iloc[j].loc["FEATURE:prediction"] for j in range(selected_mutations.shape[0])
                                                                if pd.isna(selected_mutations.iloc[j].loc["FEATURE:prediction"]) == False]

                if len(scores) > 0:
                    extracted_mutations["FEATURE:escape"][i]     = self.count_by_threshold(scores, self.params["score_threshold"]["FEATURE:escape"]) / len(scores)
                    extracted_mutations["FEATURE:prediction"][i] = np.mean(scores)
                    extracted_mutations["FEATURE:target"][i]     = self.count_by_threshold(scores, self.params["score_threshold"]["FEATURE:target"]) / len(scores)

                elif self.params["zeros_included"] == True:
                    extracted_mutations["FEATURE:escape"][i]     = 0
                    extracted_mutations["FEATURE:target"][i]     = 0

                
            elif self.params["zeros_included"] == False:
                # default values should be consistent with settings for analyze_targets
                extracted_mutations["FEATURE:frameshift_mutations"][i] = 0

            elif self.params["zeros_included"] == True:
                extracted_mutations["FEATURE:escape"][i]               = 0
                extracted_mutations["FEATURE:expression"][i]           = 0
                extracted_mutations["FEATURE:frameshift"][i]           = 0
                extracted_mutations["FEATURE:frameshift_mutations"][i] = 0
                extracted_mutations["FEATURE:target"][i]               = 0

            extracted_mutations["FEATURE:ptc_mutations"][i] = selected_mutations.shape[0]
            
            bar.next()
        bar.finish()

        for feature in extracted_mutations:
            if None not in extracted_mutations[feature]: print(feature, np.mean(extracted_mutations[feature]))
            sample_data.insert(sample_data.shape[1], feature, extracted_mutations[feature])

        print("< total values", sample_data.shape[0],
              ", existing frameshifts", sample_data[~sample_data["FEATURE:frameshift"].isna()].shape[0],
              ", existing prediction values", sample_data[~sample_data["FEATURE:prediction"].isna()].shape[0])

        return sample_data
    

    def _prepare_analysis_from_ptcs(self, sample_data, ptcs):
        # rename columns
        updated_cols = {}
        for col in sample_data.columns:
            if col in self.params["label_ids"]:         updated_cols[col] = "LABEL:"+col
            elif col in self.params["feature_targets"]: updated_cols[col] = "FEATURE:"+col
            else:                                       updated_cols[col] = "ID:"+col

        sample_data = sample_data.rename(columns=updated_cols)

        # add categories
        ptcs["index"] = [i for i in range(ptcs.shape[0])]
        mapped_index  = append_df_with_mapping([sample_data, ptcs], "ID:SAMPLE_ID", "ID:sample id", "index",
                                                set_bar("mapping MSK data"), non_redundant=False, reverse=True)
        
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
                    extracted_mutations["FEATURE:escape"][i]     = self.count_by_threshold(scores, self.params["score_threshold"]["FEATURE:escape"]) / len(scores)
                    extracted_mutations["FEATURE:prediction"][i] = np.mean(scores)
                    extracted_mutations["FEATURE:target"][i]     = self.count_by_threshold(scores, self.params["score_threshold"]["FEATURE:target"]) / len(scores)

                elif self.params["zeros_included"] == True:
                    extracted_mutations["FEATURE:escape"][i]     = 0
                    extracted_mutations["FEATURE:target"][i]     = 0

                
            elif self.params["zeros_included"] == False:
                # default values should be consistent with settings for analyze_targets
                extracted_mutations["FEATURE:frameshift_mutations"][i] = 0

            elif self.params["zeros_included"] == True:
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
            sample_data.insert(sample_data.shape[1], feature, extracted_mutations[feature])

        print("< total values", sample_data.shape[0],
              ", existing frameshifts", sample_data[~sample_data["FEATURE:frameshift"].isna()].shape[0],
              ", existing prediction values", sample_data[~sample_data["FEATURE:prediction"].isna()].shape[0])

        return sample_data


    # preparation of data from folder: tmb_mskcc_2018, files: data_clinical_sample_edited.txt, data_mutations_appended.txt, and TCGA_expressions.txt
    def prepare_analysis(self, sample_data, ptcs):
        if self.params["cancer_wise"] == True:
            cancer_types     = sample_data.drop_duplicates(subset="CANCER_TYPE")["CANCER_TYPE"].tolist()
            full_sample_data = sample_data

            for cancer_type in cancer_types:
                sample_data           = full_sample_data[full_sample_data["CANCER_TYPE"] == cancer_type]
                processed_sample_data = self._prepare_analysis_from_ptcs(full_sample_data, ptcs)
                processed_sample_data = get_random_cohorts(processed_sample_data, self.params["cohorts"], targets=["ID:PATIENT_ID"])
                processed_sample_data.to_csv(path_or_buf=self.params["data_dir"]+self.params["os_sep"]+self.params["outfname"].split(".")[0]+"_"+cancer_type.lower().replace(" ", "_")+".txt", sep=",", index=False)

        else:
            processed_sample_data = self._prepare_analysis_from_ptcs(sample_data, ptcs)
            processed_sample_data = get_random_cohorts(processed_sample_data, self.params["cohorts"], targets=["ID:PATIENT_ID"])
            processed_sample_data.to_csv(path_or_buf=self.params["data_dir"]+self.params["os_sep"]+self.params["outfname"], sep=",", index=False)

        return processed_sample_data


    def test_patient_data(self, processed_patient_data, sample_data, mutations, treatment):
        features       = ["escape", "expression", "frameshift", "frameshift_mutations", "ptc_mutations", "target"]
        treatment_cols = ["ID:Chemo", "ID:Investigational", "ID:Immuno", "ID:Hormone", "ID:Bone Treatment", "ID:Biologic", "ID:Targeted", "ID:Other"]

        bar = IncrementalBar(set_bar("testing patient data"), max=processed_patient_data.shape[0])
        for i in range(processed_patient_data.shape[0]):
            selected_mutations   = mutations[mutations["Tumor_Sample_Barcode"] == processed_patient_data.iloc[i].loc["ID:SAMPLE_ID"]]
            selected_mutations   = selected_mutations[selected_mutations["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"])]
            selected_sample_data = sample_data[sample_data["PATIENT_ID"] == processed_patient_data.iloc[i].loc["ID:PATIENT_ID"]]
            selected_treatment   = treatment[treatment["PATIENT_ID"] == processed_patient_data.iloc[i].loc["ID:PATIENT_ID"]]

            if selected_sample_data.shape[0] == 0:
                print("< patient not found @", i, processed_patient_data.iloc[i].loc["PATIENT_ID"])

            # check validity of filtering
            if "Primary" in selected_sample_data["SAMPLE_TYPE"].tolist() and processed_patient_data.iloc[i].loc["ID:SAMPLE_TYPE"] != "Primary":
                print("< filtering error occurred @", i, processed_patient_data.iloc[i].loc["PATIENT_ID"])
            
            if "Primary" in selected_sample_data["SAMPLE_TYPE"].tolist():
                filtered_sample_data = selected_sample_data[selected_sample_data["SAMPLE_TYPE"] == "Primary"]

            else:
                filtered_sample_data = selected_sample_data
            
            if i < 10: print("i", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"], selected_mutations.shape[0], selected_sample_data.shape[0],
                       selected_treatment.shape[0], filtered_sample_data.shape[0])

            # check cancer type
            if processed_patient_data.iloc[i].loc["ID:CANCER_TYPE"] != filtered_sample_data.iloc[0].loc["CANCER_TYPE"]:
                print("< cancer type error occurred @test_patient_data", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"],
                      processed_patient_data.iloc[i].loc["ID:CANCER_TYPE"], "/", filtered_sample_data.iloc[0].loc["CANCER_TYPE"])

            # check sample type
            if processed_patient_data.iloc[i].loc["ID:SAMPLE_TYPE"] != filtered_sample_data.iloc[0].loc["SAMPLE_TYPE"]:
                print("< cancer type error occurred @test_patient_data", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"],
                      processed_patient_data.iloc[i].loc["ID:SAMPLE_TYPE"], "/", filtered_sample_data.iloc[0].loc["SAMPLE_TYPE"])

            # check treatment features
            for treatment_col in treatment_cols:
                if processed_patient_data.iloc[i].loc[treatment_col] == 0 and treatment_col.replace("ID:", "") in selected_treatment["SUBTYPE"].tolist():
                    print("< treatment assignment error1 occurred @", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"], processed_patient_data.iloc[i].loc[treatment_col])

                if processed_patient_data.iloc[i].loc[treatment_col] == 1 and treatment_col.replace("ID:", "") not in selected_treatment["SUBTYPE"].tolist():
                    print("< treatment assignment error2 occurred @", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"], processed_patient_data.iloc[i].loc[treatment_col])

            # check features
            feature_values                  = {feature: None for feature in features}
            feature_values["ptc_mutations"] = selected_mutations.shape[0]

            if selected_mutations.shape[0] == 0 and self.params["zeros_included"] == True:
                for key in feature_values:
                    if key != "prediction":
                        feature_values[key] = 0

            if selected_mutations.shape[0] > 0:
                feature_values["expression"]           = selected_mutations["expression"].mean()
                values                                 = [get_frameshift(selected_mutations.iloc[j].loc["HGVSp"]) for j in range(selected_mutations.shape[0])
                                                          if type(selected_mutations.iloc[j].loc["HGVSp"]) == str and "?" not in selected_mutations.iloc[j].loc["HGVSp"]
                                                          and len(selected_mutations.iloc[j].loc["HGVSp"].split("Ter")) == 2]
                values                                 = [value for value in values if pd.isna(value) == False]
                feature_values["frameshift"]           = np.mean(values)
                feature_values["frameshift_mutations"] = len([value for value in values if value > 0])

                scores                                 = [selected_mutations.iloc[j].loc["FEATURE:prediction"] for j in range(selected_mutations.shape[0])
                                                          if pd.isna(selected_mutations.iloc[j].loc["FEATURE:prediction"]) == False]

                if len(scores) > 0:
                    feature_values["escape"]     = len([score for score in scores if score <= self.params["score_threshold"]["FEATURE:escape"][1]])/len(scores)
                    feature_values["prediction"] = selected_mutations["FEATURE:prediction"].mean()
                    feature_values["target"]     = len([score for score in scores if score >= self.params["score_threshold"]["FEATURE:target"][0]])/len(scores)

                elif self.params["zeros_included"] == True:
                    feature_values["escape"]     = 0
                    feature_values["target"]     = 0


            for key in feature_values:
                if pd.isna(processed_patient_data.iloc[i].loc["FEATURE:"+key]) == True and pd.isna(feature_values[key]) == True:
                    pass

                elif round(processed_patient_data.iloc[i].loc["FEATURE:"+key], 5) != round(feature_values[key], 5):
                    print("<", key, "mismatch @", i, processed_patient_data.iloc[i].loc["ID:PATIENT_ID"], processed_patient_data.iloc[i].loc["FEATURE:"+key], "/", feature_values[key])
                
            bar.next()
        bar.finish()


    def test_ptcs(self, mutations, specimen, treatment, patient_col, sample_col):
        bar = IncrementalBar(set_bar("testing PTC data"), max=mutations.shape[0])
        for i in range(mutations.shape[0]):
            selected_specimen = specimen[specimen["SAMPLE_ID"] == mutations.iloc[i].loc[sample_col]]
            selected_specimen = selected_specimen[selected_specimen["START_DATE"] == 0]

            if selected_specimen.shape[0] > 0:
                selected_treatment = treatment[treatment["PATIENT_ID"] == mutations.iloc[i].loc[patient_col]]
                selected_treatment = selected_treatment[selected_treatment["START_DATE"] < 0]

                if selected_treatment.shape[0] > 0:
                    if "Immuno" in selected_treatment["SUBTYPE"].tolist() and mutations.iloc[i].loc["ID:subtype"] != "Immuno":
                        print("< error1 occurred @test_ptcs", i, mutations.iloc[i].loc[sample_col])

                    if "Immuno" not in selected_treatment["SUBTYPE"].tolist() and mutations.iloc[i].loc["ID:subtype"] == "Immuno":
                        print("< error2 occurred @test_ptcs", i, mutations.iloc[i].loc[sample_col])

            if pd.isna(mutations.iloc[i].loc["ID:subtype"]) == False and mutations.iloc[i].loc[patient_col] not in treatment["PATIENT_ID"].tolist():
                print("< error3 occurred @", i, mutations.iloc[i].loc[sample_col])

            if pd.isna(mutations.iloc[i].loc["ID:subtype"]) == True and mutations.iloc[i].loc[patient_col] in treatment["PATIENT_ID"].tolist():
                print("< error4 occurred @", i, mutations.iloc[i].loc[sample_col])

            bar.next()
        bar.finish()