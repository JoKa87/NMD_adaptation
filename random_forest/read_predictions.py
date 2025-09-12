from read_predictions_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    params = {
             "create_newdir"                     : False,
             "data_dir"                          : parent_dir+r"\data",
             # list of errors specified in create_genome_predictions ("cds_mismatch" should better be name "lindeboom_cds_mismatch" as checking Lindeboom cds size)
             # ['cds_mismatch', 'chromosome_not_found', 'last_ejc_unknown',
             # 'lindeboom_cds_size_error', 'lindeboom_exon_index_error', 'lindeboom_exon_size_error', 'lindeboom_not_found', 'lindeboom_size_test_error',
             # 'total_exon_size_unknown']
             "error_filter"                      : ["chromosome_not_found", "total_exon_size_unknown"], 
             "os_sep"                            : "\\",
             "prediction_fnames"                 : ["hg38_NMD_susceptibilities.txt"],
             "prediction_transcript_identifier"  : "transcript id",
             "target_dir"                        : parent_dir+r"\data",
             "target_fname"                      : "tcga_variants.txt",
             "target"                            : {"position_identifier"       : "FEATURE:ptc cds position",
                                                    #"prediction_identifier"     : "FEATURE:prediction", # if identifier is provided, the read NMD score is compared to the input NMD score (for testing)
                                                    "transcript_identifier"     : "ID:transcript id"},
             "variant_id"                        : "ID:variant id", # "ID:variant id", should be None for Cuomo data to avoid wrong assignment of variants (with different transcript ids)
            }
    
    if params["variant_id"] == None:
        print("< warning. variant_id set to 'None' which is currently reserved for Cuomo data")

    print("<", params["target_fname"])
    # load predictions
    rpu = Read_predictions_utils(params)
    rpu.load(params["prediction_fnames"], mode="predictions")
    rpu.predictions[params["prediction_transcript_identifier"]] = [rpu.predictions.iloc[i].loc[params["prediction_transcript_identifier"]].split(".")[0]
                                                                   for i in range(rpu.predictions.shape[0])]

    
    # load targets
    if params["target_fname"] != None:
        targets = pd.read_csv(params["target_dir"]+params["os_sep"]+params["target_fname"], delimiter=",", index_col=False)

        # duplicates are removed so that mapping must be conducted afterwards (tools -> append_predictions)
        if params["variant_id"] != None and params["variant_id"] in targets.columns: targets = targets.drop_duplicates(subset=params["variant_id"])

        if params["target"]["transcript_identifier"] in targets.columns and params["target"]["position_identifier"] in targets.columns:
            targets = targets.sort_values(by=[params["target"]["transcript_identifier"]])
            targets = targets.reset_index()
            scores  = rpu.read_scores(targets, rpu.predictions)

            if "prediction_identifier" in params["target"] and params["target"]["prediction_identifier"] in targets.columns:
                rpu.compare_scores(scores, targets, rpu.predictions)

            rpu.print_scores(scores, targets)

        else:
            print("< specified identifiers", params["target"]["transcript_identifier"], "and/or", params["target"]["position_identifier"], "not found. read-out stopped.")


if __name__ == '__main__':
    main()