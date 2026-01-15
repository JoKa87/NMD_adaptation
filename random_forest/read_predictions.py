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
                                                    # if identifier is provided, the read NMD score is compared to the input NMD score (for testing)
                                                    #"prediction_identifier"     : "FEATURE:prediction", 
                                                    "transcript_identifier"     : "ID:transcript id"},
             # "ID:variant id", should be None for Cuomo data to avoid wrong assignment of variants (with different transcript ids)
             "variant_id"                        : "ID:variant id", 
            }
    
    if params["variant_id"] == None:
        print("< warning. variant_id set to 'None' which is currently reserved for Cuomo data")

    print("<", params["target_fname"])
    # load predictions
    rpu = Read_predictions_utils(params)
    rpu.load(params["prediction_fnames"], mode="predictions")
    rpu.predictions[params["prediction_transcript_identifier"]] = [rpu.predictions.iloc[i].loc[params["prediction_transcript_identifier"]].split(".")[0]
                                                                   for i in range(rpu.predictions.shape[0])]

    # comment regarding prediction of Cuomo data
    # mapping of ucids to transcript ids for genome build 19 is ambiguous
    # multiple occur due to mismatching CDS sizes, check on ensembl.org (250319) shows that in many cases, alternative transcript ids with matching cds can be found
    # examples (not complete) are shown below

    """
    variants with mismatching CDS size for which alternative transcript exists on ensembl.org (250319) with matching CDS size:
    ENST00000397356 -> ENST00000455446
    ENST00000393310 -> ENST00000409441
    ENST00000314759 -> ENST00000400018
    ENST00000425015 -> ENST00000334586
    ENST00000394773 -> ENST00000278845
    ENST00000539352 -> ENST00000252071
    ENST00000452576 -> ENST00000275358
    ENST00000243776 -> ENST00000535926
    ENST00000280330 -> ENST00000510047
    ENST00000463345 -> ENST00000542638
    ENST00000295981 -> ENST00000413608
    ENST00000334241 -> ENST00000591539
    ENST00000360916 -> ENST00000502374
    ENST00000392153 -> ENST00000589392
    ENST00000324411 -> ENST00000392153

    variants with mismatching CDS size for which alternative transcript do not exist on ensembl.org (250319) with matching CDS size:
    ENST00000391898, ENST00000434651, ENST00000236495, ENST00000537469, ENST00000262031, ENST00000252444

    yet to be determined:
    'ENST00000323813', 'ENST00000269829', 'ENST00000404249', 'ENST00000373652', 'ENST00000453321', 'ENST00000360586', 'ENST00000267430', 'ENST00000373345',
    'ENST00000530611', 'ENST00000282096', 'ENST00000455858', 'ENST00000407712', 'ENST00000359246', 'ENST00000420772', 'ENST00000458141', 'ENST00000525723',
    'ENST00000289528', 'ENST00000288985', 'ENST00000358495', 'ENST00000319622', 'ENST00000360422', 'ENST00000397902', 'ENST00000361847', 'ENST00000305233',
    'ENST00000355552', 'ENST00000549091', 'ENST00000372040', 'ENST00000301096', 'ENST00000536937', 'ENST00000541738', 'ENST00000290765', 'ENST00000359595',
    'ENST00000360004', 'ENST00000535617', 'ENST00000346416', 'ENST00000361936', 'ENST00000322764', 'ENST00000395023', 'ENST00000397293', 'ENST00000409509',
    'ENST00000409037', 'ENST00000274008', 'ENST00000263578', 'ENST00000486759', 'ENST00000307439', 'ENST00000444616', 'ENST00000415822', 'ENST00000264553',
    'ENST00000312989', 'ENST00000313843', 'ENST00000256594', 'ENST00000356683', 'ENST00000574993', 'ENST00000531224', 'ENST00000325083', 'ENST00000453153',
    'ENST00000393915', 'ENST00000368708', 'ENST00000441350', 'ENST00000441779', 'ENST00000277458'
    """
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