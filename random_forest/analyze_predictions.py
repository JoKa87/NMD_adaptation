import copy
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\plots")

from analyze_predictions_utils import *
from plot_utils import *
from profiles import *


# parameters
params = {
         "analysis_steps"               : 2,
         "apply_3mers"                  : True,
         "apply_cutoff"                 : False, # if True NMD escape and target cutoffs are applied for selection analysis, if False internal reference (avg. model prediction) is used
         "apply_global_predictions"     : False, # if True, randomized subset of global predictions is used for testing, "apply_cutoff" auto-set to True nmd_escape_cutoff and nmd_target_cutoff set to subset average
         "apply_mutation_stats"         : True,
         "apply_mutation_weights"       : False,
         "appris_priorities"            : ["PRINCIPAL:1", "PRINCIPAL:2", "PRINCIPAL:3", "PRINCIPAL:4", "PRINCIPAL:5", "ALTERNATIVE:1", "ALTERNATIVE:2", "MINOR"],
         "appris_selection"             : False,
         "average_indel_probabilities"  : True,
         "block"                        : "ID:sample id", # "ID:cancer type" "ID:project" None
         "block_targets"                : ["ID:project", "ID:case id"],
         "bootstrapping_steps"          : 100000,
         "create_newdir"                : True,
         "data_dir"                     : parent_dir+r"\data",
         "error_filter"                 : ["chromosome_not_found", "total_exon_size_unknown"], #["chromosome_not_found", "total_exon_size_unknown"], 
         "exon_resolution"              : True,
         "file_tag"                     : "control_projectwise",
         "filters"                      : [], #
         "global_dir"                   : parent_dir+r"\data",
         "global_predictions_path"      : r"", # fill in path to pre-calculated stats file
         "harmonize"                    : True, # if True predictions and variants are harmonized based on variant_filter_col
         "hist_steps"                   : 100,
         "id_filter"                    : "custom", # "CPTAC3" "Cuomo" "MSK" "TCGA" "Teran" 
         "knowngene_fname"              : "hg38_knownGene_appended.txt",
         "load_lindeboom"               : False,
         "lindeboom_fnames"             : ["hg38_lindeboom_predictions.txt"], # "h38_lindeboom_predictions.txt" "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"
         "mapping_target"               : "transcript id", # "gene id" "transcript id"
         "masks"                        : ["frameshift", "nonsense"], # "frameshift" "nonsense" [] <- added on 250528
         # "analyze_adaptation" "analyze_blocks" "analyze_predictions" "analyze_selection" "analyze_stop_codons" "convert_lindeboom_predictions" "show_stats"
         "min_patients"                 : 10, # added on 251024
         "mode"                         : "analyze_predictions",
         "mutation_stats_del_scale"     : 1, # <- added on 250614 to adjust indel probabilities to ratio of observed to projected frameshift mutations with scale 1
         "mutation_stats_fname"         : "mutation_stats.json", # <- added on 250603
         "mutation_stats_gene_target"   : "total", # <- added on 250603, in project-wise mode, None forces values to be replaced by current block
         "mutation_stats_ins_scale"     : 1, # <- added on 250614 to adjust indel probabilities to ratio of observed to projected frameshift mutations with scale 1
         "mutation_stats_non_scale"     : 1, # <- added on 250624 to adjust indel probabilities to ratio of observed to projected nonsense mutations with scale 1
         "mutation_stats_pair_target"   : "total", # <- added on 250603, in project-wise mode, None forces values to be replaced by current block
         "mutation_stats_ptc_weights"   : False, # <- added on 250821, allows application of PTC weights during masking (overwrites empirical gene factors)
         "mutation_stats_scale"         : 1, 
         "nmd_escape_cutoff"            : 0.57,
         "nmd_target_cutoff"            : 0.64,
         "os_sep"                       : "\\",
         "plot_projections"             : True,
         "polynomial_degree"            : 9,
         "prediction_fnames"            : ["hg38_NMD_scores_appended_appended.txt"], # "hg38_NMD_scores_appended_appended.txt" contains exon correction
         "profile"                      : "tcga", # if None this field is ignored <- added on 250709
         "projection_cutoff"            : 20,
         "random_sample_size"           : None, # if None this field is ignored <- added on 251030
         "remove_inconsistencies"       : True, # <- added on 250704
         "transformation_steps"         : 100,
         "selection_cols"               : {"gene symbol": "ID:gene symbol"}, #"project": "ID:project", "case id": "ID:case id", "sample id": "ID:sample id"},
         #"selection_cols"               : {"gene symbol": "ID:gene symbol", "sample id": "ID:sample id", "cancer type": "ID:cancer type"},
         "selection_fname"              : "tcga_frameshift_projectwise",
         "selected_genes"               : ["ENSG00000134086"],
         "selection_method"             : "binomial", # "binomial" "bootstrapping" "fishers exact" "fishers exact escape" "fishers exact target"
         "selection_minsize"            : 10, # <- added on 251021 to allow exclusion of low-counting variants, "None" to switch off
         "selection_mode"               : "relative", # "absolute" "relative"
         "selection_output_type"        : "GOrilla", # output type: "gProfiler" (http://biit.cs.ut.ee/gprofiler/gost), "GOrilla" (https://cbl-gorilla.cs.technion.ac.il/)
         "selected_transcripts"         : ["ENST00000269305.9", "ENST00000257430.9"], # ENST00000256474.3
         "seq_fname"                    : "hg38_seqs_appended.txt", # <- added on 251111
         "shared_variants"              : False,
         "shared_variants_cases"        : False,
         "shared_variants_filter"       : 2,
         "stats_mode"                   : "values", # "means" "values"
         "status_fname"                 : "filtered_status_ptc.json",
         "tcga_filter"                  : "biallelic_only", # "biallelic_only" "off"
         "use_variant_filter"           : False, # variant filter specified in variant_filter_fnames is loaded and applied to predictions and variants
         "variant_filter_col"           : "ID:transcript id",
         "variant_filter_fnames"        : ["tcga.txt"],
         # if "_EXCLUDED" is added to key, all values but indicated one are included 
         "variant_filter_value_filter"  : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False]}, # "FEATURE:frameshift_IDENTICAL": 0}, # {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0]}, #}, # "ID:cancer type": ["Non-Small Cell Lung Cancer"]
         "variant_fnames"               : ["control_filtered_appended.txt"],
         "variants_only"                : False, # if True only variant data are projected
         "variant_targets"              : {
                                           "FEATURE:prediction"              : "FEATURE:abs. ptc index",
                                           "FEATURE:ptc cds position"        : "FEATURE:abs. ptc index",
                                          },
         "variant_test_cutoff"          : 1,
         "variant_test_targets"         : ["FEATURE:prediction"], # features used to test selection (i.e. different distribution of real variants vs. theoretical variants)
         "variant_value_filter"         : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Frame_Shift_Del", "Frame_Shift_Ins"]},# "FEATURE:frameshift_IDENTICAL": 0} # , "FEATURE:frameshift_IDENTICAL": 0, "ID:variant classification": ["Nonsense_Mutation"]}, # minimum expression for variants to be included in the analysis "ID:ptc reads": 0
         }


if params["profile"] != None:
    params = apply_profile(params)

if params["mode"] == "analyze_predictions" and params["variant_test_cutoff"] != 1:
    params["variant_test_cutoff"] = 1
    print("< for analyze_predictions mode, significance cutoff is set to 1")

if params["load_lindeboom"] == True and len(params["variant_targets"]) > 0:
    print("< variant targets are incompatible with lindeboom options due to deviating absolute positions.")
    exit()

if len(params["masks"]) == 0 and params["apply_mutation_stats"] == True:
    print("< mask must be specified in order to apply mutation stats")
    exit()

# initialize class
print(params["file_tag"])

if (params["mode"] == "analyze_predictions" or params["mode"] == "analyze_adaptation") and (params["apply_mutation_stats"] == True or params["apply_mutation_weights"] == True): # recheck
    with open(params["data_dir"]+params["os_sep"]+params["mutation_stats_fname"], "r") as f:
        mutation_stats = json.load(f)

else:
    mutation_stats = {}

apu = Analyze_predictions_utils(params, mutation_stats)

if params["apply_global_predictions"] == True:
    params["apply_cutoff"]      = True
    apu.load([params["global_predictions_path"]], mode="global_predictions")
    params["nmd_escape_cutoff"] = np.mean(apu.global_predictions)
    params["nmd_cutoff"]        = np.mean(apu.global_predictions)
    params["nmd_target_cutoff"] = np.mean(apu.global_predictions)
    print("< mean", np.mean(apu.global_predictions), len(apu.global_predictions))


def _adjust_probabilities(mutation_stats, target):
    # apply empirical correction to adjust frameshift probabilities
    for key1 in ["del-1", "del-1_3mers"]:
        if target == "total":
            for key2 in mutation_stats[key1]:
                for key3 in mutation_stats[key1][key2]:
                    mutation_stats[key1][key2][key3] *= params["mutation_stats_del_scale"]

        else:
            for key3 in mutation_stats[key1][target]:
                mutation_stats[key1][target][key3] *= params["mutation_stats_del_scale"]

    for key1 in ["ins+1", "ins+1_3mers"]:
        if target == "total":
            for key2 in mutation_stats[key1]:
                for key3 in mutation_stats[key1][key2]:
                    mutation_stats[key1][key2][key3] *= params["mutation_stats_ins_scale"]

        else:
            for key3 in mutation_stats[key1][target]:
                mutation_stats[key1][target][key3] *= params["mutation_stats_ins_scale"]

    # apply empirical correction to adjust nonsense probabilities
    for key1 in ["pairs", "pairs_3mers"]:
        if target == "total":
            for key2 in mutation_stats[key1]:
                for key3 in mutation_stats[key1][key2]:
                    mutation_stats[key1][key2][key3] *= params["mutation_stats_non_scale"]

        else:
            for key3 in mutation_stats[key1][target]:
                mutation_stats[key1][target][key3] *= params["mutation_stats_non_scale"]

    return mutation_stats


def _analyze_predictions():
    if params["harmonize"] == True:
        #apu.predictions = apu.apply_variant_filter(apu.predictions, apu.params["variant_filter_col"].replace("ID:", ""), apu.params["variant_filter_col"], mode="variants") # <- removed on 250706
        apu.predictions = apu.apply_variant_filter(apu.predictions, params["variant_filter_col"].replace("ID:", ""), params["variant_filter_col"], mode="variants") # <- added on 250706

    # show data shapes (should be identical if perfect matching is indended)
    print("<", apu.variants.drop_duplicates(subset=params["variant_filter_col"]).shape[0], "mutated genes and",
               apu.predictions.drop_duplicates(subset=params["variant_filter_col"].replace("ID:", "")).shape[0], "reference genes are analyzed")

    if params["analysis_steps"] > 1:
        predictions = apu.predictions.applymap(copy.deepcopy)

    for step in range(params["analysis_steps"]):
        if step > 0:
            apu.masking_stats = pd.DataFrame()
            apu.predictions   = predictions.applymap(copy.deepcopy)
            apu.stats         = None

        # map variants to predictions
        apu.map_variants()

        if step == params["analysis_steps"]-1:
            # conduct selection tests for all variants combined
            apu.analyze_total_selection()

            # store test results
            apu.store_test()

            # remove entry for combined statistics
            apu.predictions = apu.predictions[apu.predictions["block id"] != "total"]

            if params["plot_projections"] == True:
                # initialize statistics
                apu.init_stats()

                # calculate statistics
                bar = IncrementalBar(set_bar("calculating prediction stats"), max=apu.predictions.shape[0])
                for i in range(apu.predictions.shape[0]):
                    apu.get_prediction_stats(i)
                    bar.next()

                bar.finish()

                # calculate projections
                if params["variants_only"] == False: apu.get_projection()

                # save statistics
                apu.save_stats()

        # <- if-condition added on 250612 
        if params["mode"] == "analyze_predictions" and params["apply_mutation_stats"] == True:
            apu.masking_stats["gene symbol"]         = [apu.variants[apu.variants["ID:transcript id"] == apu.masking_stats.iloc[i].loc["transcript id"]].iloc[0].loc["ID:gene symbol"]
                                                        for i in range(apu.masking_stats.shape[0])]
            frameshift_variants                      = apu.variants[apu.variants["ID:variant classification"].isin(["frame_shift_del", "frame_shift_ins"])]
            nonsense_variants                        = apu.variants[apu.variants["ID:variant classification"] == "nonsense"]
            deletion_variants                        = apu.variants[apu.variants["ID:variant classification"].isin(["frame_shift_del"])]
            insertion_variants                       = apu.variants[apu.variants["ID:variant classification"].isin(["frame_shift_ins"])]

            apu.masking_stats["observed frameshift"] = [frameshift_variants[frameshift_variants["ID:transcript id"] == apu.masking_stats.iloc[i].loc["transcript id"]].shape[0]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["observed del-1"]      = [deletion_variants[deletion_variants["ID:transcript id"] == apu.masking_stats.iloc[i].loc["transcript id"]].shape[0]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["observed ins+1"]      = [insertion_variants[insertion_variants["ID:transcript id"] == apu.masking_stats.iloc[i].loc["transcript id"]].shape[0]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["observed nonsense"]   = [nonsense_variants[nonsense_variants["ID:transcript id"] == apu.masking_stats.iloc[i].loc["transcript id"]].shape[0]
                                                        for i in range(apu.masking_stats.shape[0])]
                    
            apu.masking_stats["frameshift distance"] = [apu.masking_stats.iloc[i].loc["observed frameshift"]-apu.masking_stats.iloc[i].loc["proj. frameshift"]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["del-1 distance"]      = [apu.masking_stats.iloc[i].loc["observed del-1"]-apu.masking_stats.iloc[i].loc["proj. del-1"]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["ins+1 distance"]      = [apu.masking_stats.iloc[i].loc["observed ins+1"]-apu.masking_stats.iloc[i].loc["proj. ins+1"]
                                                        for i in range(apu.masking_stats.shape[0])]
            apu.masking_stats["nonsense distance"]   = [apu.masking_stats.iloc[i].loc["observed nonsense"]-apu.masking_stats.iloc[i].loc["proj. nonsense"]
                                                        for i in range(apu.masking_stats.shape[0])]
            
            masking_stats = pd.DataFrame({"block"               : [params["mutation_stats_gene_target"]],
                                          "transcript id"       : [params["mutation_stats_gene_target"]],
                                          "gene symbol"         : [params["mutation_stats_gene_target"]],
                                          "nonsense"            : [apu.masking_stats["nonsense"].mean()],
                                          "frameshift"          : [apu.masking_stats["frameshift"].mean()],
                                          "del-1"               : [apu.masking_stats["del-1"].mean()],
                                          "ins+1"               : [apu.masking_stats["ins+1"].mean()],
                                          "proj. nonsense"      : [apu.masking_stats["proj. nonsense"].sum()],
                                          "proj. frameshift"    : [apu.masking_stats["proj. frameshift"].sum()],
                                          "proj. del-1"         : [apu.masking_stats["proj. del-1"].sum()],
                                          "proj. ins+1"         : [apu.masking_stats["proj. ins+1"].sum()],
                                          "observed nonsense"   : [apu.masking_stats["observed nonsense"].sum()],
                                          "observed frameshift" : [apu.masking_stats["observed frameshift"].sum()],
                                          "observed del-1"      : [apu.masking_stats["observed del-1"].sum()],
                                          "observed ins+1"      : [apu.masking_stats["observed ins+1"].sum()],
                                          "nonsense distance"   : [apu.masking_stats["nonsense distance"].mean()],
                                          "frameshift distance" : [apu.masking_stats["frameshift distance"].mean()],
                                          "del-1 distance"      : [apu.masking_stats["del-1 distance"].mean()],
                                          "ins+1 distance"      : [apu.masking_stats["ins+1 distance"].mean()],
                                          })
            
            apu.masking_stats = pd.concat([apu.masking_stats, masking_stats]).sort_values(by="observed nonsense", ascending=False)
            params["mutation_stats_non_scale"] = masking_stats.iloc[0].loc["observed nonsense"]/masking_stats.iloc[0].loc["proj. nonsense"]
            params["mutation_stats_del_scale"] = masking_stats.iloc[0].loc["observed del-1"]/masking_stats.iloc[0].loc["proj. del-1"]
            params["mutation_stats_ins_scale"] = masking_stats.iloc[0].loc["observed ins+1"]/masking_stats.iloc[0].loc["proj. ins+1"]
            print("< scales", params["mutation_stats_non_scale"], "/", params["mutation_stats_del_scale"], "/", params["mutation_stats_ins_scale"])

            pd.set_option('display.max_columns', None)
            print(apu.masking_stats)

            if step == params["analysis_steps"]-1:
                apu.masking_stats.to_csv(apu.newdir+params["os_sep"]+params["file_tag"]+"_mutation_stats.txt", index=False)

            # adjust probabilites if multiple steps are selected
            else:
                apu.mutation_stats = _adjust_probabilities(apu.mutation_stats, params["mutation_stats_pair_target"])
                apu.print_params()

        if len(params["masks"]) > 0 and step == params["analysis_steps"]-1:
            print({key: len(apu.errors[key]) for key in apu.errors})
            with open(apu.newdir+params["os_sep"]+params["file_tag"]+"_errors.json", "w") as f:
                f.write(json.dumps(apu.errors, indent=4))


def _analyze_selection(block=None):   
    if block == None: fname = params["selection_fname"]
    if block != None: fname = params["selection_fname"]+"_"+block+"_selection_stats.txt"

    if os.path.isfile(params["data_dir"]+params["os_sep"]+fname):
        genes_under_selection = pd.read_csv(params["data_dir"]+params["os_sep"]+fname, delimiter=",", low_memory=False)
        apu.analyze_selection(genes_under_selection, block=block)

    else:
        print("<", params["data_dir"]+params["os_sep"]+fname, "not found.")
        exit()


if params["mode"] == "analyze_blocks":
    with open(params["data_dir"]+params["os_sep"]+params["status_fname"], "r") as f:
        status = json.load(f)

    apu.load(params["variant_fnames"], mode="variants")
    apu.analyze_blocks(status)


if params["mode"] == "analyze_adaptation" or params["mode"] == "analyze_predictions" or params["mode"] == "analyze_selection":
    # load data
    if params["mode"] != "analyze_selection":
        apu.seqs       = pd.read_csv(params["data_dir"]+params["os_sep"]+params["seq_fname"])
        apu.seqs       = apu.seqs[[True if len(apu.seqs.iloc[i].loc["cds"]) % 3 == 0 else False for i in range(apu.seqs.shape[0])]]
        apu.seqs.index = [apu.seqs.iloc[i].loc["transcript id"].split(".")[0] for i in range(apu.seqs.shape[0])]

        # <- if-condition added on 250603
        if params["apply_mutation_stats"] == True:
            if params["average_indel_probabilities"] == True:
                for key1 in apu.mutation_stats:
                    if ("del" in key1 or "ins" in key1) and "full" not in key1:
                        for key2 in apu.mutation_stats[key1]:
                            # exception required because for ins+1, total is the sum of base statistics
                            for key3 in apu.mutation_stats[key1][key2]:
                                if key1 == "ins+1" and key3 != "total":
                                    apu.mutation_stats[key1][key2][key3] = apu.mutation_stats[key1.split("_")[0]+"_full"][key2]/4
                                
                                else:
                                    apu.mutation_stats[key1][key2][key3] = apu.mutation_stats[key1.split("_")[0]+"_full"][key2]

            if params["mutation_stats_pair_target"] != None:
                apu.mutation_stats = _adjust_probabilities(apu.mutation_stats, params["mutation_stats_pair_target"])
                    
        if params["use_variant_filter"] == True: apu.load(params["variant_filter_fnames"], mode="variant_filter")
        if params["load_lindeboom"] == False:    apu.load(params["prediction_fnames"], mode="predictions")
        if params["load_lindeboom"] == True:     apu.load(params["lindeboom_fnames"], mode="predictions")
        apu.load(params["variant_fnames"], mode="variants")


    if params["mode"] != "analyze_adaptation" and params["block"] != None:
        file_tag = params["file_tag"]

        set_gene_target = False; set_pair_target = False
        if params["mutation_stats_gene_target"] == None: set_gene_target = True
        if params["mutation_stats_pair_target"] == None: set_pair_target = True

        # copy attribute state
        if params["mode"] != "analyze_selection":
            if params["harmonize"] == True:
                apu.predictions = apu.apply_variant_filter(apu.predictions, params["variant_filter_col"].replace("ID:", ""), params["variant_filter_col"], mode="variants")
            
            #predictions = apu.predictions # <- removed
            predictions    = apu.predictions.applymap(copy.deepcopy)
            mutation_stats = copy.deepcopy(apu.mutation_stats)
            
            if "patientwise" in params["profile"]:
                sample_ids          = np.unique(apu.variants[params["block"]])
                filtered_sample_ids = [sample_id for sample_id in sample_ids
                                       if apu.variants[apu.variants[params["block"]] == sample_id].shape[0] >= params["min_patients"]]
                apu.variants        = apu.variants[apu.variants[params["block"]].isin(filtered_sample_ids)]
                print("<", len(filtered_sample_ids), "patients selected.")
            
            variants    = apu.variants
            blocks      = apu.variants.drop_duplicates(subset=params["block"])[params["block"]].tolist()
            print(blocks)

        if params["mode"] == "analyze_selection":
            paths  = os.listdir(params["data_dir"])
            blocks = [path.split("_")[len(path.split("_"))-3] for path in paths if "selection_stats" in path]

        for block in blocks:
            params["file_tag"] = file_tag + "_" + block

            # reinitialize class attributes
            apu.base_dict      = None
            apu.means          = {"gene id": [], "mean": []}

            if params["mode"] != "analyze_selection":
                apu.predictions    = predictions.applymap(copy.deepcopy)
                apu.stats          = None
                apu.variants       = variants[variants[params["block"]] == block]
                apu.mutation_stats = copy.deepcopy(mutation_stats)

                if "projectwise" in params["profile"]:
                    apu.masking_stats = pd.DataFrame()
                    if set_gene_target == True: params["mutation_stats_gene_target"] = block
                    if set_pair_target == True: params["mutation_stats_pair_target"] = block
                
                if "patientwise" in params["profile"]:
                    apu.masking_stats = pd.DataFrame()

                    if params["id_filter"] == "MSK":
                        if set_gene_target == True: params["mutation_stats_gene_target"] = apu.variants.iloc[0].loc["ID:cancer type"]
                        if set_pair_target == True: params["mutation_stats_pair_target"] = apu.variants.iloc[0].loc["ID:cancer type"]

                    elif params["id_filter"] == "CPTAC3" or params["id_filter"] == "TCGA":
                        if set_gene_target == True: params["mutation_stats_gene_target"] = apu.variants.iloc[0].loc["ID:project"]
                        if set_pair_target == True: params["mutation_stats_pair_target"] = apu.variants.iloc[0].loc["ID:project"]

                    else:
                        print("< unknown datatype", params["datatype"])
                        exit()

            if params["mode"] == "analyze_adaptation" or params["mode"] == "analyze_predictions":
                _analyze_predictions()
                if params["mode"] == "analyze_adaptation": apu.analyze_adaptation()

            if params["mode"] == "analyze_selection":
                _analyze_selection(block)

    else:
        if params["mode"] == "analyze_adaptation" or params["mode"] == "analyze_predictions":
            _analyze_predictions()
            if params["mode"] == "analyze_adaptation": apu.analyze_adaptation()

        if params["mode"] == "analyze_selection":
             _analyze_selection()


if params["mode"] == "analyze_stop_codons":
    if params["use_variant_filter"] == True: apu.load(params["variant_filter_fnames"], mode="variant_filter")
    apu.load(params["variant_fnames"], mode="variants")
    apu.analyze_stop_codons()


if params["mode"] == "convert_lindeboom_predictions":
    # load lindeboom file
    with open(params["data_dir"]+params["os_sep"]+params["lindeboom_fnames"][0], "r") as f:
        lines = f.readlines()

    # create gene dictionary
    apu.convert_lindeboom_predictions(lines)


if params["mode"] == "show_stats":
    apu.seqs       = pd.read_csv(params["data_dir"]+params["os_sep"]+params["seq_fname"])
    apu.seqs       = apu.seqs[[True if len(apu.seqs.iloc[i].loc["cds"]) % 3 == 0 else False for i in range(apu.seqs.shape[0])]]
    apu.seqs.index = [apu.seqs.iloc[i].loc["transcript id"].split(".")[0] for i in range(apu.seqs.shape[0])]

    pu = Plot_utils()
    apu.load(params["prediction_fnames"], mode="predictions", loading_key={"transcript id": params["selected_transcripts"]})
    apu.load(params["variant_fnames"], mode="variants", loading_key={"ID:transcript id": [transcript.split(".")[0] for transcript in params["selected_transcripts"]]})
    apu.map_variants()

    gene_symbols = {"ENST00000269305": "TP53", "ENST00000257430": "APC"}
    for transcript_id in params["selected_transcripts"]:
        fig, ax = plt.subplots(1, 2)

        transcript_ids = apu.predictions[apu.predictions["transcript id"] == transcript_id.split(".")[0]]["transcript id"].tolist()
        for transcript_id in transcript_ids:
            print("<", transcript_id)
            data = {"model predictions"     : apu.predictions.loc[transcript_id].loc["predictions"],
                    "exon predictions"      : apu.predictions.loc[transcript_id].loc["predictions_by_exon"],
                    "real positions"        : apu.predictions.loc[transcript_id].loc["variant targets"]["FEATURE:ptc cds position"],
                    "real predictions"      : apu.predictions.loc[transcript_id].loc["variant targets"]["FEATURE:prediction"]}

            if len(params["masks"]) > 0: remove_placeholders = True
            else:                        remove_placeholders = False

            ax[0] = pu.plot_coordinates(ax[0], data, features={"label"               : gene_symbols[transcript_id],
                                                               "remove_placeholders" : remove_placeholders,
                                                               "stats"               : True,
                                                               "xlabel"              : "coordinate",
                                                               "ylabel"              : "NMD score",
                                                               "yrange"              : (0.5, 0.8)})

            ax[1] = pu.plot_gene_histogram(ax[1], data, features={"bins"                : 20,
                                                                  "label"               : gene_symbols[transcript_id],
                                                                  "remove_placeholders" : remove_placeholders,
                                                                  "prediction_mode"     : "histogram",
                                                                  "stats"               : True,
                                                                  "xlabel"              : "NMD score",
                                                                  "ylabel"              : "rel. counts"})

        plt.show()
        fig.savefig(apu.newdir+params["os_sep"]+gene_symbols[transcript_id]+".svg", dpi=600)

        with open(apu.newdir+params["os_sep"]+gene_symbols[transcript_id]+".json", "w") as f:
            f.write(json.dumps(data))