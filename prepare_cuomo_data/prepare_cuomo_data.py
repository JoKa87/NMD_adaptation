import os
import pandas as pd
from progress.bar import IncrementalBar
# marked (<-) added on 250427 to allow integration of hg38 in combination with liftover
from pyliftover import LiftOver # <-
import sys
import time

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\prepare_data")

from prepare_data_utils import *
        

# parameters
params = {
         "data_dir1"             : r"C:\Programming\Translational_genomics\NMD_analysis\data",
         "data_dir2"             : r"C:\Programming\Translational_genomics\NMD_analysis\cuomo_et_al\ase_aggregated_by_donor_open_access_lines",
         # <- added on 250616 from prepare_data
         "error_handling"        : { # defines handling of errors during sequence information extraction (True means excluding errors)
                                    "chromosome_incompatible"         : True, # in absolute mode, alternative chromosomes cannot be used
                                    "chromosome_not_found"            : True,
                                    # this error is taken out because some entries are outdated
                                    "inconsistent_cds_size"           : True, # only applied if value is provided
                                    "motif_start_not_found"           : True,
                                    "mutated_wt_base_not_found"       : False,
                                    # marked (<-) added to account for no cds (relevant for calculation of "FEATURE:cds EJC density" and "FEATURE:downstream cds EJC density")
                                    "no_cds"                          : True, # <-
                                    "no_exons"                        : True,
                                    "no_start"                        : True,
                                    "no_stop"                         : True,
                                    "ptc_cds_index_not_found"         : True,
                                    "ptc_equals_stop_codon"           : True,
                                    "ptc_exon_index_not_found"        : True,
                                    "ptc_not_in_cds"                  : True,
                                    "ptc_not_in_exon"                 : True,
                                    "ptc_in_utr"                      : True
                                    },
         "evaluation_mode"       : "sums", # "separation", "sums"
         "features"              : ["ID:variant id", "ID:gene id", "ID:transcript id", "ID:gene symbol", "ID:uc id", "ID:cell id", "ID:sample id", # removed "ID:alt. variant id", "ID:isoform count",
                                    "ID:appris annotation", "ID:strand",
                                    "LABEL:alt counts", "LABEL:total counts",

                                    "FEATURE:abs. ptc index", "FEATURE:total exon size", "FEATURE:total cds size", "FEATURE:5'utr-size", "FEATURE:3'utr-size", 
                                    "FEATURE:cds", "FEATURE:exons", "FEATURE:upstream exons", "FEATURE:downstream exons",
                                    # marked (<-) features added on 250425 to account for cds exons separately
                                    "FEATURE:cds exons", "FEATURE:upstream cds exons", "FEATURE:downstream cds exons", # <- added
                                    "FEATURE:EJC density", "FEATURE:downstream EJC density",
                                    # marked (<-) features added on 250425 to account for cds exons separately
                                    "FEATURE:cds EJC density", "FEATURE:downstream cds EJC density", # <- added
                                    "FEATURE:ptc upstream distance", "FEATURE:ptc downstream distance", "FEATURE:dist. from last EJC",
                                    "FEATURE:ptc cds position", "FEATURE:ptc exon position", "FEATURE:ptc-wt stop codon distance",
                                    # marked (<-) features added on 250425 to integrate kim features (+ cds specific features)
                                    "FEATURE:last exon", "FEATURE:last cds exon", "FEATURE:50 nt to last EJC", "FEATURE:50 nt to last cds EJC", "FEATURE:ptc exon size", # <- added

                                    "FEATURE:total start count cds", "FEATURE:upstream start ptc count cds", "FEATURE:downstream start ptc count cds",
                                    "FEATURE:first upstream start ptc distance cds", "FEATURE:first downstream start ptc distance cds",
                                    "FEATURE:total start count ptc cds", "FEATURE:upstream start ptc count ptc cds", "FEATURE:downstream start ptc count ptc cds",
                                    "FEATURE:first upstream start ptc distance ptc cds", "FEATURE:first downstream start ptc distance ptc cds",

                                    "FEATURE:total GC count exons", "FEATURE:upstream GC ptc count exons", "FEATURE:downstream GC ptc count exons",
                                    "FEATURE:total GC count ptc exon", "FEATURE:upstream GC ptc count ptc exon", "FEATURE:downstream GC ptc count ptc exon",

                                    "FEATURE:5'utr", "FEATURE:ptc", "FEATURE:5'ptc", "FEATURE:3'ptc", "FEATURE:5'ejc", "FEATURE:3'ejc",
                                    "FEATURE:ptc cds", "FEATURE:3'utr", "FEATURE:all cds", "FEATURE:all exons"],
         "file_tag"              : "ptc_variants",
         "filter"                : {
                                    "appris"            : ["-"], # None
                                     # "exclusive" "inclusive"
                                     # exclusive: appris annotation containing elements of "appris" list are removed (exact match)
                                     # inclusive: only appris annotation containing elements of "appris" list are kept ("appris" keyword contained in appris annotation)
                                    "appris_mode"       : "exclusive",
                                     # if block filtering is chosen, only the highest ranking appris priorities are kept
                                    "appris_priorities" : ["PRINCIPAL:1", "PRINCIPAL:2", "PRINCIPAL:3", "PRINCIPAL:4", "PRINCIPAL:5", "ALTERNATIVE:1", "ALTERNATIVE:2", "MINOR"],
                                    "appris_redundant"  : False, # if True, highest ranking appris terms are all exported, if False, one is selected randomly
                                    "block_targets"     : ["ID:cell id", "ID:sample id", "ID:variant id"], # "ID:variant id",
                                    "duplicate_targets" : ["ID:transcript id"], # entries are used for duplicate removal in blocks
                                    "hg38"              : True, # if True, entries are kept only if present in hg38_knownGene database
                                    "min_blocksize"     : 1,
                                    "min_reads"         : 8,
                                    "variant_target"    : "ID:variant id"
                                   },
         "hg_build"              : "hg19", # "hg19" "hg38"
         "max_procs"             : 1,
         "min_reads"             : 1, # refers to minimum total reads (>=)
         "mode"                  : "evaluation", # "calculation", "evaluation", "filtering"
         "motifs"                : {
                                    "start": ["ATG"],
                                    "GC":    ["G", "C"]
                                   },
         "motif_distance"        : {"start": True, "GC": False},
         "motif_inframe"         : {"start": True, "GC": False},
         "motif_relative"        : {"start": False, "GC": True},
         "os_sep"                : "\\",
         "overwrite"             : False,
         "pattern_cutoff"        : 5,
         "ptc_mode"              : "absolute",
         "run_dir"               : r"C:\Programming\Translational_genomics\NMD_analysis\cuomo_et_al\ptc_variants_hg19",
         # if selected_features = [], only specified features are stored (relevant for non-ptc variant extraction)
         "selected_features"     : ["ID:variant id", "ID:gene id", "ID:transcript id", "ID:gene symbol", "ID:uc id", "ID:cell id", "ID:sample id", # removed "ID:alt. variant id", "ID:isoform count",
                                    "ID:appris annotation", "ID:mutation type",
                                    "LABEL:alt counts", "LABEL:total counts",

                                    "FEATURE:abs. ptc index", "FEATURE:total cds size", "FEATURE:3'utr-size", "FEATURE:downstream exons", "FEATURE:downstream EJC density",
                                    "FEATURE:ptc downstream distance", "FEATURE:dist. from last EJC", "FEATURE:ptc cds position", "FEATURE:ptc exon position", "FEATURE:ptc-wt stop codon distance",
                                    "FEATURE:downstream GC ptc count ptc exon"],
         "stepsize"              : 50000,
         "target_mutations"      : "ptc" # "ptc" "non-ptc"
         }


def _run(pu, selected_paths, proc_index, params, lock, thread_id):
    ptc_count  = 0
    totalcount = 0
    for i in range(proc_index[0], proc_index[1]+1, 1):
        outfname = selected_paths[i].split(".")[0] + "_" + params["file_tag"]
        # check if outpath already exists
        if (params["overwrite"] == False and os.path.isfile(params["run_dir"]+params["os_sep"]+outfname) == False) or params["overwrite"] == True:
            print(selected_paths[i])
            ase     = pd.read_csv(params["data_dir2"]+params["os_sep"]+selected_paths[i], delimiter="\t", index_col=False)

            if params["target_mutations"] == "ptc": mutants = {col: [] for col in params["features"]}
            else:                                   mutants = {col: [] for col in params["selected_features"]}
            
            # marked (<-) added / removed on 250616
            # mutants = pu.assign_mutants(mutants, ase, lock, outfname, thread_id) # <- removed
            mutants, error_report = pu.assign_mutants(mutants, ase, lock, outfname, thread_id) # <- added

            # load respective path with "altcount" tag and read count data
            ase     = pd.read_csv(params["data_dir2"]+params["os_sep"]+selected_paths[i].replace("totalcount", "altcount"), delimiter="\t", index_col=False)
            mutants = pd.DataFrame(mutants)
            mutants = pu.read_counts(mutants, ase)

            ptc_count  += mutants.shape[0]
            totalcount += ase.shape[0]

            lock.acquire()
            try:
                print("thread", thread_id, "mutants", mutants.shape[0], "of", totalcount, "rate", ptc_count/totalcount)

            finally:
                lock.release()

            # print to file
            mutants.to_csv(path_or_buf=params["run_dir"]+params["os_sep"]+outfname, sep=",", index=False)

            with open(params["run_dir"]+params["os_sep"]+outfname+"_error_report.json", "w") as f:
                f.write(json.dumps(error_report, indent=4))
    

def main():
    if params["mode"] == "calculation":
        # load ASE data and conduct analysis
        paths = os.listdir(params["data_dir2"])

        # initialize class for preparation
        pu = Prepare_cuomo_data_utils(params)
        
        selected_paths = []
        for path in paths:
            if "totalcount" in path:
                outfname = path.split(".")[0] + "_" + params["file_tag"]
                if os.path.isfile(params["run_dir"]+params["os_sep"]+outfname) == False:
                    selected_paths.append(path)

        proc_index = split_index(len(selected_paths), params["max_procs"])
        lock       = Lock()

        procs = []
        for i in range(len(proc_index)):
            proc = Process(target=_run, args=(pu, selected_paths, proc_index[i], params, lock, i))
            procs.append(proc)
            proc.start()

        for i, proc in enumerate(procs):
            proc.join()


    if params["mode"] == "evaluation":
        eu          = Evaluation_utils(params)
        paths       = os.listdir(params["run_dir"])
        all_mutants = []
        for path in paths:
            if params["file_tag"] in path and os.path.isfile(params["run_dir"]+params["os_sep"]+path) == True:
                mutants                 = pd.read_csv(params["run_dir"]+params["os_sep"]+path, delimiter=",", index_col=False)
                mutants["ID:source id"] = [path.split("_")[0] for _ in range(mutants.shape[0])]
                all_mutants.append(mutants)

        all_mutants = pd.concat(all_mutants, ignore_index=True)
        all_mutants = all_mutants.sort_values(by=["ID:gene id", "ID:variant id"])
        
        if params["evaluation_mode"] == "separation": all_mutants = eu.separate_counts(all_mutants)
        if params["evaluation_mode"] == "sums":       all_mutants = eu.sum_counts(all_mutants)
        all_mutants.to_csv(path_or_buf=params["run_dir"]+params["os_sep"]+params["file_tag"]+"_"+params["evaluation_mode"], sep=",", index=False)


    if params["mode"] == "filtering":
        eu = Evaluation_utils(params)

        if params["evaluation_mode"] == "separation":
            all_mutants = pd.read_csv(params["data_dir1"]+params["os_sep"]+"ptc_variants_separation", delimiter=",", index_col=False)

        if params["evaluation_mode"] == "sums":
            all_mutants = pd.read_csv(params["data_dir1"]+params["os_sep"]+"ptc_variants_sums", delimiter=",", index_col=False)

        all_mutants = eu.apply_filter(all_mutants)
        all_mutants.to_csv(path_or_buf=params["data_dir1"]+params["os_sep"]+params["file_tag"]+"_"+params["evaluation_mode"]+"_filtered", sep=",", index=False)


if __name__ == '__main__':
    main()