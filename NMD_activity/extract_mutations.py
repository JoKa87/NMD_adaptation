import platform
from multiprocessing import Process, Lock
import time

from extract_mutations_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


# parameters
params = {
            "aggregate_dataset"   : False,
            "calculate_complements": True, # <- added on 250623
            "cnv_threshold"       : -1, # -1 to exclude misses (CNV data contained in RNA_c)
            "create_calculations" : False,
            "create_dataset"      : False,
            "create_mutation_stats": True, # <- added on 250603
            "data_dir"            : parent_dir+r"\data", # TCGA must be located here
            # added on 251027 to allow alternative project keys (required for CPTAC-3), None to switch off
            "external_project_path": None, # r"C:\Programming\Translational_genomics\NMD_analysis\data\cptac3_clinical.tsv",
            "extraction_stages"   : ["ptc_variants"], # "ptc_variants" "rnaseq" "calculations"
            # full_extraction allows extraction of all mutations with columns specified in extraction_targets (empty list means all values are extracted)
            "extraction_targets"  : {"Hugo_Symbol": [], "Gene": [], "Transcript_ID": [], "case_id": [], # identifiers
                                     "Chromosome": [], "Start_Position": [], "End_Position": [], "Strand": [], # localization
                                     "SIFT": [], "PolyPhen": [], # mutational outcome prediction
                                     "HGVSp": [], # mutational info
                                     "Variant_Classification": ["In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins",
                                                                "Missense_Mutation", "Nonsense_Mutation", "Silent",
                                                                "Splice_Site", "In_Frame_Del", "Nonstop_Mutation"]},
            "fname"               : "mutation_stats_test.txt",
            "get_mutation_info"   : False, # if True all mutations for the given PTC variant are extracted
            "ids"                 : [
                                    "variant_id", # unique variant identifier according to Teran et al.
                                    "project",
                                    "cluster_key"
                                    ],
            "info"                : {
                                     "rna": [
                                             "RNASEQ_noptc_gene_id", "RNASEQ_noptc_unstranded", "RNASEQ_noptc_stranded_first", "RNASEQ_noptc_stranded_second", "RNASEQ_noptc_tpm_unstranded",
                                             "RNASEQ_noptc_fpkm_unstranded", "RNASEQ_noptc_fpkm_uq_unstranded",
                                             "RNASEQ_noptc_cnv_total", "RNASEQ_noptc_cnv_minor", "RNASEQ_noptc_cnv_avg",
                                             "RNASEQ_ptc_gene_id", "RNASEQ_ptc_unstranded", "RNASEQ_ptc_stranded_first", "RNASEQ_ptc_stranded_second", "RNASEQ_ptc_tpm_unstranded",
                                             "RNASEQ_ptc_fpkm_unstranded", "RNASEQ_ptc_fpkm_uq_unstranded",
                                             "RNASEQ_ptc_cnv_total", "RNASEQ_ptc_cnv_minor", "RNASEQ_ptc_cnv_avg" # from rnaseq data
                                            ]
                                    },
            "masks"               : ["nonsense", "frameshift", "missense", "silent"], # "frameshift" "nonsense" [] <- added on 250618
            "max_procs"           : 1,
            # full_extraction allows extraction of all mutations with columns specified in extraction_targets
            # with use_targets=True, ptc_extraction still allows extraction of specified mutations
            # full_extraction does not involve further processing steps of PTC mutation info
            "mode"                : "ptc_extraction", # "full_extraction" "ptc_extraction"
            "mutation_stats_ptc_weights": False, # <- added on 251025
            "mutation_targets"    : {},
            "os_sep"              : "//",
            "projects"            : "all", # ["TCGA-BRCA"], # "all" or list of projects ["TCGA-CHOL", "TCGA-LAML"]
            "randomize_cluster"   : True, # <- added on 251024
            "rna_calculation"     : ["avg", "avg", "avg"], # "median" in Lindeboom et al. 2016
            "rna_criterion"       : "tpm_unstranded",
            "rna_noptc_criterion" : "noptc_tpm_unstranded",
            "rna_threshold"       : 0, # None if no threshold should be applied
            "rna_values"          : [
                                     "RNASEQ_ptc_unstranded", "RNASEQ_ptc_stranded_first", "RNASEQ_ptc_stranded_second", "RNASEQ_ptc_tpm_unstranded",
                                     "RNASEQ_ptc_fpkm_unstranded", "RNASEQ_ptc_fpkm_uq_unstranded",
                                     "RNASEQ_noptc_unstranded", "RNASEQ_noptc_stranded_first", "RNASEQ_noptc_stranded_second", "RNASEQ_noptc_tpm_unstranded",
                                     "RNASEQ_noptc_fpkm_unstranded", "RNASEQ_noptc_fpkm_uq_unstranded",
                                     "RNASEQ_noptc_cnv_total", "RNASEQ_noptc_cnv_minor", "RNASEQ_noptc_cnv_avg",
                                     "RNASEQ_ptc_cnv_total", "RNASEQ_ptc_cnv_minor", "RNASEQ_ptc_cnv_avg" # RNASEQ_cnv_total is cnv read-out, -1 means cnv status unknown
                                    ],
            "run_dir"             : parent_dir+r"\data",
            "seq_fname"           : "hg38_seqs_appended.txt", # <- added on 250618
            "status_fname"        : "filtered_status_ptc.json", # "filtered_status_cpatc3.json", "filtered_status_with_sclc.json", "filtered_extracted_status.json" "filtered_status_ptc.json" 
            "target_fname"        : "test_genes.txt",
            "target_identifier"   : {"rna": "gene_id", "wxs": "Gene"},
            "targets"             : [],
            "transform_type"      : "RNA_c_ptc", # if not None, transformed rna-seq data are used
            "use_targets"         : False,
            "wxs_cols"            : [],
            "wxs_identifier"      : "Transcript_ID" # wxs column that is used to create variant id and conduct ptc checks
        }

if platform.system() == "Windows": params["os_sep"] = "\\"

'''
the following WXS file columns that are currently not used in post-processing trigger pandas warning currently ignored for better loading performance (250613) 
74 Consequence
78 Amino_acids
112 gnomAD_AFR_AF
114 gnomAD_ASJ_AF
120 MAX_AF
132 gnomAD_non_cancer_SAS_AF
137 PUBMED
155 GDC_FILTER
'''

def get_clusters_(cluster, params):
    clusters = []
    for project_key in cluster:
        #if "TCGA" in project_key:
        for cluster_key in cluster[project_key]:
            if (project_key in params["projects"] or params["projects"] == "all"):
                if project_key == "TCGA-PAAD" and cluster_key == "cnmf_not_assigned":
                    pass
                    
                else:
                    clusters.append({"project_key": project_key, "cluster_key": cluster_key})
        
    print("<", len(clusters), "clusters selected.")
    return clusters


def extract_(cluster):
    clusters = get_clusters_(cluster, params)

    rand_index = np.arange(len(clusters))
    np.random.shuffle(rand_index)
    clusters   = [{"project_key": clusters[rand_i]["project_key"], "cluster_key": clusters[rand_i]["cluster_key"]} for rand_i in rand_index]
    proc_index = split_index(len(clusters), params["max_procs"])
    lock       = Lock()
    procs      = []

    for i in range(len(proc_index)):
        proc = Process(target=extract, args=(cluster, clusters, proc_index[i], params, lock, i))
        procs.append(proc)
        proc.start()

    for i, proc in enumerate(procs):
        proc.join()

    
    return clusters


def main():
    start_time = time.time()
    if params["create_dataset"] == True:
        for extraction_stage in params["extraction_stages"]:
            if os.path.isdir(params["data_dir"]+params["os_sep"]+extraction_stage) == True:
                print("< warning.", extraction_stage, "folder already exists.")

        # load status data
        status  = load_status(params, fname=params["status_fname"])

        # initialize cluster dictionary
        cluster = init_cluster(status, params, cluster_key="firehouse")
        cluster = filter_cluster(cluster, filter_mode="rna_and_wxs_present")
        # <- added on 251024
        if params["randomize_cluster"] == True: cluster = randomize_cluster(cluster, params["max_procs"])

        # load target filenames (if selected)
        if params["use_targets"] == True:
            params["targets"] = load_df(params["run_dir"]+params["os_sep"]+params["target_fname"], delimiter=",")
            params["targets"] = params["targets"].drop_duplicates(subset=[params["target_identifier"]["rna"]])[params["target_identifier"]["rna"]].tolist()
            print(params["targets"])

        # extract ptc variants from wxs files, integrating cnv info
        # extract rna seq info for ptc variants
        clusters = extract_(cluster)
        print("< dataset calculation time:", round(time.time()-start_time, 4))


    if params["aggregate_dataset"] == True:
         # load status data
        status  = load_status(params, fname=params["status_fname"])

        # initialize cluster dictionary
        cluster = init_cluster(status, params, cluster_key="Firehouse")
        # <- added on 251024 (max_procs must match the no. of randomized clusters during calculation)
        if params["randomize_cluster"] == True: cluster = randomize_cluster(cluster, params["max_procs"])

        # create cluster list
        clusters = get_clusters_(cluster, params)
                    
        # loading from memory for aggregation
        cluster  = aggregate(cluster, clusters, params, folder=params["extraction_stages"][-1])


    if params["create_calculations"] == True:
        #  load dataframe
        df = load_df(params["run_dir"]+params["os_sep"]+params["fname"])

        # conduct calculations
        df = calculate_expression(df, params)

        # print calculated dataframe
        df.to_csv(params["run_dir"]+params["os_sep"]+params["fname"].replace("TCGA_", "TCGA_calculated_"), index=False, sep="\t") #TCGA_calculated_ptc_variants2.txt


    # marked (<-) if-statement added on 250603
    if params["create_mutation_stats"] == True:
        # load status data
        status = load_status(params, fname=params["status_fname"])

        # load cds sequences
        seqs                  = pd.read_csv(params["run_dir"]+params["os_sep"]+params["seq_fname"], sep=",")
        seqs["transcript id"] = [transcript_id.split(".")[0] for transcript_id in seqs["transcript id"]]
        seqs.index            = seqs["transcript id"]

        # initialize cluster dictionary
        cluster = init_cluster(status, params, cluster_key="Firehouse")

        # create cluster list
        clusters = get_clusters_(cluster, params)

        # create cluster list
        mutation_stats = create_mutation_stats(cluster, clusters, seqs, params)
        
        #with open(params["data_dir"]+params["os_sep"]+params["fname"].replace(".txt", ".json")) as f:
        #     mutation_stats = json.load(f)

        evaluate_mutation_stats(mutation_stats, seqs, params)


    # store parameters of current session
    params_fname = params["fname"].split(".")[0] + "_params.json"
    with open(params["run_dir"]+params["os_sep"]+params_fname, "w") as f:
        f.write(json.dumps(params, indent=4))


if __name__ == '__main__':
    main()