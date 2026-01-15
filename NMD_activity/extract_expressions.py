
import platform
from multiprocessing import Process, Lock
import time

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

from analyze_targets_utils import *
from extract_expressions_utils import *


# parameters
params = {
            "cluster_key"                : None, # "Firehouse", None # if None, calculation is conducted project-wise only and not cluster-resolved
            "data_dir"                   : parent_dir+r"data", # here, TCGA data should be located
            "file_tag"                   : "_TCGA_expressions",
            "max_procs"                  : 9,
            "mode"                       : "calculate_full_stats", # "calculate_full_stats" "create_dataset"
            "os_sep"                     : "//",
            "prediction_fnames"          : ["hg38_NMD_scores_appended.txt"],
            "processed_fname"            : "TCGA_expressions.txt",
            "projects"                   : "all", # "all"
            "rna_identifier"             : "gene_id",
            "rna_targets"                : ["fpkm_unstranded", "tpm_unstranded"], # "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", , "fpkm_uq_unstranded"
            "rna_threshold"              : 0,
            "run_dir"                    : parent_dir+r"data",
            "status_fname"               : "filtered_status_ptc.json", # "filtered_status_with_sclc.json", "filtered_extracted_status.json"
            "template_fname"             : "TCGA_rna_example",
            "transform_type"             : "RNA_ccorr_ptc",
            "verbose"                    : True
        }

if platform.system() == "Windows": params["os_sep"] = "\\"


def extract_(cluster, template):
    clusters = []
    for project_key in cluster:
        if "TCGA" in project_key:
            for cluster_key in cluster[project_key]:
                if (project_key in params["projects"] or params["projects"] == "all"):
                    clusters.append({"project_key": project_key, "cluster_key": cluster_key})

    print("<", len(clusters), "selected.")

    rand_index = np.arange(len(clusters))
    np.random.shuffle(rand_index)
    clusters   = [{"project_key": clusters[rand_i]["project_key"], "cluster_key": clusters[rand_i]["cluster_key"]} for rand_i in rand_index]
    proc_index = split_index(len(clusters), params["max_procs"])
    lock  = Lock()
    procs = []

    for i in range(len(proc_index)):
        proc = Process(target=extract_expressions, args=(cluster, clusters, template, proc_index[i], params, lock, i))
        procs.append(proc)
        proc.start()

    for i, proc in enumerate(procs):
        proc.join()

    return clusters


def main():
    if params["mode"] == "create_dataset":
        # save parameters
        params_json = json.dumps(params, indent=4)
        with open(params["run_dir"]+params["os_sep"]+params["processed_fname"].split(".")[0]+"_params.json", "w") as file:
            file.write(params_json)

        start_time = time.time()
        print("< calculation of", params["file_tag"].replace("_", " "))

        # load status data
        status  = load_status(params, fname=params["status_fname"])

        # initialize cluster dictionary
        cluster = init_cluster(status, params, cluster_key=params["cluster_key"])

        # filter cluster dictionary
        cluster = filter_cluster(cluster, filter_mode="rna_present")

        # load example file
        template = load_rna(params["run_dir"]+params["os_sep"]+params["template_fname"], params)

        # extract expression data
        clusters = extract_(cluster, template)
        print("< dataset calculation time:", round(time.time()-start_time, 4))

        # loading from storage for aggregation
        cluster = aggregate(cluster, clusters, params, axis=1, delimiter=",", folder="rnaseq"+params["file_tag"],
                            path=params["run_dir"]+params["os_sep"]+params["processed_fname"], redundant_cols=["gene_id", "gene_name", "gene_type"])


    if params["mode"] == "calculate_full_stats":
        stats       = {**{"mean_"+rna_target: [] for rna_target in params["rna_targets"]}, **{"median_"+rna_target: [] for rna_target in params["rna_targets"]}}
        expressions = pd.read_csv(params["run_dir"]+params["os_sep"]+params["processed_fname"], delimiter=",", low_memory=False)

        bar = IncrementalBar(set_bar("calculating stats"), max=expressions.shape[0]*len(params["rna_targets"])/100)
        for rna_target in params["rna_targets"]:
            for i in range(expressions.shape[0]):
                means   = [expressions.iloc[i].loc[col] for col in expressions.columns if "mean" in col and rna_target in col and pd.isna(expressions.iloc[i].loc[col]) == False]
                medians = [expressions.iloc[i].loc[col] for col in expressions.columns if "median" in col and rna_target in col and pd.isna(expressions.iloc[i].loc[col]) == False]

                if len(means) > 0:   stats["mean_"+rna_target].append(np.mean(means))
                else:                stats["mean_"+rna_target].append(None)
                if len(medians) > 0: stats["median_"+rna_target].append(np.mean(medians))
                else:                stats["median_"+rna_target].append(None)

                #print("i", i, len(means), len(medians), [col for col in data[0].columns if "median" in col])
                if i % 100 == 0: bar.next()

        for key in stats:
            expressions[key] = stats[key]

        bar.finish()
        expressions.to_csv(path_or_buf=params["run_dir"]+params["os_sep"]+params["processed_fname"].split(".")[0]+"_appended.txt", sep=",", index=False)


if __name__ == '__main__':
    main()