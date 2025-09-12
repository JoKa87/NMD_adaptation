import os

from plot_load import *
from plot_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    data = [
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                               "bins"       : 20,
                               "density"    : True,
                               "labels"     : ["GTex8 data (training & validation)", "single-cell data (testing)"],
                               "x"          : "LABEL:NMD score",
                               "xcol"       : [["LABEL:NMD score"], ["LABEL:NMD score"]],
                               "xlabel"     : "NMD susceptibility",
                               "ylabel"     : "rel. counts"
                              },
            "paths"         : [parent_dir+r"\data\model_training_variants.txt",
                               parent_dir+r"\data\model_test_variants.txt"],
            "separators"    : [",", ","],
            "type"          : "step_histogram"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                               "stats"      : "pearson2",
                               "x"          : "LABEL:NMD score",
                               "xcol"       : [[None]],
                               "y"          : "LABEL:avg. NMD score",
                               "xlabel"     : "NMD susceptibility of indidivual variant",
                               "ylabel"     : "averaged NMD susceptibility \n of identical variants"
                              },
            "paths"         : [parent_dir+r"\data\model_training_variants_full_block_correlation.txt"],
            "separators"    : [","],
            "type"          : "correlation"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "x"           : "average",
                              "xcol"        : [[None]],
                              "ylabel"      : "feature importance",# "mean decrease in impurity",
                              "z"           : "sem",
                              "zcutoff"     : (0, 0.01),
                              "zlabel"      : "sem"
                              },
            "paths"         : [parent_dir+r"\data\nmdelphi\forest_importances.txt"], # insert model folder name here
            "separators"    : [","],
            "type"          : "random forest coefficients"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                               "stats"      : "pearson2",
                               "x"          : "LABEL:NMD score",
                               "xcol"       : [[None]],
                               "y"          : "FEATURE:prediction",
                               "xlabel"     : "NMD susceptibility of indidivual variant",
                               "ylabel"     : "predicted NMD susceptibility"
                              },
            "paths"         : [parent_dir+r"\data\nmdelphi\single_unbalanced\all_preds.txt"],
            "separators"    : [","],
            "type"          : "correlation"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : ["balanced_roc"],
            "features"      : {
                               "x"           : "fpr",
                               "xcol"        : [[None]],
                               "y"           : "tpr",
                               "xlabel"      : "false positive rate",
                               "ylabel"      : "true positive rate"
                              },
            "paths"         : [parent_dir+r"\data\nmdelphi\balanced"],
            "separators"    : [","],
            "type"          : "auroc"
            }
        ]
    

    dims           = [3, 4]
    resolution     = 600
    run_dir        = parent_dir+r"\data"

    # loading section
    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)


    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4))
    plt.subplots_adjust(wspace=0.15, hspace=0.25)

    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 0:6]))
    subplots.append(fig.add_subplot(gs[0, 6:12]))
    subplots.append(fig.add_subplot(gs[1:3, 0:6]))
    subplots.append(fig.add_subplot(gs[1, 6:12]))
    subplots.append(fig.add_subplot(gs[2, 6:12]))

    for i in range(len(data)):
        if data[i]["type"] == "auroc":
            subplots[i] = pu.plot_auroc(subplots[i], data[i]["data"], data[i]["features"])

        if data[i]["type"] == "correlation":
            subplots[i] = pu.plot_correlation(subplots[i], data[i]["data"][0], data[i]["features"])

        if data[i]["type"] == "random forest coefficients":
            data[i]["data"][0].index = data[i]["data"][0][data[i]["data"][0].columns[0]]
            subplots[i]              = pu.plot_coefficients(subplots[i], data[i]["data"][0], data[i]["features"])

        if data[i]["type"] == "step_histogram":
            subplots[i] = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])
        
        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()
    fig.savefig(run_dir+"\\figures\\Fig2.svg", dpi=resolution, transparent=True)
    return


if __name__ == '__main__':
    main()