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
             "convert_label" : False,
             "data"          : [],
             "datatype"      : ["pandas"],
             "extensions"    : [None],
             "features"      : {
                                "boxtype"      : "compact",
                                "colors"       : ["royalblue", "crimson"],
                                "item"         : "",
                                "meanline"     : False,
                                "reverse"      : False,
                                "scalebar"     : True,
                                "showfliers"   : False,
                                "showmeans"    : True,
                                "tag"          : "target",
                                "xcol"         : [[None]],
                                "ycol"         : ["class 1", "class 2"],
                                "xmute"        : False,
                                "ylabel"       : "log2FC (SMG1i/DMSO)",
                                "ylabels"      : ["escape / \n target"]
                               },
             "filter"        : {
                                "allelic_fraction":     (0.4, 1),
                                "expression_threshold": 100,
                                "pred_threshold":       0.57
                               },
             "off"           : False,
             "paths"         : [parent_dir+r"\data\nmd_inhibition_variants.txt"],
             "separators"    : [","],
             "type"          : "box_plot"
            }
        ]
    

    dims       = [4, 4]
    resolution = 600
    run_dir    = parent_dir+r"\data"

    # loading section
    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)


    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4))
    plt.subplots_adjust(hspace=0.25)

    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 0:6]))
    subplots.append(fig.add_subplot(gs[1, 6:12]))

    for i in range(len(data)):
        if data[i]["type"] == "box_plot":
            # data filtering
            data[i]["data"][0] = data[i]["data"][0][data[i]["data"][0]["ID:allelic fraction"] >= data[i]["filter"]["allelic_fraction"][0]]
            data[i]["data"][0] = data[i]["data"][0][data[i]["data"][0]["ID:allelic fraction"] <= data[i]["filter"]["allelic_fraction"][1]]
            data[i]["data"][0] = data[i]["data"][0][data[i]["data"][0]["ID:base mean"] >= data[i]["filter"]["expression_threshold"]]

            if data[i]["convert_label"] == True:
                data[i]["data"][0]["LABEL:NMD score"] = convert_labels(data[i]["data"][0]["LABEL:NMD score"].tolist(), input_format="log2FC", output_format="ASE")

            class1 = data[i]["data"][0][data[i]["data"][0]["FEATURE:prediction"] <= data[i]["filter"]["pred_threshold"]]["LABEL:NMD score"].tolist()
            class2 = data[i]["data"][0][data[i]["data"][0]["FEATURE:prediction"] > data[i]["filter"]["pred_threshold"]]["LABEL:NMD score"].tolist()
            mw     = scipy.stats.mannwhitneyu(class1, class2)
            welch  = scipy.stats.ttest_ind(class1, class2, equal_var=False)
            print("< mw", round(mw.statistic, 4), round(mw.pvalue, 4), "welch", round(welch.statistic, 4), round(welch.pvalue, 4), "size", len(class1), "/", len(class2))
            
            temp_data                   = {"class 1": [class1], "class 2": [class2], "description": data[i]["features"]["ylabels"]}
            data[i]["features"]["xcol"] = "description"
            data[i]["features"]["ycol"] = ["class 1", "class 2"]
            subplots[i]                 = pu.box_plot(subplots[i], pd.DataFrame(temp_data), data[i]["features"])
            subplots[i].set_yticks([-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])

        if data[i]["type"] == "step_histogram":
            subplots[i] = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])
        
        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()
    fig.savefig(run_dir+"\\figures\\Fig2.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()