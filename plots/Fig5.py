import os
import statsmodels.api as sm

from plot_load import *
from plot_utils import *
from plot_tools import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    boxplot_colors = []
    [boxplot_colors.extend(["royalblue", "crimson"]) for _ in range(4)]

    # params #
    data = [
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "colors"      : ["royalblue", "crimson"],
                              "label_col"   : "class",
                              "size_col"    : "size",
                              "xcol"        : [[None]],
                              "xcols"       : ["NMD not mut.", "NMD mut."],
                              "xlabel"      : "NMD not mut. / \n NMD mut."
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\tcga_class_sizes.txt"],
            "separators"    : [","],
            "type"          : "pie"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "colors"      : ["royalblue", "crimson"],
                              "label_col"   : "class",
                              "size_col"    : "size",
                              "xcol"        : [[None]],
                              "xcols"       : ["low PTC", "high PTC"],
                              "xlabel"      : "low / \n high PTC"
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\tcga_class_sizes.txt"],
            "separators"    : [","],
            "type"          : "pie"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "colors"      : ["royalblue", "crimson"],
                              "label_col"   : "class",
                              "size_col"    : "size",
                              "xcol"        : [[None]],
                              "xcols"       : ["low frameshift", "high frameshift"],
                              "xlabel"      : "low / \n high frameshift"
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\tcga_class_sizes.txt"],
            "separators"    : [","],
            "type"          : "pie"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "colors":     {"data": "white", "class1": "dimgray", "class2": "lightgray"},
                              "edgecolors": {"data": "dimgray", "class1": "black", "class2": "black"},
                              "switch_axes": True,
                              "xcol":       [[None]],
                              "xlabel":     "PTC mutations",
                              "ylabel":     "frameshifts"
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_plot.txt"],
            "separators"    : [","],
            "type"          : "class_selection"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                               "features"           : ["FEATURE:ptc_mutations", "fpkm_unstranded", "FEATURE:prediction", "ID:cnv total"],
                               "normalization_mode" : "min-max",
                               "xcol"               : [[None]],
                               "ylabel"             : "distance",
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\cancer_scores_TCGA_NMD_targets_analysis_FPKM_exp_ccorr"],
            "separators"    : [","],
            "type"          : "dendrogram"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas", "pandas"],
            "extensions"    : [None],
            "features"      : {
                              "bar_off"     : True,
                              "include_empty": False,
                              "reverse"     : True,
                              "scale"       : {"binary": True, "ycol": [0.05, 0.01, 0.001], "zcol": [1], "colors": ["crimson", "royalblue"], "sizes": [0, 0.5, 1, 1.5]},
                              "tag"         : "target",
                              "xcol"        : [[None], [None], [None], [None]],
                              "xrange"      : (-0.5, 32.5),
                              "xticks"      : False,
                              "ylabels"     : ["NMD not mut. / \n nmd mut.", "low / \n high PTC", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],
            "separators"    : [",", ",", ",", ","],
            "type"          : "project_matrix"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas" for _ in range(4)],
            "extensions"    : [None for _ in range(4)],
            "features"      : {
                              "boxtype"     : "compact",
                              "colors"      : boxplot_colors,
                              "item"        : "NMD susceptibility",
                              "reverse"     : False,
                              "showfliers"  : False,
                              "tag"         : "target",
                              "xcol"        : [[None] for _ in range(4)],
                              "xmute"       : True,
                              "ylabels"     : ["NMD not mut. / \n NMD mut.", "low / \n high PTC", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],
            "separators"    : ["," for _ in range(4)],
            "type"          : "class_boxplot"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas" for _ in range(4)],
            "extensions"    : [None for _ in range(4)],
            "features"      : {
                              "boxtype"     : "compact",
                              "colors"      : boxplot_colors,
                              "item"        : "copy number",
                              "reverse"     : False,
                              "scalebar"    : True,
                              "showfliers"  : False,
                              "tag"         : "target",
                              "xcol"        : [[None] for _ in range(4)],
                              "xmute"       : True,
                              "ylabels"     : ["NMD not mut. / \n NMD mut.", "low / \n high PTC", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],
            "separators"    : ["," for _ in range(4)],
            "type"          : "class_boxplot"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas" for _ in range(3)],
            "extensions"    : [None for _ in range(3)],
            "features"      : {
                              "boxtype"     : "compact",
                              "colors"      : boxplot_colors,
                              "item"        : "NMD activity",
                              "reverse"     : False,
                              "showfliers"  : False,
                              "tag"         : "target",
                              "xcol"        : [[None] for _ in range(3)],
                              "xmute"       : False,
                              "ylabels"     : ["NMD not mut. / \n NMD mut.", "low / \n high PTC", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],
            "separators"    : ["," for _ in range(3)],
            "type"          : "class_boxplot"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas" for _ in range(12)],
            "extensions"    : [None for _ in range(12)],
            "features"      : {
                              "bar_off"     : True,
                              "layer"       : {0: {"marker": "o", "markersize": 12}},# 1: {"marker": "o", "markersize": 5}},
                              "reverse"     : False,
                              "scale"       : {"binary": True, "ycol": [0.05, 0.01, 0.001], "zcol": [0], "colors": ["royalblue", "crimson"], "sizes": [0, 0.5, 1, 1.5]},
                              "tag"         : "target",
                              "xcol"        : [[None] for _ in range(12)],
                              "ylabels"     : ["NMD not mut.", "NMD mut.", "low PTC", "high PTC", "low \n frameshift", "high \n frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class2_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class2_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class2_smoothing\stats_summary.txt"],
            "separators"    : ["," for _ in range(12)],
            "type"          : "correlation_matrix"
            }
        ]
    
    item_dict = {
                 "cnv total"                            : "copy number",
                 "frameshift_mutations"                 : "frameshift mutations",
                 "ptc_mutations2"                       : "PTC mutations",
                 "fpkm_unstranded"                      : "NMD activity",
                 "prediction"                           : "NMD susceptibility",
                 "ptc_mutations2_cnv total"             : "PTC mutations-copy number",
                 "ptc_mutations2_fpkm_unstranded"       : "PTC mutations-NMD activity",
                 "ptc_mutations2_prediction"            : "PTC mutations-NMD susceptibility",
                 }


    dims       = (6, 4)
    g          = (20, 10)
    margin     = (0.8, 0.2)
    resolution = 600
    run_dir    = parent_dir+r"\data"

    pu = Plot_utils()


    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(constrained_layout=False)
    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 0:1]))
    subplots.append(fig.add_subplot(gs[0, 1:2]))
    subplots.append(fig.add_subplot(gs[0, 2:3]))
    subplots.append(fig.add_subplot(gs[1, 0:3]))
    subplots.append(fig.add_subplot(gs[0, 3:6]))
    subplots.append(fig.add_subplot(gs[1, 3:6]))
    subplots.append(fig.add_subplot(gs[2, 0:2]))
    subplots.append(fig.add_subplot(gs[2, 2:4]))
    subplots.append(fig.add_subplot(gs[2, 4:6]))
    subplots.append(fig.add_subplot(gs[4:5, 2:6]))

    step = 0
    for i in range(len(data)):
        if data[i]["type"] == "class_boxplot":
            dimension_test = []
            for j in range(len(data[i]["data"])):
                # exclude project-specific items (individual TCGA projects and total for all TCGA projects)
                data[i]["data"][j] = data[i]["data"][j][[False if "TCGA" in data[i]["data"][j].iloc[k].loc["item"] or data[i]["data"][j].iloc[k].loc["item"] == "total" else True
                                                         for k in range(data[i]["data"][j].shape[0])]]
                dimension_test.append(data[i]["data"][j].shape[0])
            
            if len(np.unique(dimension_test)) != 1: print("< dimension error @class_boxplot")

            # initialize container for class item values, removing datatype specifiers (ID: or FEATURE:)
            # marked (<-) added / removed on 250516 to allow selection of features
            #temp_data = {item_dict[data[i]["data"][0].iloc[j].loc["item"].replace("FEATURE:", "").replace("ID:", "")]:
            #             {"class 1": [], "class 2": [], "description": data[i]["features"]["ylabels"]} for j in range(data[i]["data"][0].shape[0])} # <- removed
            temp_data = {item_dict[data[i]["data"][0].iloc[j].loc["item"].replace("FEATURE:", "").replace("ID:", "")]:
                         {"class 1": [], "class 2": [], "description": data[i]["features"]["ylabels"]} for j in range(data[i]["data"][0].shape[0])
                         if data[i]["data"][0].iloc[j].loc["item"].replace("FEATURE:", "").replace("ID:", "") in item_dict} # <- added

            # fill container with class item values for each class test
            for j in range(len(data[i]["data"])):
                # marked (<-) added / removed on 250516 to allow selection of features
                # data[i]["data"][j]["item"] = [item_dict[data[i]["data"][j].iloc[k].loc["item"].replace("FEATURE:", "").replace("ID:", "")] # <- removed from here
                #                              for k in range(data[i]["data"][j].shape[0])]
                # [temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 1"].append(data[i]["data"][j].iloc[k].loc["class 1"]) for k in range(data[i]["data"][j].shape[0])]
                # [temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 2"].append(data[i]["data"][j].iloc[k].loc["class 2"]) for k in range(data[i]["data"][j].shape[0])]
                # <- removed until here
                for k in range(data[i]["data"][j].shape[0]): # <- added from here
                    if data[i]["data"][j].iloc[k].loc["item"].replace("FEATURE:", "").replace("ID:", "") in item_dict:
                        data[i]["data"][j].at[data[i]["data"][j].index[k], "item"] = item_dict[data[i]["data"][j].iloc[k].loc["item"].replace("FEATURE:", "").replace("ID:", "")] 
                        temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 1"].append(data[i]["data"][j].iloc[k].loc["class 1"])
                        temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 2"].append(data[i]["data"][j].iloc[k].loc["class 2"])
                # <- added until here

            data[i]["features"]["xcol"]   = "description"
            data[i]["features"]["ycol"]   = ["class 1", "class 2"]
            data[i]["features"]["ylabel"] = data[i]["features"]["item"]
            subplots[i]                   = pu.box_plot(subplots[i], pd.DataFrame(temp_data[data[i]["features"]["item"]]), data[i]["features"])


        if data[i]["type"] == "class_selection":
            subplots[i] = pu.plot_class_selection(subplots[i], data[i]["data"][0], data[i]["features"])


        if data[i]["type"] == "correlation_matrix":
            for j in range(len(data[i]["data"])):
                # exclude project-specific item pairs (individual TCGA projects and total for all TCGA projects)
                data[i]["data"][j]         = data[i]["data"][j][[False if "TCGA" in data[i]["data"][j].iloc[k].loc["pair"]
                                                                 else True for k in range(data[i]["data"][j].shape[0])]]

                # modify item pair names to not contain type specifiers
                data[i]["data"][j]["pair"] = [data[i]["data"][j].iloc[k].loc["pair"].replace("FEATURE:", "").replace("ID:", "") for k in range(data[i]["data"][j].shape[0])]
                # select item pairs present in the item dictionary
                data[i]["data"][j]         = data[i]["data"][j][data[i]["data"][j]["pair"].isin(list(item_dict.keys()))]

                # replace and select item pair names based on the item dictionary
                data[i]["data"][j]["pair"] = [item_dict[data[i]["data"][j].iloc[k].loc["pair"]] if data[i]["data"][j].iloc[k].loc["pair"] in item_dict else None
                                              for k in range(data[i]["data"][j].shape[0])]
                data[i]["data"][j]         = data[i]["data"][j][~data[i]["data"][j]["pair"].isna()]                

                data[i]["data"][j]         = data[i]["data"][j].sort_values(by="pair", ascending=False)
                data[i]["data"][j].index   = data[i]["data"][j]["pair"]

            
            # create features only for one half of the data representing one layer
            data[i]["features"]["xcol"] = ["pair" for _ in range(int(len(data[i]["data"])/2))]
            # changed on 250317 from padj to p
            #data[i]["features"]["ycol"] = ["spearman-padj" for _ in range(int(len(data[i]["data"])/2))]
            data[i]["features"]["ycol"] = ["spearman-p" for _ in range(int(len(data[i]["data"])/2))]
            data[i]["features"]["zcol"] = ["spearman-r" for _ in range(int(len(data[i]["data"])/2))]

            # first layer: analysis w/o smoothing
            subplots[i] = pu.dotplot(subplots[i], data[i]["data"][0:int(len(data[i]["data"])/2)], data[i]["features"], layer=0)
 
            # second layer: analysis w smoothing
            if 1 in data[i]["features"]["layer"]: subplots[i] = pu.dotplot(subplots[i], data[i]["data"][int(len(data[i]["data"])/2)::], data[i]["features"], layer=1)
            subplots[i].grid(axis="y", color="lightgray", linewidth=3)


        if data[i]["type"] == "dendrogram":
            data[i]["data"], data[i]["features"]["projects"] = calculate_linkage(data[i]["data"][0], data[i]["features"])
            # sorted_projects are used for matching appearance of project matrix
            subplots[i], sorted_projects                     = pu.plot_dendrogram(subplots[i], data[i]["data"], data[i]["features"])


        if data[i]["type"] == "pie":
            data[i]["data"][0] = data[i]["data"][0][data[i]["data"][0][data[i]["features"]["label_col"]].isin(data[i]["features"]["xcols"])]
            subplots[i]        = pu.plot_pie(subplots[i], data[i]["data"][0], data[i]["features"])

        
        if data[i]["type"] == "project_matrix":
            for j in range(len(data[i]["data"])):
                # exclude item pairs other than individual TCGA projects and "total" for all TCGA projects
                data[i]["data"][j]         = data[i]["data"][j][[True if "TCGA" in data[i]["data"][j].iloc[k].loc["item"] else False for k in range(data[i]["data"][j].shape[0])]]
                data[i]["data"][j]["item"] = [data[i]["data"][j].iloc[k].loc["item"].replace("TCGA-", "") for k in range(data[i]["data"][j].shape[0])]
                data[i]["data"][j].index   = data[i]["data"][j]["item"]
                # sort projects according to dendrogram
                data[i]["data"][j]         = data[i]["data"][j].loc[sorted_projects]
            
            data[i]["features"]["xcol"] = ["item", "item", "item", "item"]
            data[i]["features"]["ycol"] = ["padj", "padj", "padj", "padj"]
            data[i]["features"]["zcol"] = ["statistic", "statistic", "statistic", "statistic"]
            subplots[i]                 = pu.dotplot(subplots[i], data[i]["data"], data[i]["features"])

        if i != 1 and i != 2:
            subplots[i].text(-0.1, 1.1, string.ascii_lowercase[step], transform=subplots[i].transAxes, size=14, weight='bold')
            step += 1

    plt.show()
    fig.savefig(run_dir + "\\Fig5.svg", dpi=resolution)
    return


if __name__ == '__main__':
    main()