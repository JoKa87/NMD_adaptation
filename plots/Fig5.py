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
                              "scale"       : {"binary": True, "ycol": [0.05, 0.01, 0.001], "zcol": [1], "colors": ["crimson", "royalblue"], "sizes": [0, 0.5, 0.75, 1]},
                              "tag"         : "target",
                              "xcol"        : [[None], [None], [None], [None]],
                              "xrange"      : (-0.5, 32.5),
                              "xticks"      : False,
                              "ylabels"     : ["low / \n high PTC", "NMD not mut. / \n nmd mut.", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
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
                              "ylabels"     : ["low / \n high PTC", "NMD not mut. / \n NMD mut.", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
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
                              "ylabels"     : ["low / \n high PTC", "NMD not mut. / \n NMD mut.", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],            "separators"    : ["," for _ in range(4)],
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
                              "ylabels"     : ["low / \n high PTC", "NMD not mut. / \n NMD mut.", "low / \n high frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\class_test.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\class_test.txt"],            "separators"    : ["," for _ in range(3)],
            "type"          : "class_boxplot"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas" for _ in range(12)],
            "extensions"    : [None for _ in range(12)],
            "features"      : {
                              "bar_off"     : True,
                              "layer"       : {0: {"marker": "o", "markersize": 7}},# 1: {"marker": "o", "markersize": 5}},
                              "reverse"     : False,
                              "scale"       : {"binary": True, "ycol": [0.05, 0.01, 0.001], "zcol": [0], "colors": ["royalblue", "crimson"], "sizes": [0, 0.5, 1, 1.5]},
                              "tag"         : "target",
                              "xcol"        : [[None] for _ in range(12)],
                              "ylabels"     : ["low PTC", "high PTC", "NMD not mut.", "NMD mut.", "low \n frameshift", "high \n frameshift"]
                              },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class2\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class2_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class2_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1_smoothing\stats_summary.txt",
                               parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class2_smoothing\stats_summary.txt"],
            "separators"    : ["," for _ in range(12)],
            "type"          : "correlation_matrix",
            },
            {
            "data"      : [],
            "datatype"  : ["pandas"],
            "extensions": [None, None],
            "features"  : {
                           "bins"               : 40,
                           "colors"             : ["forestgreen"],
                           "density"            : False,
                           "labels"             : ["PTC mutations"],
                           "xcol"               : [[None]],
                           "xlabel"             : "NMD adaptation score",
                           "x_mute"             : False,
                           "y_mute"             : True,
                           "xrange"             : (-0.5, 0.6),
                           "ylabel"             : "counts"
                        },
            "paths"     : [r"C:\Programming\Translational_genomics\NMD_analysis\data\prediction_analysis_tcga\2025-10-24_18-34-19_tcga_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas"],
            "extensions": [None],
            "features"  : {
                           "bins"               : 40,
                           "colors"             : ["forestgreen"],
                           "density"            : False,
                           "labels"             : ["PTC mutations"],
                           "xcol"               : [[None]],
                           "xlabel"             : "NMD adaptation score",
                           "x_mute"             : False,
                           "xrange"             : (-0.5, 0.6),
                           "ylabel"             : "counts"
                          },
            "paths"     : [r"C:\Programming\Translational_genomics\NMD_analysis\data\prediction_analysis_tcga\2025-10-24_18-34-19_tcga_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 0.1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas"],
            "extensions": [None, None],
            "features"  : {
                           "bins"               : 40,
                           "colors"             : ["forestgreen"],
                           "density"            : False,
                           "labels"             : ["PTC mutations"],
                           "xcol"               : [[None]],
                           "xlabel"             : "NMD adaptation score",
                           "x_mute"             : False,
                           "xrange"             : (-0.5, 0.6),
                           "ylabel"             : "counts"
                        },
            "paths"     : [r"C:\Programming\Translational_genomics\NMD_analysis\data\prediction_analysis_msk\2025-11-14_19-11-03_msk_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas"],
            "extensions": [None],
            "features"  : {
                           "bins"               : 30,
                           "colors"             : ["forestgreen"],
                           "density"            : False,
                           "labels"             : ["PTC mutations"],
                           "xcol"               : [[None]],
                           "xlabel"             : "NMD adaptation score",
                           "x_mute"             : False,
                           "xrange"             : (0, 0.6),
                           "ylabel"             : "counts"
                          },
            "paths"     : [r"C:\Programming\Translational_genomics\NMD_analysis\data\prediction_analysis_msk\2025-11-14_19-11-03_msk_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 0.1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram"
            },
        ]
    
    item_dict = {
                 "cnv total"                            : "copy number",
                 "escape"                               : "escape",
                 "frameshift"                           : "frameshift",
                 "frameshift_mutations"                 : "frameshift mutations",
                 #"Immune_score"                         : "immune editing",
                 "ptc_mutations"                        : "PTC mutations alt.",
                 "ptc_mutations2"                       : "PTC mutations",
                 "target"                               : "NMD targets",
                 "fpkm_unstranded"                      : "NMD activity",
                 "prediction"                           : "NMD susceptibility",
                 "escape"                                : "NMD escape",
                 #"total_hla"                            : "immune editing",
                 #"frameshift_cnv total"                 : "frameshift-copy number",
                 "ptc_mutations2_cnv total"             : "PTC mutations-copy number",
                 "ptc_mutations2_fpkm_unstranded"       : "PTC mutations-NMD activity",
                 #"ptc_mutations2_escape"                : "PTC mutations-NMD escape",
                 "ptc_mutations2_prediction"            : "PTC mutations-NMD susceptibility",
                 #"ptc_mutations2_prediction"            : "PTC mutations-NMD susceptibility",
                 #"ptc_mutations2_total_hla"             : "PTC mutations-immune editing",
                 #"ptc_mutations2_Immune_score"          : "PTC mutations-immune editing",
                 #"expression_cnv total"                 : "expression-copy number",
                 #"expression_fpkm_unstranded"           : "expression-nmd activity",
                 #"expression_frameshift"                : "expression-frameshift",
                 #"expression_Immune_score"              : "expression-immune score",
                 #"expression_ptc_mutations2"            : "basal expression-PTC mutations",
                 #"expression_target"                    : "expression-nmd targets",
                 #"fpkm_unstranded_cnv total"            : "nmd activity-copy number",
                 #"fpkm_unstranded_Immune_score"         : "nmd activity-immune score",
                 #"Immune_score_cnv total"               : "immune score-copy number",
                 #"target_cnv total"                     : "NMD targets-copy number",
                 #"target_fpkm_unstranded"               : "NMD targets-NMD activity",
                 }


    dims       = (6, 4)
    resolution = 600
    run_dir    = parent_dir+r"\data"

    pu = Plot_utils()


    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 0:1]))
    subplots.append(fig.add_subplot(gs[0, 1:2]))
    subplots.append(fig.add_subplot(gs[0, 2:3]))
    subplots.append(fig.add_subplot(gs[1, 0:4]))
    subplots.append(fig.add_subplot(gs[0, 6:12]))
    subplots.append(fig.add_subplot(gs[1, 6:12]))
    subplots.append(fig.add_subplot(gs[2, 0:3]))
    subplots.append(fig.add_subplot(gs[2, 4:7]))
    subplots.append(fig.add_subplot(gs[2, 8:11]))
    subplots.append(fig.add_subplot(gs[4:5, 5:12]))
    subplots.append(fig.add_subplot(gs[5, 0:2]))
    subplots.append(fig.add_subplot(gs[5, 3:5]))
    subplots.append(fig.add_subplot(gs[5, 6:8]))
    subplots.append(fig.add_subplot(gs[5, 9:11]))

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
            temp_data = {item_dict[data[i]["data"][0].iloc[j].loc["item"].replace("FEATURE:", "").replace("ID:", "")]:
                         {"class 1": [], "class 2": [], "description": data[i]["features"]["ylabels"]} for j in range(data[i]["data"][0].shape[0])
                         if data[i]["data"][0].iloc[j].loc["item"].replace("FEATURE:", "").replace("ID:", "") in item_dict} # <- added

            # fill container with class item values for each class test
            for j in range(len(data[i]["data"])):
                for k in range(data[i]["data"][j].shape[0]):
                    if data[i]["data"][j].iloc[k].loc["item"].replace("FEATURE:", "").replace("ID:", "") in item_dict:
                        data[i]["data"][j].at[data[i]["data"][j].index[k], "item"] = item_dict[data[i]["data"][j].iloc[k].loc["item"].replace("FEATURE:", "").replace("ID:", "")] 
                        temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 1"].append(data[i]["data"][j].iloc[k].loc["class 1"])
                        temp_data[data[i]["data"][j].iloc[k].loc["item"]]["class 2"].append(data[i]["data"][j].iloc[k].loc["class 2"])

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

        
        if data[i]["type"] == "step_histogram":
            selection_col               = list(data[i]["selection"].keys())[0]
            data[i]["data"][0]          = data[i]["data"][0][data[i]["data"][0][selection_col] <= data[i]["selection"][selection_col]][[data[i]["target_col"]]]
            print(data[i]["data"][0])
            data[i]["features"]["xcol"] = [[data[i]["target_col"]]]
            subplots[i]                 = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])

        if i != 1 and i != 2:
            subplots[i].text(-0.1, 1.1, string.ascii_lowercase[step], transform=subplots[i].transAxes, size=9, weight='bold')
            step += 1

    plt.show()
    fig.savefig(run_dir + "\\Fig5.svg", dpi=resolution)


if __name__ == '__main__':
    main()