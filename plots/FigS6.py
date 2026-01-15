import os

from plot_load import *
from plot_utils import *
from plot_tools import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    boxplot_colors = []
    [boxplot_colors.extend(["crimson", "royalblue"]) for _ in range(34)]

    # params #
    data = [
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                               "colors":         [["dimgray"]],
                               "edgecolors":     [["black"]],
                               "labels":         [["all ptcs"]],
                               "significance":   0.01,
                               "xcol":           [[None], [None]],
                               "xmute":          True,
                               "xrange":         (-0.25, 33.4),
                               "ylabel":         "PTCs \n per sample"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-03_13-17-59_tcga_blocks\tcga_blocks_outer_block_stats.txt"]],
            "separators"    : [","],
            "type"          : "barplot1"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas", "pandas", "pandas", "pandas", "pandas"],
            "extensions"    : ["selection_stats", None, "selection_stats", None, "selection_stats", None],
            "features"      : {
                               "colors":         [["crimson", "lightcoral", "darkred"]],
                               "edgecolors":     [["black", "black", "black"]],
                               "labels":         [["nonsense+frameshift", "nonsense", "frameshift"]],
                               "line":           0,
                               "no_zeros":       False,
                               "xcol":           [[None], [None], [None], [None], [None], [None]],
                               "xmute":          True,
                               "xrange":         (-0.25, 33.4),
                               "xlabels":        [],
                               "ylabel":         "NMD adaptation",
                               "yrange":         (-0.22, 0.17)
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_19-29-11_tcga_mw_projectwise\2025-11-11_11-05-16_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_18-37-17_tcga_mw\2025-11-11_10-43-56_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_20-33-03_tcga_mw_nonsense_projectwise\2025-11-11_11-05-47_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_19-35-15_tcga_mw_nonsense\2025-11-11_10-46-29_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_20-51-18_tcga_mw_frameshift_projectwise\2025-11-11_11-06-14_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_19-58-09_tcga_mw_frameshift\2025-11-11_10-49-58_binomial\selection_stats.txt"]],
            "separators"    : [",", ",", ",", ",", ",", ","],
            "type"          : "barplot2"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "Blues",
                              "tag"         : "escape",
                              "xcol"        : [[None]],
                              "xmute"       : True,
                              "ymute"       : True,
                              "ycol"        : "block id"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_18-37-17_tcga_mw\2025-11-11_10-43-56_binomial\selection_stats.txt"]],
            "separators"    : [","],
            "type"          : "matrix"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "twilight",
                              "scale"       : (0, 0.5),
                              "tag"         : "escape",
                              "xcol"        : [[None], [None]],
                              "xmute"       :  True,
                              "ycol"        : "block id"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_18-37-17_tcga_mw\2025-11-11_10-43-56_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\TCGA.PanCancer.onco.genes.OncoVar.tsv"]],
            "separators"    : [",", "\t"],
            "target_cols"   : ["OncoKB", "2020Rule", "CTAT"],
            "type"          : "matrix2"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "Reds",
                              "tag"         : "target",
                              "xcol"        : [[None]],
                              "ycol"        : "block id",
                              "ymute"       : True
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_18-37-17_tcga_mw\2025-11-11_10-43-56_binomial\selection_stats.txt"]],
            "separators"    : [","],
            "type"          : "matrix"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "twilight",
                              "scale"       : (0, 0.5),
                              "tag"         : "target",
                              "xcol"        : [[None], [None]],
                              "xmute"       :  True,
                              "ycol"        : "block id"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-11-10_18-37-17_tcga_mw\2025-11-11_10-43-56_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\TCGA.PanCancer.onco.genes.OncoVar.tsv"]],
            "separators"    : [",", "\t"],
            "target_cols"   : ["OncoKB", "2020Rule", "CTAT"],
            "type"          : "matrix2"
            }
        ]

    dims       = (6, 4)
    resolution = 600
    run_dir    = parent_dir+r"\data\figures"

    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)

    #fig, subplots = apply_grid(fig, grid, dims, g, margin)
    gs = fig.add_gridspec(dims[0], 3*dims[1]+1, wspace=2)

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 1:]))
    subplots.append(fig.add_subplot(gs[1, 1:]))
    subplots.append(fig.add_subplot(gs[2:4, 1:]))
    subplots.append(fig.add_subplot(gs[2:4, 0:1]))
    subplots.append(fig.add_subplot(gs[4:6, 1:]))
    subplots.append(fig.add_subplot(gs[4:6, 0:1]))

    step = 0
    for i in range(len(data)):
        if data[i]["type"] == "barplot1":
            # initialize dictionary with ptc mutations and shared variants
            ptc_mutations = {}

            # calculate average ptc / shared ptc count per sample
            for j in range(len(data[i]["data"])):
                ptc_mutations[data[i]["features"]["labels"][0][j]+"_cases"] = data[i]["data"][j]["cases"].tolist() # no. of cases per project
                ptc_mutations["project"]                                    = data[i]["data"][j]["ID:project"].tolist()

                # same reference is used for both labels (currently, only samples with >= 1 ptc mutation are considered)
                # ratio = avg. ptc count per case per project
                ratios                                             = [data[i]["data"][j].iloc[k].loc["count"] / data[i]["data"][0].iloc[k].loc["cases"]
                                                                      if data[i]["data"][j].iloc[k].loc["cases"] > 0 else 0 for k in range(data[i]["data"][j].shape[0])]
                ptc_mutations[data[i]["features"]["labels"][0][j]] = ratios

                # add average values
                ptc_mutations[data[i]["features"]["labels"][0][j]+"_cases"].append(np.mean(ptc_mutations[data[i]["features"]["labels"][0][j]+"_cases"]))
                ptc_mutations["project"].append("all")
                ptc_mutations[data[i]["features"]["labels"][0][j]].append(np.mean(ratios))

            ptc_mutations                  = pd.DataFrame(ptc_mutations).sort_values(by=data[i]["features"]["labels"][0][0])
            data[i]["features"]["xcol"]    = data[i]["features"]["labels"]
            data[i]["features"]["xlabels"] = ptc_mutations["project"].tolist()

            if len(data[i]["features"]["labels"][0]) == 1:
                subplots[i] = pu.bar_plot(subplots[i], [ptc_mutations[[data[i]["features"]["labels"][0][0]]]], data[i]["features"])

            elif len(data[i]["features"]["labels"][0]) == 2:
                subplots[i] = pu.bar_plot(subplots[i], [ptc_mutations[[data[i]["features"]["labels"][0][0], data[i]["features"]["labels"][0][1]]]], data[i]["features"]) # <- added

            else:
                print("< unexpected label size")
                exit()


        if data[i]["type"] == "barplot2" or data[i]["type"] == "barplot3":
            # initialize dataframe based on sorted projects (barplot1)
            temp_data = pd.DataFrame({**{label:    [0 for _ in ptc_mutations["project"]] for label in data[i]["features"]["labels"][0]},
                                         "escape": [0 for _ in ptc_mutations["project"]],
                                         "target": [0 for _ in ptc_mutations["project"]]}, index=ptc_mutations["project"])
            
            xlabels = []
            for j in range(len(data[i]["data"])):
                selected_data = data[i]["data"][j]
                
                names = data[i]["features"]["xlabels"][j].split("\\")
                path  = "".join(names[k]+"\\" if k < len(names)-2 else names[k] for k in range(len(names)-1))
                fname = names[-1].split("_")[-1].split(".")[0]             
                data[i]["features"]["xlabels"].append(fname)


                if path in np.ravel(data[i]["paths"]).tolist(): index = np.ravel(data[i]["paths"]).tolist().index(path)
                # create entries for combined data
                else:                                           index = np.ravel(data[i]["paths"]).tolist().index(data[i]["features"]["xlabels"][j])

                # check whether gene summary is last entry ('total')
                if selected_data.iloc[-1].loc["block id"] == "total":
                    if data[i]["type"] == "barplot2":
                        if index % 2 == 0:
                            temp_data.at[fname, data[i]["features"]["labels"][0][(int(index/2))]] = -1*(selected_data.iloc[-1].loc["binomial-statistic FEATURE:prediction"])
                        
                        if index == 1 or index % 2 == 1:
                            temp_data.at["all", data[i]["features"]["labels"][0][(int(index/2))]] = -1*(selected_data.iloc[-1].loc["binomial-statistic FEATURE:prediction"])

                    if data[i]["type"] == "barplot3":
                        if index == 0: temp_data.at[fname, "escape"] = selected_data.iloc[-1].loc["fishers exact escape-statistic FEATURE:prediction"]
                        if index == 1: temp_data.at["all", "escape"] = selected_data.iloc[-1].loc["fishers exact escape-statistic FEATURE:prediction"]
                        if index == 2: temp_data.at[fname, "target"] = selected_data.iloc[-1].loc["fishers exact target-statistic FEATURE:prediction"]
                        if index == 3: temp_data.at["all", "target"] = selected_data.iloc[-1].loc["fishers exact target-statistic FEATURE:prediction"]

            data[i]["features"]["xlabels"] = xlabels
            data[i]["features"]["xcol"]    = data[i]["features"]["labels"]

            if data[i]["type"] == "barplot2": subplots[i] = pu.bar_plot(subplots[i], [temp_data[data[i]["features"]["xcol"][0]]], data[i]["features"])
            if data[i]["type"] == "barplot3": subplots[i] = pu.bar_plot(subplots[i], [temp_data[["escape", "target"]]], data[i]["features"])
            subplots[i].set_yticks([-0.2, -0.1, 0, 0.1])


        if data[i]["type"] == "matrix" or data[i]["type"] == "matrix2":
            # aggregated value contains 
            if data[i]["features"]["tag"] == "escape":
                selected_data = data[i]["data"][0][[True if data[i]["data"][0].iloc[j].loc["binomial-statistic FEATURE:prediction"] > 0 else False
                                                    for j in range(data[i]["data"][0].shape[0])]]
                selected_data = selected_data[selected_data["block id"] != "total"]
                selected_data = selected_data.sort_values(by=data[i]["features"]["tag"], ascending=False).iloc[0:15]
            
            if data[i]["features"]["tag"] == "target":
                selected_data = data[i]["data"][0][[True if data[i]["data"][0].iloc[j].loc["binomial-statistic FEATURE:prediction"] < 0 else False
                                                    for j in range(data[i]["data"][0].shape[0])]]
                selected_data = selected_data[selected_data["block id"] != "total"]
                selected_data = selected_data.sort_values(by=data[i]["features"]["tag"], ascending=False).iloc[0:15]

            selected_data.index = [selected_data.iloc[i].loc["block id"] for i in range(selected_data.shape[0])]
            selected_cols       = [col for col in selected_data.columns if data[i]["features"]["tag"] in col
                                    and "relative" not in col and "TCGA" in col]
            selected_data       = selected_data[selected_cols]

            # add mean value
            selected_data["all"] = selected_data.mean(axis=1)
            selected_data        = selected_data.rename(columns={col: col.split("_")[1] for col in selected_cols})

            # sort data according to ptc mutations
            data[i]["data"][0]   = selected_data[ptc_mutations["project"]]

            # check correctness of sorting
            test = [data[i]["data"][0].columns[j] for j in range(data[i]["data"][0].shape[1]) if data[i]["data"][0].columns[j] != ptc_mutations.iloc[j].loc["project"]]
            
            if len(test) > 0:
                print("< wrong sorting of projects.")
                exit()

            # normalize data based on ptc mutations (escape/target per average ptc count per case per project)
            ptc_mutations.index = ptc_mutations["project"]

            for col in data[i]["data"][0].columns:
                data[i]["data"][0][col] = [data[i]["data"][0].iloc[j].loc[col] / ptc_mutations.loc[col].loc["all ptcs_cases"] for j in range(data[i]["data"][0].shape[0])]

            data[i]["data"][0] = data[i]["data"][0].rename(columns={col: col.replace("TCGA-", "") for col in data[i]["data"][0].columns})

            if data[i]["type"] == "matrix":
                subplots[i] = pu.plot_matrix(subplots[i], data[i]["data"][0], data[i]["features"])

            if data[i]["type"] == "matrix2":
                data[i]["data"][1].index = data[i]["data"][1]["Gene_symbol"]
                data[i]["data"][0]       = pd.concat(data[i]["data"], axis=1)
                data[i]["data"][0]       = data[i]["data"][0].loc[selected_data.index]
                data[i]["data"][0]       = data[i]["data"][0].drop(columns=[col for col in data[i]["data"][0].columns if col not in data[i]["target_cols"]])
                data[i]["data"][0]       = data[i]["data"][0].replace({np.nan: 0, "Y": 0.3, "N": 0, "Pssible_Oncogene": 0.4, "Pssible_TSG": 0.2, "Oncogene": 0.4, "Oncogene/TSG": 0.3, "TSG": 0.2}).astype(float)
                subplots[i] = pu.plot_matrix(subplots[i], data[i]["data"][0], data[i]["features"])

        if data[i]["type"] != "matrix2":
            subplots[i].text(-0.1, 1.1, string.ascii_lowercase[step], transform=subplots[i].transAxes, size=9, weight='bold')
            step += 1

    plt.show()
    fig.savefig(run_dir + "\\FigS6.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()