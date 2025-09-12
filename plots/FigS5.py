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
                               "yrange":         (-0.02, 0.15)
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_15-44-58_tcga_mo_projectwise\2025-07-11_09-31-29_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-52-18_tcga_mo\2025-07-11_08-57-55_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_18-08-29_tcga_mo_nonsense_projectwise\2025-07-11_09-31-58_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_12-24-16_tcga_mo_nonsense\2025-07-11_08-58-13_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_18-29-39_tcga_mo_frameshift_projectwise\2025-07-11_09-32-27_binomial"],
                               [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_13-28-27_tcga_mo_frameshift\2025-07-11_09-29-23_binomial\selection_stats.txt"]],
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
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-52-18_tcga_mo\2025-07-11_08-57-55_binomial\selection_stats.txt"]],
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
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-52-18_tcga_mo\2025-07-11_08-57-55_binomial\selection_stats.txt"],
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
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-52-18_tcga_mo\2025-07-11_08-57-55_binomial\selection_stats.txt"]],
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
            "paths"         : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-52-18_tcga_mo\2025-07-11_08-57-55_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\TCGA.PanCancer.onco.genes.OncoVar.tsv"]],
            "separators"    : [",", "\t"],
            "target_cols"   : ["OncoKB", "2020Rule", "CTAT"],
            "type"          : "matrix2"
            },
            {
            "data"          : [],
            "datatype"      : ["json"],
            "extensions"    : [None],
            "features"      : {
                               "color"              : "royalblue",
                               "label"              : "APC",
                               "remove_placeholders": True,
                               "xcol"               : [[None]],
                               "xlabel"             : "sequence coordinate",
                               "ylabel"             : "NMD efficiency",
                               "yrange"             : (0.5, 0.8),
                               "yticks"             : [0.5, 0.6, 0.7, 0.8]
                               },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_12-26-54_tcga_traces\APC.json"],
            "separators"    : [","],
            "to_json"       : True,
            "type"          : "coordinates"
            },
            {
            "data"          : [],
            "datatype"      : ["json"],
            "extensions"    : [None],
            "features"      : {
                               "color"              : "royalblue",
                               "label"              : "observed",
                               "prediction_mode"    : "histogram",
                               "remove_placeholders": True,
                               "xcol"               : [[None]],
                               "xlabel"             : "NMD efficiency",
                               "xrange"             : (0.5, 0.8),
                               "ylabel"             : "rel. count"
                               },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_12-26-54_tcga_traces\APC.json"],
            "separators"    : [","],
            "to_json"       : True,
            "type"          : "gene_histogram"
            },
            {
            "data"          : [],
            "datatype"      : ["json"],
            "extensions"    : [None],
            "features"      : {
                               "color"              : "crimson",
                               "label"              : "TP53",
                               "remove_placeholders": True,
                               "xcol"               : [[None]],
                               "xlabel"             : "sequence coordinate",
                               "ylabel"             : "NMD efficiency",
                               "yrange"             : (0.5, 0.8),
                               "yticks"             : [0.5, 0.6, 0.7, 0.8],
                               },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_12-26-54_tcga_traces\TP53.json"],
            "separators"    : [","],
            "to_json"       : True,
            "type"          : "coordinates"
            },
            {
            "data"          : [],
            "datatype"      : ["json"],
            "extensions"    : [None],
            "features"      : {
                               "color"              : "crimson",
                               "label"              : "observed",
                               "prediction_mode"    : "histogram",
                               "remove_placeholders": True,
                               "xcol"               : [[None]],
                               "xlabel"             : "NMD efficiency",
                               "xrange"             : (0.5, 0.8),
                               "ylabel"             : "rel. count"
                               },
            "off"           : False,
            "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_12-26-54_tcga_traces\TP53.json"],
            "separators"    : [","],
            "to_json"       : True,
            "type"          : "gene_histogram"
            },
        ]

    dims            = (7, 4)
    resolution      = 600
    run_dir        = parent_dir+r"\data\figures"

    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs = fig.add_gridspec(dims[0], 3*dims[1]+1, wspace=2)

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 1:]))
    subplots.append(fig.add_subplot(gs[1, 1:]))
    subplots.append(fig.add_subplot(gs[2:4, 1:]))
    subplots.append(fig.add_subplot(gs[2:4, 0:1]))
    subplots.append(fig.add_subplot(gs[4:6, 1:]))
    subplots.append(fig.add_subplot(gs[4:6, 0:1]))
    subplots.append(fig.add_subplot(gs[6, 1:4]))
    subplots.append(fig.add_subplot(gs[6, 4:7]))
    subplots.append(fig.add_subplot(gs[6, 7:10]))
    subplots.append(fig.add_subplot(gs[6, 10:]))

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

            # marked (<-) added / removed on 250523
            # subplots[i] = pu.bar_plot(subplots[i], [ptc_mutations[[data[i]["features"]["labels"][0][0], data[i]["features"]["labels"][0][1]]]], data[i]["features"]) # <- removed
            if len(data[i]["features"]["labels"][0]) == 1: # <- added
                subplots[i] = pu.bar_plot(subplots[i], [ptc_mutations[[data[i]["features"]["labels"][0][0]]]], data[i]["features"]) # <- added

            elif len(data[i]["features"]["labels"][0]) == 2: # <- added
                subplots[i] = pu.bar_plot(subplots[i], [ptc_mutations[[data[i]["features"]["labels"][0][0], data[i]["features"]["labels"][0][1]]]], data[i]["features"]) # <- added

            else: # <- added
                print("< unexpected label size") # <- added
                exit() # <- added


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
                        #if index == 0: temp_data.at[fname, "relative"] = (selected_data.iloc[-1].loc["avg. real prediction"] / selected_data.iloc[-1].loc["avg. model prediction"])-1
                        #if index == 1: temp_data.at["all", "relative"] = (selected_data.iloc[-1].loc["avg. real prediction"] / selected_data.iloc[-1].loc["avg. model prediction"])-1

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


        if data[i]["type"] == "coordinates":
            subplots[i] = pu.plot_coordinates(subplots[i], data[i]["data"][0], data[i]["features"])


        if data[i]["type"] == "gene_histogram":
            subplots[i] = pu.plot_gene_histogram(subplots[i], data[i]["data"][0], data[i]["features"])


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
    fig.savefig(run_dir + "\\FigS4_1.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()