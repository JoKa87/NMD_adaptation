import os

from plot_load import *
from plot_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    data = [
            {
            "data"      : [],
            "datatype"  : ["json", "json"],
            "extensions": [None, None],
            "features"  : {
                           "density"            : True,
                           "fit"                : [False, False],
                           "fill_area"          : {
                                                   "color"      : "crimson",
                                                   "range"      : (0.69, 0.77),
                                                  },
                           "labels"             : ["expected NMD susceptibilities", "NMD susceptibilities of observed PTCs"],
                           "xcol"               : [["values"], ["FEATURE:prediction_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "x_mute"             : False,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                        },
            "paths"     : [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-10-35_cptac3\cptac3_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
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
                               "no_zeros":       True,
                               "xcol":           [[None], [None], [None], [None], [None], [None]],
                               "xmute":          True,
                               "xrange":         (-0.25, 6),
                               "xlabels":        [],
                               "ylabel":         "NMD adaptation",
                               "yrange":         (-0.01, 0.06)
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_16-36-22_cptac3_projectwise\2025-11-15_14-16-31_binomial"],
                               [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-10-35_cptac3\2025-11-15_11-23-50_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_16-42-49_cptac3_nonsense_projectwise\2025-11-15_14-17-40_binomial"],
                               [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-11-45_cptac3_nonsense\2025-11-15_11-24-41_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_17-05-45_cptac3_frameshift_projectwise\2025-11-15_14-20-44_binomial"],
                               [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-46-16_cptac3_frameshift\2025-11-15_11-25-26_binomial\selection_stats.txt"]],
            "separators"    : [",", ",", ",", ",", ",", ","],
            "type"          : "barplot2"
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
            "paths"     : [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-05-59_cptac3_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram2"
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
                           "xrange"             : (0, 0.6),
                           "ylabel"             : "counts"
                          },
            "paths"     : [parent_dir+r"\data\prediction_analysis_cptac3\2025-11-14_15-05-59_cptac3_patientwise\assembled_selection_stats.txt"],
            "separators": [","],
            "selection" : {"binomial padj": 0.1},
            "target_col": "binomial statistic",
            "type"      : "step_histogram2"
            }
        ]

    dims       = [4, 4]
    resolution = 600
    run_dir    = parent_dir+r"\data"
    data       = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]


    pu = Plot_utils()

    # loading section
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    plt.subplots_adjust(wspace=0.30, hspace=0.50)

    gs = fig.add_gridspec(4*dims[0], 2*dims[1])

    # <- subplots in the current form adapted to usage of duplicate_axis=False
    subplots = []
    subplots.append(fig.add_subplot(gs[0:4, 0:4]))
    subplots.append(fig.add_subplot(gs[0:4, 4:8]))
    subplots.append(fig.add_subplot(gs[4:8, 0:4]))
    subplots.append(fig.add_subplot(gs[4:8, 4:8]))

    for i in range(len(data)):
        if data[i]["type"] == "barplot2":
            # targets with >= 100 PTC mutations
            targets = ["Adenomas and Adenocarcinomas", "Ductal and Lobular Neoplasms", "Epithelial Neoplasms, NOS", "Gliomas", "Squamous Cell Neoplasms", "all"]
            temp_data = pd.DataFrame({"relative": [0 for _ in targets],
                                      "escape": [0 for _ in targets],
                                      "target": [0 for _ in targets]}, index=targets)

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
                    if index % 2 == 0:
                        temp_data.at[fname, data[i]["features"]["labels"][0][(int(index/2))]] = -1*(selected_data.iloc[-1].loc["binomial-statistic FEATURE:prediction"])
                    
                    if index == 1 or index % 2 == 1:
                        temp_data.at["all", data[i]["features"]["labels"][0][(int(index/2))]] = -1*(selected_data.iloc[-1].loc["binomial-statistic FEATURE:prediction"])
                    print(fname, -1*(selected_data.iloc[-1].loc["binomial-statistic FEATURE:prediction"]), selected_data.iloc[-1].loc["binomial-pvalue FEATURE:prediction"])

            data[i]["features"]["xlabels"] = xlabels
            data[i]["features"]["xcol"]    = data[i]["features"]["labels"]

            subplots[i] = pu.bar_plot(subplots[i], [temp_data[data[i]["features"]["xcol"][0]]], data[i]["features"])
            subplots[i].set_yticks([0, 0.025, 0.05])


        if data[i]["type"] == "step_histogram":
            subplots[i] = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])
        
        
        if data[i]["type"] == "step_histogram2":
            selection_col               = list(data[i]["selection"].keys())[0]
            data[i]["data"][0]          = data[i]["data"][0][data[i]["data"][0][selection_col] <= data[i]["selection"][selection_col]][[data[i]["target_col"]]]
            print(data[i]["data"][0])
            data[i]["features"]["xcol"] = [[data[i]["target_col"]]]
            subplots[i]                 = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])

        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()

    fig.savefig(run_dir + "\\FigS12.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()