import os
import statsmodels.api as sm

from plot_load import *
from plot_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    boxplot_colors = []
    [boxplot_colors.extend(["royalblue", "crimson"]) for _ in range(5)]

    # params #
    data = [
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "Breast Cancer",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "Colorectal Cancer",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "Non-Small Cell Lung Cancer",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "Pancreatic Cancer",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "Prostate Cancer",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas", "pandas"],
            "extensions": [None, None],
            "features"  : {
                           "colors"        : ["royalblue", "crimson"],
                           "labels"        : ["lower NMD susceptibility", "higher NMD susceptibility"],
                           "xcol"          : [[None], [None]],
                           "xlabel"        : "timeline / months",
                           "ylabel"        : "survival probability"
                          },
            "km_filter" : {"FEATURE:prediction": None},
            "label_ids" : ["LABEL:OS_STATUS", "LABEL:OS_MONTHS"],
            "paths"     : [parent_dir+r"\data\msk_chord_survival_analysis_kaplan_meier.txt",
                           parent_dir+r"\data\msk_chord_survival_analysis.txt"],
            "separators": [",", ","],
            "target_col": "total",
            "type"      : "kaplan_meier"
            },
            {
            "data"      : [],
            "datatype"  : ["pandas"],
            "extensions": [None, None],
            "features"  : {
                           "density"            : True,
                           "labels"             : ["PTC mutations"],
                           "line"               : [{20: "black"}],
                           "xcol"               : [[None]],
                           "xlabel"             : "PTC mutations",
                           "x_mute"             : False,
                           "xrange"             : (0, 200),
                           "ylabel"             : "rel. counts"
                        },
            "paths"     : [parent_dir+r"\data\msk_chord_variants.txt"],
            "separators": [","],
            "target_col1": {"ID:cancer type": "Colorectal Cancer"},
            "target_col2": "FEATURE:ptc_mutations",
            "type"      : "step_histogram2"
            },
            {
            "data"      : [],
            "datatype"  : ["json"],
            "extensions": [None],
            "features"  : {
                           "density"            : True,
                           "labels"             : ["NMD susceptibilities of expected PTCs", "NMD susceptibilities of observed PTCs"],
                           "xcol"               : [["values"], ["FEATURE:prediction_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "x_mute"             : False,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                        },
            "paths"     : [parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-27-29_msk_colorectal_low_ptc\msk_colorectal_low_ptc_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["json"],
            "extensions": [None],
            "features"  : {
                           "density"            : True,
                           "labels"             : ["NMD susceptibilities of expected PTCs", "NMD susceptibilities of observed PTCs"],
                           "xcol"               : [["values"], ["FEATURE:prediction_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "x_mute"             : False,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                        },
            "paths"     : [parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-49-34_msk_colorectal_high_ptc\msk_colorectal_high_ptc_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            }
        ]

    dims           = (4, 4)
    resolution     = 600
    run_dir        = parent_dir+r"\data\figures"
    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]


    pu = Plot_utils()

    # loading section
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0, 0:4]))
    subplots.append(fig.add_subplot(gs[0, 4:8]))
    subplots.append(fig.add_subplot(gs[0, 8:12]))
    subplots.append(fig.add_subplot(gs[1, 0:4]))
    subplots.append(fig.add_subplot(gs[1, 4:8]))
    subplots.append(fig.add_subplot(gs[1, 8:12]))
    subplots.append(fig.add_subplot(gs[2, 0:4]))
    subplots.append(fig.add_subplot(gs[2, 4:8]))
    subplots.append(fig.add_subplot(gs[2, 8:12]))

    pvalues = {"cancer type": [], "pvalue": [], "counts1": [], "counts2": []}

    for i in range(len(data)):
        if data[i]["type"] == "kaplan_meier":
            selected_stats = data[i]["data"][0][data[i]["data"][0]["class"] == data[i]["target_col"]]
            cutoff         = selected_stats.iloc[0].loc["cutoff"]
            pvalue         = selected_stats.iloc[0].loc["pvalue"]
            if data[i]["target_col"] != "total": selected_data = data[i]["data"][1][data[i]["data"][1]["ID:CANCER_TYPE"] == data[i]["target_col"]]
            else:                                selected_data = data[i]["data"][1]

            selected_data1 = pd.DataFrame({data[i]["label_ids"][1]: selected_data[selected_data[list(data[i]["km_filter"].keys())[0]] <= cutoff][data[i]["label_ids"][1]],
                                           data[i]["label_ids"][0]: selected_data[selected_data[list(data[i]["km_filter"].keys())[0]] <= cutoff][data[i]["label_ids"][0]]})
            selected_data2 = pd.DataFrame({data[i]["label_ids"][1]: selected_data[selected_data[list(data[i]["km_filter"].keys())[0]] > cutoff][data[i]["label_ids"][1]],
                                           data[i]["label_ids"][0]: selected_data[selected_data[list(data[i]["km_filter"].keys())[0]] > cutoff][data[i]["label_ids"][0]]})

            pvalues["cancer type"].append(data[i]["target_col"])
            pvalues["pvalue"].append(pvalue)
            pvalues["counts1"].append(selected_data1.shape[0])
            pvalues["counts2"].append(selected_data2.shape[0])

            subplots[i] = pu.plot_kaplan_meier(subplots[i], [selected_data1, selected_data2], data[i]["features"], pvalue=pvalue)


        if data[i]["type"] == "step_histogram":
            subplots[i] = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])


        if data[i]["type"] == "step_histogram2":
            target_col                  = list(data[i]["target_col1"].keys())[0]
            data[i]["data"][0]          = data[i]["data"][0][data[i]["data"][0][target_col] == data[i]["target_col1"][target_col]]
            data[i]["features"]["xcol"] = [[data[i]["target_col2"]]]
            subplots[i] = pu.plot_step_histogram(subplots[i], data[i]["data"], data[i]["features"])


        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    pvalues["pvalue"] = [*sm.stats.fdrcorrection(pvalues["pvalue"][0:len(pvalues["pvalue"])-1], alpha=0.05)[1], pvalues["pvalue"][len(pvalues["pvalue"])-1]]
    print("< FDR-corrected pvalues")
    print(pd.DataFrame(pvalues))

    if np.sum(pvalues["counts1"][0:len(pvalues)-1]) != pvalues["counts1"][-1] or np.sum(pvalues["counts2"][0:len(pvalues)-1]) != pvalues["counts2"][-1]:
        print("< dimension error occurred.")

    plt.show()

    fig.savefig(run_dir + "\\Fig6.svg", dpi=resolution, transparent=True)
    return


if __name__ == '__main__':
    main()