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
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "Blues",
                              "tag"         : "escape",
                              "xcol"        : [[None], [None]],
                              "xmute"       : True,
                              "ymute"       : True,
                              "ycol"        : "block id"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-21-27_msk\2025-07-10_15-46-32_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\msk_chord_survival_analysis.txt"]],
            "separators"    : [",", ","],
            "type"          : "matrix"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "twilight",
                              "tag"         : "escape",
                              "xcol"        : [[None], [None]],
                              "xmute"       :  True,
                              "ycol"        : "block id"
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-21-27_msk\2025-07-10_15-46-32_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\TCGA.PanCancer.onco.genes.OncoVar.tsv"]],
            "separators"    : [",", "\t"],
            "target_cols"   : ["OncoKB", "2020Rule", "CTAT"],
            "type"          : "matrix2"
            },
            {
            "data"          : [],
            "datatype"      : ["pandas", "pandas"],
            "extensions"    : [None, None],
            "features"      : {
                              "bar_label"   : "sample \n fraction",
                              "cmap"        : "Reds",
                              "tag"         : "target",
                              "xcol"        : [[None], [None]],
                              "ycol"        : "block id",
                              "ymute"       : True
                               },
            "off"           : False,
            "paths"         : [[parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-21-27_msk\2025-07-10_15-46-32_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\msk_chord_survival_analysis.txt"]],
            "separators"    : [",", ","],
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
            "paths"         : [[parent_dir+r"\data\prediction_analysis_msk\2025-07-10_10-21-27_msk\2025-07-10_15-46-32_binomial\selection_stats.txt"],
                               [parent_dir+r"\data\TCGA.PanCancer.onco.genes.OncoVar.tsv"]],
            "separators"    : [",", "\t"],
            "target_cols"   : ["OncoKB", "2020Rule", "CTAT"],
            "type"          : "matrix2"
            },
        ]

    dims       = (7, 4)
    resolution = 600
    run_dir    = parent_dir+r"\data\figures"

    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs  = fig.add_gridspec(dims[0], 3*dims[1]+1)

    subplots = []
    subplots.append(fig.add_subplot(gs[0:2, 1:4]))
    subplots.append(fig.add_subplot(gs[0:2, 0:1]))
    subplots.append(fig.add_subplot(gs[2:4, 1:4]))
    subplots.append(fig.add_subplot(gs[2:4, 0:1]))
    for i in range(3):
        subplots.append(fig.add_subplot(gs[4+i, :]))

    step = 0
    for i in range(len(data)):
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
                                    and "relative" not in col and "Cancer" in col]
            selected_data       = selected_data[selected_cols]

            # add mean value
            selected_data["all"] = selected_data.mean(axis=1)
            selected_data        = selected_data.rename(columns={col: col.split("_")[1] for col in selected_cols})
            data[i]["data"][0]   = selected_data            


            if data[i]["type"] == "matrix":
                # determine cancer type-specific counts
                cancer_types   = np.unique(data[i]["data"][1]["ID:CANCER_TYPE"])
                patient_counts = {**{cancer_type: data[i]["data"][1][data[i]["data"][1]["ID:CANCER_TYPE"] == cancer_type].shape[0] for cancer_type in cancer_types},
                                  **{"all": data[i]["data"][1].shape[0]}}
                                
                for col in data[i]["data"][0].columns:
                    data[i]["data"][0][col] = [data[i]["data"][0].iloc[j].loc[col] / patient_counts[col] for j in range(data[i]["data"][0].shape[0])]

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
    fig.savefig(run_dir + "\\FigS8.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()