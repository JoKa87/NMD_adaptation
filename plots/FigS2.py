import os

from plot_load import *
from plot_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    data = [{
             "data"          : [],
             "datatype"      : ["pandas"],
             "extensions"    : [None],
             "features"      : {
                                "labels"     : ["nonsense"],
                                "x"          : ["sequence coordinate"],
                                "xcol"       : [[None]],
                                "y"          : ["nonsense"],
                                "xlabel"     : "sequence coordinate",
                                "ylabel"     : "norm. count"
                               },
             "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\stop_codon_distribution.txt"],
             "separators"    : [","],
             "type"          : "xy_plot"
            },
            {
             "data"          : [],
             "datatype"      : ["pandas"],
             "extensions"    : [None],
             "features"      : {
                                "labels"     : ["frameshift (+1)"],
                                "x"          : ["sequence coordinate"],
                                "xcol"       : [[None]],
                                "y"          : ["frameshift_1"],
                                "xlabel"     : "sequence coordinate",
                                "ylabel"     : "norm. count"
                               },
             "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\stop_codon_distribution.txt"],
             "separators"    : [","],
             "type"          : "xy_plot"
            },
            {
             "data"          : [],
             "datatype"      : ["pandas"],
             "extensions"    : [None],
             "features"      : {
                                "labels"     : ["frameshift (+2)"],
                                "x"          : ["sequence coordinate"],
                                "xcol"       : [[None]],
                                "y"          : ["frameshift_2"],
                                "xlabel"     : "sequence coordinate",
                                "ylabel"     : "norm. count"
                               },
             "paths"         : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\stop_codon_distribution.txt"],
             "separators"    : [","],
             "type"          : "xy_plot"
            },
            {
            "data"           : [],
            "datatype"       : ["pandas"],
            "extensions"     : [None],
            "features"       : {
                               "bar_label"   : "pvalue",
                               "cmap"        : "viridis_r",
                                "scale"       : (0, 0.1),
                               "xcol"        : [[None]],
                               "xmute"       : False,
                               "ymute"       : False,
                               },
            "off"            : False,
            "paths"          : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\codon_efficiency_comparison.txt"],
            "separators"     : [","],
            "targets"        : ["ks", "mw"],
            "type"           : "matrix"
            },
            {
            "data"           : [],
            "datatype"       : ["pandas"],
            "extensions"     : [None],
            "features"       : {
                                "bar_label"   : "pvalue",
                                "cmap"        : "viridis_r",
                                "scale"       : (0, 0.1),
                                "xcol"        : [[None]],
                                "xmute"       : False,
                                "ymute"       : False,
                               },
            "off"            : False,
            "paths"          : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\codon_efficiency_comparison_non_redundant.txt"],
            "separators"     : [","],
            "targets"        : ["ks", "mw"],
            "type"           : "matrix"
            },
            {
            "data"           : [],
            "datatype"       : ["pandas"],
            "extensions"     : [None],
            "features"       : {
                                "bar_label"   : "norm. count",
                                "cmap"        : "viridis_r",
                                "xcol"        : [[None]],
                                "ylabel"      : "sequence position / %",
                                "yticks"      : [0, 25, 50, 75, 100],
                               },
            "off"            : False,
            "paths"          : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-09_17-24-03_tcga_stop_codons\codon_distribution.txt"],
            "separators"     : [","],
            "type"           : "3d_plot"
            }
        ]

    dims       = [7, 4]
    resolution = 600
    run_dir    = parent_dir+r"\data\figures"
    data       = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]


    pu = Plot_utils()

    # loading section
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs = fig.add_gridspec(dims[0], 3*dims[1], wspace=8)

    subplots = []
    subplots.append(fig.add_subplot(gs[0:2, 0:4]))
    subplots.append(fig.add_subplot(gs[0:2, 4:8]))
    subplots.append(fig.add_subplot(gs[0:2, 8:12]))
    subplots.append(fig.add_subplot(gs[2, 0:12]))
    subplots.append(fig.add_subplot(gs[3, 0:12]))
    subplots.append(fig.add_subplot(gs[4:6, 0:12]))

    for i in range(len(data)):
        if data[i]["type"] == "3d_plot":
            for col in data[i]["data"][0].columns:
                max_value               = data[i]["data"][0][col].max()
                data[i]["data"][0][col] = data[i]["data"][0][col].div(max_value)

            subplots[i] = pu.plot_3d(subplots[i], data[i]["data"][0], data[i]["features"])

        if data[i]["type"] == "matrix":
            cols               = data[i]["data"][0][data[i]["data"][0].columns[0]]
            data[i]["data"][0] = data[i]["data"][0][data[i]["targets"]].transpose()
            data[i]["data"][0] = data[i]["data"][0].rename(columns={data[i]["data"][0].columns[j]: cols[j] for j in range(data[i]["data"][0].shape[1])})
            subplots[i]        = pu.plot_matrix(subplots[i], data[i]["data"][0], data[i]["features"])

        if data[i]["type"] == "xy_plot":
            data[i]["data"][0]["sequence coordinate"] = np.arange(data[i]["data"][0].shape[0])
            
            for col in data[i]["features"]["y"]:
                max_value               = data[i]["data"][0][col].max()
                data[i]["data"][0][col] = [value/max_value for value in data[i]["data"][0][col]]

            subplots[i] = pu.xy_plot(subplots[i], data[i]["data"][0].iloc[1:data[i]["data"][0].shape[0]-1], data[i]["features"])
        
        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()

    fig.savefig(run_dir + "\\FigS2.svg", dpi=resolution)
    return


if __name__ == '__main__':
    main()