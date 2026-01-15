import os

from plot_load import *
from plot_utils import *
from plot_tools import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    data = [
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

    dims       = (1, 4)
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
    subplots.append(fig.add_subplot(gs[0, 1:4]))
    subplots.append(fig.add_subplot(gs[0, 4:7]))
    subplots.append(fig.add_subplot(gs[0, 7:10]))
    subplots.append(fig.add_subplot(gs[0, 10:]))

    step = 0
    for i in range(len(data)):
        if data[i]["type"] == "coordinates":
            subplots[i] = pu.plot_coordinates(subplots[i], data[i]["data"][0], data[i]["features"])


        if data[i]["type"] == "gene_histogram":
            subplots[i] = pu.plot_gene_histogram(subplots[i], data[i]["data"][0], data[i]["features"])

    plt.show()
    fig.savefig(run_dir + "\\FigS8.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()