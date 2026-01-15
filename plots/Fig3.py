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
                           "x_mute"             : True,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                        },
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json", "json"],
            "extensions": [None, None, None],
            "features"  : {
                           "bar_items"          : [None, "FEATURE:ptc cds position_projected_mean_hist"],
                           "duplicate_axis"     : False,
                           "labels"             : ["NMDetectiveA", "NMDelphi"],
                           "offset_correction"  : True,
                           "offset_items"       : ["projected_mean_hist", "projected_mean_hist"],
                           "xcol"               : [["projected_mean_hist", "FEATURE:ptc cds position_projected_mean_hist"]],
                           "xlabel"             : "sequence position / %",
                           "x_mute"             : True,
                           "ylabel"             : "NMD susceptibility",
                           "ylabel2"            : "rel. counts",
                           "yrange"             : (0.5, 0.8)
                          },
            "off"       : False,
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_13-10-54_tcga_lindeboom\tcga_lindeboom_stats.json",
                           parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "to_json"   : True,
            "type"      : "projections"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json"],
            "extensions": [None, None],
            "features"  : {
                           "density"            : True,
                           "labels"             : ["expected NMD susceptibilities", "NMD susceptibilities of observed PTCs"],
                           "xcol"               : [["first_exon_values"], ["FEATURE:prediction_first_exon_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "x_mute"             : True,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                          },
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json", "json"],
            "extensions": [None, None, None],
            "features"  : {
                           "bar_items"          : [None, "FEATURE:ptc cds position_projected_first_exon_mean_hist"],
                           "duplicate_axis"     : False,
                           "labels"             : ["NMDetectiveA", "NMDelphi"],
                           "offset_correction"  : True,
                           "offset_items"       : ["projected_first_exon_mean_hist", "projected_first_exon_mean_hist"],
                           "xcol"               : [["projected_first_exon_mean_hist", "FEATURE:ptc cds position_projected_first_exon_mean_hist"]],
                           "xlabel"             : "sequence position / %",
                           "ylabel"             : "NMD susceptibility",
                           "ylabel2"            : "rel. counts",
                           "x_mute"             : True,
                           "yrange"             : (0.55, 0.65)
                           },
            "off"       : False,
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_13-10-54_tcga_lindeboom\tcga_lindeboom_stats.json",
                           parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "to_json"   : True,
            "type"      : "projections"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json"],
            "extensions": [None, None],
            "features"  : {
                           "density"            : True,
                           "fill_area"          : {
                                                   "color"      : "crimson",
                                                   "range"      : (0.69, 0.77),
                                                  },
                           "labels"             : ["expected NMD susceptibilities", "NMD susceptibilities of observed PTCs"],
                           "xcol"               : [["middle_exon_values"], ["FEATURE:prediction_middle_exon_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "x_mute"             : True,
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                          },
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json", "json"],
            "extensions": [None, None, None],
            "features"  : {
                           "bar_items"          : [None, "FEATURE:ptc cds position_projected_middle_exon_mean_hist"],
                           "duplicate_axis"     : False,
                           "labels"             : ["NMDetectiveA", "NMDelphi"],
                           "offset_correction"  : True,
                           "offset_items"       : ["projected_middle_exon_mean_hist", "projected_middle_exon_mean_hist"],
                           "xcol"               : [["projected_middle_exon_mean_hist", "FEATURE:ptc cds position_projected_middle_exon_mean_hist"]],
                           "xlabel"             : "sequence position / %",
                           "ylabel"             : "NMD susceptibility",
                           "ylabel2"            : "rel. counts",
                           "x_mute"             : True,
                           "yrange"             : (0.65, 0.67)
                          },
            "off"       : False,
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_13-10-54_tcga_lindeboom\tcga_lindeboom_stats.json",
                           parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "to_json"   : True,
            "type"      : "projections"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json"],
            "extensions": [None, None],
            "features"  : {
                           "density"            : True,
                           "labels"             : ["expected NMD susceptibilities", "NMD susceptibilities of observed PTCs"], 
                           "xcol"               : [["last_exon_values"], ["FEATURE:prediction_last_exon_values"]],
                           "xlabel"             : "NMD susceptibility",
                           "xrange"             : (0.5, 0.8),
                           "ylabel"             : "rel. counts"
                          },
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "type"      : "step_histogram"
            },
            {
            "data"      : [],
            "datatype"  : ["json", "json", "json"],
            "extensions": [None, None, None],
            "features"  : {
                           "bar_items"          : [None, "FEATURE:ptc cds position_projected_last_exon_mean_hist"],
                           "duplicate_axis"     : False,
                           "labels"             : ["NMDetectiveA", "NMDelphi"],
                           "offset_correction"  : True,
                           "offset_items"       : ["projected_last_exon_mean_hist", "projected_last_exon_mean_hist"],
                           "xcol"               : [["projected_last_exon_mean_hist", "FEATURE:ptc cds position_projected_last_exon_mean_hist"]],
                           "xlabel"             : "sequence position / %",
                           "ylabel"             : "NMD susceptibility",
                           "ylabel2"            : "rel. counts",
                           "yrange"             : (0.52, 0.57)
                          },
            "off"       : False,
            "paths"     : [parent_dir+r"\data\prediction_analysis_tcga\2025-07-11_13-10-54_tcga_lindeboom\tcga_lindeboom_stats.json",
                           parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\tcga_stats.json"],
            "separators": [","],
            "to_json"   : True,
            "type"      : "projections"
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

    gs = fig.add_gridspec(4*dims[0], 2*dims[1]) # , wspace=0.7, hspace=3.5

    # <- subplots in the current form adapted to usage of duplicate_axis=False
    subplots = []
    for i in range(dims[0]):
        subplots.append(fig.add_subplot(gs[4*i:4*(i+1), 0:4]))
        subplots.append(fig.add_subplot(gs[4*i:(4*i+1), 4:8]))
        subplots.append(fig.add_subplot(gs[(4*i+1):4*(i+1), 4:8]))


    step = 0; step_size = 0
    for i in range(len(data)):
        if data[i]["type"] == "projections":
            if "duplicate_axis" in data[i]["features"] and data[i]["features"]["duplicate_axis"] == True:
                subplots[step] = pu.plot_projections(subplots[step], data[i]["data"], data[i]["features"])
                subplots[step].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[step].transAxes, size=9, weight='bold')
                step_size      = 1

            else:
                subplots[step]   = pu._plot_projections(subplots[step], data[i]["data"], data[i]["features"])
                subplots[step+1] = pu.plot_projections(subplots[step+1], data[i]["data"], data[i]["features"])
                subplots[step+1].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[step+1].transAxes, size=9, weight='bold') 
                step_size        = 2

        if data[i]["type"] == "step_histogram":
            subplots[step] = pu.plot_step_histogram(subplots[step], data[i]["data"], data[i]["features"])
            subplots[step].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[step].transAxes, size=9, weight='bold')
            step_size      = 1
        
        step += step_size

    plt.show()

    fig.savefig(run_dir + "\\Fig3.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()