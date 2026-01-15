import copy
import os

from plot_load import *
from plot_utils import *
from plot_tools import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    template1 = {
                "data"          : [],
                "datatype"      : ["pandas"],
                "extensions"    : [None],
                "features"      : {
                                  "color":      "lightgray",
                                  "edgecolor":  "dimgray",
                                  "pair":       "FEATURE:ptc_mutations2_FEATURE:prediction",
                                  "switch_axes": False,
                                  "stats":       "spearman",
                                  "xcol":        [[None]],
                                  "xlabel":      "ptc mutations",
                                  "ylabel":      "NMD susceptibility"
                                  },
                "off"           : False,
                "paths"         : [],
                "separators"    : [","],
                "type"          : "correlation"
                }
    
    template2 = {
                "data"          : [],
                "datatype"      : ["pandas"],
                "extensions"    : [None],
                "features"      : {
                                  "color":      "lightgray",
                                  "edgecolor":  "dimgray",
                                  "pair":       "FEATURE:ptc_mutations2_ID:cnv total",
                                  "switch_axes": False,
                                  "stats":       "spearman",
                                  "xcol":        [[None]],
                                  "xlabel":      "ptc mutations",
                                  "ylabel":      "CNV total"
                                  },
                "off"           : False,
                "paths"         : [],
                "separators"    : [","],
                "type"          : "correlation"
                }
    
    template3 = {
                "data"          : [],
                "datatype"      : ["pandas"],
                "extensions"    : [None],
                "features"      : {
                                  "color":      "lightgray",
                                  "edgecolor":  "dimgray",
                                  "pair":       "FEATURE:ptc_mutations2_fpkm_unstranded",
                                  "switch_axes": False,
                                  "stats":       "spearman",
                                  "xcol":        [[None]],
                                  "xlabel":      "ptc mutations",
                                  "ylabel":      "NMD activity"
                                  },
                "off"           : False,
                "paths"         : [],
                "separators"    : [","],
                "type"          : "correlation"
                }
    
    data             = [copy.deepcopy(template1) for _ in range(6)]
    data[0]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class1\stats_summary.txt"]
    data[1]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_nmd1_class2\stats_summary.txt"]
    data[2]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class1\stats_summary.txt"]
    data[3]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\ptc50_class2\stats_summary.txt"]
    data[4]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class1\stats_summary.txt"]
    data[5]["paths"] = [parent_dir+r"\data\2025-06-23_16-06-47_TCGA_NMD_targets_analysis_FPKM_exp_ccorr\class_selection_ptc50_class2\stats_summary.txt"]

    data.extend([copy.deepcopy(template2) for _ in range(6)])
    for i in range(6):
        data[6+i]["paths"] = data[i]["paths"]

    data.extend([copy.deepcopy(template3) for _ in range(6)])
    for i in range(6):
        data[12+i]["paths"] = data[i]["paths"]

    dims       = (6, 4)
    resolution = 600
    run_dir    = parent_dir+r"\data\figures"

    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4), constrained_layout=True)
    gs  = fig.add_gridspec(dims[0], 3*dims[1]+1)

    subplots = []
    for i in range(6):
        subplots.append(fig.add_subplot(gs[i, 0:4]))
        subplots.append(fig.add_subplot(gs[i, 4:8]))
        subplots.append(fig.add_subplot(gs[i, 8:]))
    

    # load data
    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)


    pu = Plot_utils()

    for i in range(len(data)):        
        if data[i]["type"] == "correlation":
            selected_data = data[i]["data"][0][data[i]["data"][0]["pair"] == data[i]["features"]["pair"]]

            if "switch_axes" not in data[i]["features"] or data[i]["features"]["switch_axes"] == False:
                data[i]["data"][0] = pd.DataFrame({data[i]["features"]["xlabel"]: json.loads(selected_data.iloc[0].loc["x"]),
                                                   data[i]["features"]["ylabel"]: json.loads(selected_data.iloc[0].loc["y"])})
                
            else:
                data[i]["data"][0] = pd.DataFrame({data[i]["features"]["xlabel"]: json.loads(selected_data.iloc[0].loc["y"]),
                                                   data[i]["features"]["ylabel"]: json.loads(selected_data.iloc[0].loc["x"])})
                
            data[i]["features"]["x"] = data[i]["features"]["xlabel"]
            data[i]["features"]["y"] = data[i]["features"]["ylabel"]
            subplots[i]              = pu.plot_correlation(subplots[i], data[i]["data"][0], data[i]["features"])

            # <- solution added on 250911 to fix wrong label sizes (rcParams are not correctly applied, although working for plot_correlation with Fig2.py)
            for tick in subplots[i].get_xticklabels():
                tick.set_fontsize(6)

            for tick in subplots[i].get_yticklabels():
                tick.set_fontsize(6)

            subplots[i].xaxis.label.set_fontsize(7)
            subplots[i].yaxis.label.set_fontsize(7) # <- ends here

        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()
    fig.savefig(run_dir + "\\FigS9.svg", dpi=resolution, transparent=False)
    return


if __name__ == '__main__':
    main()