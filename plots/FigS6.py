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
            "data"           : [],
            "datatype"       : ["pandas", "pandas"],
            "extensions"     : ["selection_stats", None],
            "features"       : {
                               "bar_dimension": ("10%", "10%"),
                               "bar_label"    : "NMD adaptation",
                               "cmap"         : "RdBu",
                               "scale"        : (-1, 1),
                               "xcol"         : [[None], [None]],
                               "xlabels"      : [],
                               "xmute"        : False,
                               "ymute"        : False,
                               },
            "off"            : False,
            "paths"          : [[parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_13-28-51_tcga_projectwise\2025-07-11_09-29-57_binomial"],
                                [parent_dir+r"\data\prediction_analysis_tcga\2025-07-10_11-04-37_tcga\2025-07-11_08-31-27_binomial\selection_stats.txt"]],
            "separators"     : [",", ","],
            "target_col"     : "binomial-statistic FEATURE:prediction",
            "type"           : "matrix"
            },
        ]

    dims       = (7, 4)
    resolution = 600
    run_dir    = parent_dir+r"\data\figures"

    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)

    # define figure
    fig = plt.figure(figsize=(180/25.4, 240/25.4), constrained_layout=True)
    gs = fig.add_gridspec(dims[0], 3*dims[1]+1)

    subplots = []
    subplots.append(fig.add_subplot(gs[0:7, 0:2]))
    subplots.append(fig.add_subplot(gs[0:7, 2:4]))
    subplots.append(fig.add_subplot(gs[0:7, 4:6]))
    subplots.append(fig.add_subplot(gs[7:, 0:2]))
    subplots.append(fig.add_subplot(gs[7:, 2:4]))
    subplots.append(fig.add_subplot(gs[7:, 4:6]))

    for i in range(len(data)):
        if data[i]["type"] == "matrix":
            # initialize dataframe based on sorted projects (barplot1)
            assembled_data = pd.DataFrame()
            
            for j in range(len(data[i]["data"])):                    
                names = data[i]["features"]["xlabels"][j].split("\\")
                path  = "".join(names[k]+"\\" if k < len(names)-2 else names[k] for k in range(len(names)-1))
                fname = names[-1].split("_")[-1].split(".")[0]             
                data[i]["features"]["xlabels"].append(fname)

                if path in np.ravel(data[i]["paths"]).tolist(): index = np.ravel(data[i]["paths"]).tolist().index(path)
                else:                                           index = np.ravel(data[i]["paths"]).tolist().index(data[i]["features"]["xlabels"][j])
                if index % 2 == 0:                              project = fname.replace("TCGA", "")            
                elif index == 1 or index % 2 == 1:              project = "all"

                for k in range(len(data[i]["data"][j])):
                    if data[i]["data"][j].iloc[k].loc["block id"] != "total":
                        # avoid duplicates (can occur for gene symbols)
                        if project+"_"+data[i]["data"][j].iloc[k].loc["block id"] not in assembled_data.index:
                            current_data = pd.DataFrame({"project":                        [project],
                                                        "gene symbol":                     [data[i]["data"][j].iloc[k].loc["block id"]],
                                                        data[i]["target_col"]:             [data[i]["data"][j].iloc[k].loc[data[i]["target_col"]]]},
                                                        index=[project+"_"+data[i]["data"][j].iloc[k].loc["block id"]])

                            if assembled_data.shape[0] > 0: assembled_data = pd.concat([assembled_data, current_data])
                            else:                           assembled_data = current_data

            unique_gene_symbols = sorted(np.unique(assembled_data["gene symbol"]))
            unique_projects     = sorted(np.unique(assembled_data["project"]))

            rearranged_data     = pd.DataFrame({unique_project: [assembled_data.at[unique_project+"_"+unique_gene_symbol, data[i]["target_col"]]
                                                                 if unique_project+"_"+unique_gene_symbol in assembled_data.index else 0
                                                                 for unique_gene_symbol in unique_gene_symbols] for unique_project in unique_projects},
                                                index=unique_gene_symbols)
            

            for j in rearranged_data.index:
                if (len([rearranged_data.loc[j].loc[col] for col in rearranged_data.columns if rearranged_data.loc[j].loc[col] < 0]) > 0
                    and len([rearranged_data.loc[j].loc[col] for col in rearranged_data.columns if rearranged_data.loc[j].loc[col] > 0])):
                    print("< ambiguous trend:", j)

            index = split_index(rearranged_data.shape[0], len(subplots))

            for k in range(len(index)):
                subplots[i+k] = pu.plot_matrix(subplots[i+k], rearranged_data.iloc[index[k][0]:index[k][1]], data[i]["features"])

    plt.show()
    fig.savefig(run_dir + "\\FigS6.svg", dpi=resolution, transparent=True)


if __name__ == '__main__':
    main()