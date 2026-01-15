import os
from lifelines import KaplanMeierFitter # <- added on 250705
from lifelines.statistics import logrank_test # <- added on 250705
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.stats as stats
from scipy.cluster import hierarchy
import string
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


class Plot_utils():
    def __init__(self):
        self.cmap = plt.get_cmap('viridis')

        rcParams["axes.labelsize"]      = 7
        rcParams['font.family']         = 'sans-serif'
        rcParams['font.sans-serif']     = ['Arial']
        rcParams['legend.fontsize']     = 7 # <- added on 250911
        rcParams['lines.markersize']    = 3 # <- added on 250911
        rcParams['svg.fonttype']        = 'none' # <- added on 250911
        rcParams["xtick.labelsize"]     = 6
        rcParams["ytick.labelsize"]     = 6

    # checked
    def bar_plot(self, ax, data, features):
        # default parameters
        if "colors" not in features:
            features["colors"]     = [[self.cmap(0.2+0.25*(2*j+k)) for k in range(len(features["labels"][j]))] for j in range(len(data))]
            features["edgecolors"] = [[self.cmap(0.1+0.2*(2*j+k)) for k in range(len(features["labels"][j]))] for j in range(len(data))]
        
        values = [[] for _ in range(len(features["labels"][0]))]

        if "line" in features:
            if "xrange" in features: ax.plot([features["xrange"][0], features["xrange"][1]], [features["line"], features["line"]], color="lightgray", linestyle="--", linewidth=2, zorder=0)
            else:                    ax.plot([0, data[0].shape[0]], [features["line"], features["line"]], color="lightgray", linestyle="--", linewidth=2, zorder=0)

        for i in range(data[0].shape[0]):
            for j in range(len(data)):
                for k in range(len(features["labels"][j])):
                    if "no_zeros" in features and features["no_zeros"] == True:
                        if data[j].iloc[i].loc[features["xcol"][j][k]] > 0:
                            if i == data[j].shape[0]-1: ax.bar(i+0.2*(2*j+k), data[j].iloc[i].loc[features["xcol"][j][k]], alpha=1,
                                                               color=features["colors"][j][k], edgecolor=features["edgecolors"][j][k], linewidth=0.5, width=0.2)
                            else:                       ax.bar(i+0.2*(2*j+k), data[j].iloc[i].loc[features["xcol"][j][k]], alpha=1,
                                                               color=features["colors"][j][k], edgecolor=features["edgecolors"][j][k], linewidth=0.5, width=0.2)

                        if pd.isna(data[j].iloc[i].loc[features["xcol"][j][k]]) == True or data[j].iloc[i].loc[features["xcol"][j][k]] == 0: # <- added on 251110
                            if i == data[j].shape[0]-1: ax.bar(i+0.2*(2*j+k), 0.0025, color="white", edgecolor=features["edgecolors"][j][k], linewidth=0.5, width=0.2)
                            else:                       ax.bar(i+0.2*(2*j+k), 0.0025, color="white", edgecolor=features["edgecolors"][j][k], linewidth=0.5, width=0.2)

                    else:
                        if i == data[j].shape[0]-1: ax.bar(i+0.2*(2*j+k), data[j].iloc[i].loc[features["xcol"][j][k]],
                                                           color=features["colors"][j][k], edgecolor=features["edgecolors"][j][k], alpha=1, linewidth=0.5, width=0.2)
                        else:                       ax.bar(i+0.2*(2*j+k), data[j].iloc[i].loc[features["xcol"][j][k]],
                                                           color=features["colors"][j][k], edgecolor=features["edgecolors"][j][k], alpha=1, linewidth=0.5, width=0.2)

                    values[k].append(data[j].iloc[i].loc[features["xcol"][j][k]])


        if "show_stats" in features and features["show_stats"] == True:
            ax.legend([features["labels"][0][i] + " (r:" + str(round(stats.pearsonr(values[i], features["xlabels"]).statistic, 4)) + " " + str(round(stats.pearsonr(values[i], features["xlabels"]).pvalue, 4))
                                                +  ", s:" + str(round(stats.spearmanr(values[i], features["xlabels"]).statistic, 4)) + " " + str(round(stats.spearmanr(values[i], features["xlabels"]).pvalue, 4)) + ")"
                       for i in range(len(features["labels"][0]))])
                    
        else:
            ax.legend(features["labels"][0])

        ax.set_ylabel(features["ylabel"])

        if "xmute" not in features or ("xmute" in features and features["xmute"] == False):
            ax.set_xticks(np.arange(data[0].shape[0]), features["xlabels"], rotation="vertical", va="center")

        else:
            ax.tick_params(axis='both', which='both', bottom=False, right=False, top=False, labelbottom=False)
        
        if "xrange" in features: ax.set_xlim(features["xrange"])
        if "yrange" in features: ax.set_ylim(features["yrange"])
        return ax


    # checked
    def box_plot(self, ax, data, features):
        # default settings
        boxtype    = "explicit"
        meanline = False; showfliers = True; showmeans = False
        if "boxtype" in features:    boxtype    = features["boxtype"]
        if "meanline" in features:   meanline   = features["meanline"]
        if "showfliers" in features: showfliers = features["showfliers"]
        if "showmeans" in features:  showmeans  = features["showmeans"]

        if type(features["ycol"]) == str: features["ycol"] = [features["ycol"]]

        for i in range(data.shape[0]):
            for j, ycol in enumerate(features["ycol"]):
                if type(data.iloc[i].loc[ycol]) == str:  y = json.loads(data.iloc[i].loc[ycol])
                if type(data.iloc[i].loc[ycol]) == list: y = data.iloc[i].loc[ycol]

                if "log_scale" in features and features["log_scale"] == True: y = [math.log10(y[j]) if y[j] > 0 else 0 for j in range(len(y))]           

                if boxtype == "compact":
                    boxplot = ax.boxplot(y, labels=[data.iloc[i].loc[features["xcol"]].replace("TCGA-", "").replace("_Solid_Tissue_Normal", " normal")], 
                                         positions=[i*len(features["ycol"])+0.5*j], showfliers=showfliers, patch_artist=True, widths=0.5, meanline=meanline, showmeans=showmeans)
                
                    for patch in zip(boxplot['boxes']):
                        patch[0].set_edgecolor("dimgray")

                        if "colors" not in features: patch[0].set_facecolor(self.cmap(float(i)/data.shape[0]))
                        else:                        patch[0].set_facecolor(mcolors.to_rgb(features["colors"][i*len(features["ycol"])+j])+(0.5,))

                    for line in zip(boxplot['caps']):
                        line[0].set_linewidth(1)
                        line[0].set_color("dimgray")

                    for line in zip(boxplot['medians']):
                        line[0].set_linewidth(1.5)
                        line[0].set_color(features["colors"][i*len(features["ycol"])+j]) # equal to box color to make it disappear
                        
                    for line in zip(boxplot['whiskers']):
                        line[0].set_linewidth(1)
                        line[0].set_color("dimgray")


                if boxtype == "explicit":
                    boxplot = ax.boxplot(y, labels=[data.iloc[i].loc[features["xcol"]].replace("TCGA-", "").replace("_Solid_Tissue_Normal", " normal")],
                                         meanline=True, notch=True, positions=[i*len(features["ycol"])+0.5*j], showfliers=showfliers, showmeans=True, patch_artist=True, widths=0.5)

                    for patch in zip(boxplot['boxes']):
                        if "colors" not in features: patch[0].set_facecolor(self.cmap(float(i)/data.shape[0]))
                        else:                        patch[0].set_facecolor(features["colors"][i*len(features["ycol"])+j])

                    for line in zip(boxplot['caps']):
                        line[0].set_linewidth(3)

                    for line in zip(boxplot['medians']):
                        line[0].set_linewidth(3)
                        line[0].set_color("dimgray")

                    for line in zip(boxplot['means']):
                        line[0].set_linestyle("--")
                        line[0].set_linewidth(3)
                        line[0].set_color("black")

                    for line in zip(boxplot['fliers']):
                        line[0].set_linewidth(1)
                        if "colors" not in features: line[0].set_markeredgecolor(self.cmap(float(i)/data.shape[0]))
                        else:                        line[0].set_markeredgecolor(features["colors"][i*len(features["ycol"])+j])

                    for line in zip(boxplot['whiskers']):
                        line[0].set_linewidth(2)

        ax.tick_params(axis='x', labelrotation=90)
        ax.locator_params(axis='y', nbins=3, min_n_ticks=3) 

        if "ylabel" in features:
            ax.set_ylabel(features["ylabel"])
        
        if "xmute" in features and features["xmute"] == True:
            ax.tick_params(axis='both', which='both', bottom=False, labelbottom=False)
        
        if "yrange" in features:
            ax.set_ylim(features["yrange"])

        return ax


    # checked (should be re-checked)
    def _dotplot(self, value, scale, category, binary=False):
        category_index = get_category_index(value, scale)

        if binary == False:
            if category_index == None: return category[0]
            if category_index != None: return category[category_index+1]

        else:
            if category_index == None: return category[1]
            if category_index != None: return category[0]


    # xcolumn defines item names, ycolumn defines size, zcolumn defines color
    # checked
    def dotplot(self, ax, data, features, layer=None):
        # default settings
        marker     = "o"
        markersize = 5

        if layer != None:
            marker     = features["layer"][layer]["marker"]
            markersize = features["layer"][layer]["markersize"]

        scale = {"binary": False, "ycol": None, "zcol": None}
        if "scale" in features: 
            for key in features["scale"]:
                scale[key] = features["scale"][key]

        # check if colors are explicit (not color map)
        colors_explicit = False
        if type(features["scale"]["colors"][0]) == str:
            colors_explicit = True
    
        ycol_max = []; ycol_min = []
        zcol_max = []; zcol_min = []

        for i in range(len(data)):
            ycol_max.append(data[i][features["ycol"][i]].max())
            ycol_min.append(data[i][features["ycol"][i]].min())
            zcol_max.append(data[i][features["zcol"][i]].max())
            zcol_min.append(data[i][features["zcol"][i]].min())

        max_y_difference = np.max(ycol_max)-np.min(ycol_min)
        max_z_difference = np.max(zcol_max)-np.min(zcol_min)

        for i in range(len(data)):
            for j in range(data[i].shape[0]):
                if scale["ycol"] == None: size_scale  = ((data[i].iloc[j].loc[features["ycol"][i]]-np.min(ycol_min)) / max_y_difference) + 1
                else:                     size_scale  = self._dotplot(data[i].iloc[j].loc[features["ycol"][i]], scale["ycol"], scale["sizes"])
                if scale["zcol"] == None: color_scale = (data[i].iloc[j].loc[features["zcol"][i]]-np.min(zcol_min)) / max_z_difference
                else:                     color_scale = self._dotplot(data[i].iloc[j].loc[features["zcol"][i]], scale["zcol"], scale["colors"], scale["binary"])

                if size_scale > 0:
                    if colors_explicit == False:
                        if features["reverse"] == False: ax.plot(i, j, color=self.cmap(color_scale*0.8), marker=marker, markerfacecolor=self.cmap(color_scale),
                                                                markersize=markersize*size_scale, alpha=0.7)
                        if features["reverse"] == True:  ax.plot(j, i, color=self.cmap(color_scale*0.8), marker=marker, markerfacecolor=self.cmap(color_scale),
                                                                markersize=markersize*size_scale, alpha=0.7)

                    if colors_explicit == True:
                        if features["reverse"] == False: ax.plot(i, j, color="black", marker=marker, markerfacecolor=color_scale,
                                                                markersize=markersize*size_scale, alpha=0.7)
                        if features["reverse"] == True:  ax.plot(j, i, color="black", marker=marker, markerfacecolor=color_scale,
                                                                markersize=markersize*size_scale, alpha=0.7)

                # if selected, plot dots with size 0 as empty dots
                elif "include_empty" in features and features["include_empty"] == True:
                    if colors_explicit == False:
                        if features["reverse"] == False: ax.plot(i, j, color=self.cmap(color_scale*0.8), fillstyle=None,
                                                                 marker=marker, markersize=markersize*scale["sizes"][1], alpha=0.5)
                        if features["reverse"] == True:  ax.plot(j, i, color=self.cmap(color_scale*0.8), fillstyle=None,
                                                                 marker=marker, markersize=markersize*scale["sizes"][1], alpha=0.5)

                    if colors_explicit == True:
                        if features["reverse"] == False: ax.plot(i, j, color=color_scale, fillstyle=None,
                                                                 marker=marker, markersize=markersize*scale["sizes"][1], alpha=0.5)
                        if features["reverse"] == True:  ax.plot(j, i, color=color_scale, fillstyle=None,
                                                                 marker=marker, markersize=markersize*scale["sizes"][1], alpha=0.5)

            if features["reverse"] == False:
                # modify x-axis
                if "xlabel" in features:  ax.set_xlabel(features["xlabel"])

                if "xticks" not in features or ("xticks" in features and features["xticks"] == True):
                    ax.set_yticks(np.arange(data[i].shape[0]), [xcol for xcol in data[i][features["xcol"][i]]], rotation="horizontal", va="center")

                else:
                    ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False, labelleft=False)
                
                # modify y-axis
                if "ylabels" in features:    ax.set_xticks(np.arange(len(data)), features["ylabels"], rotation="vertical", va="center")

                if "yrange" not in features: ax.set_ylim(-0.5, data[0].shape[0]-0.5)
                else:                        ax.set_ylim(features["yrange"][0], features["yrange"][1])

                if "ylabels" in features:    ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False)
                else:                        ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False, labelbottom=False)

            if features["reverse"] == True:
                # modify x-axis
                if "xlabel" in features: ax.set_ylabel(features["xlabel"])

                if "xticks" not in features or ("xticks" in features and features["xticks"] == True):
                    ax.set_xticks(np.arange(data[i].shape[0]), [xcol for xcol in data[i][features["xcol"][i]]], rotation="vertical", va="center")

                else:
                    ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False, labelbottom=False)
                
                # modify y-axis
                if "ylabels" in features:    ax.set_yticks(np.arange(len(data)), features["ylabels"], rotation="horizontal", va="center")
                
                if "xrange" not in features: ax.set_xlim(-0.5, data[0].shape[0]-0.5)
                else:                        ax.set_xlim(features["xrange"][0], features["xrange"][1])

                if "yrange" not in features: ax.set_ylim(-0.5, len(data)-0.5)
                else:                        ax.set_ylim(features["yrange"][0], features["yrange"][1])

                if "ylabels" in features:    ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False)
                else:                        ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False, labelleft=False)
                
        # size legend
        if scale["ycol"] != None:
            for i in range(1, len(scale["sizes"])):
                if features["reverse"] == False: ax.plot(0, len(data)-0.5*i, marker, color="white", fillstyle=None, markeredgecolor="black",
                                                         markersize=5*scale["sizes"][i], label=str(scale["ycol"][i-1]))
                if features["reverse"] == True:  ax.plot(0, len(data)-0.5*i, marker, color="white", fillstyle=None, markeredgecolor="black",
                                                         markersize=5*scale["sizes"][i], label=str(scale["ycol"][i-1]))

            if "legend_off" not in features or features["legend_off"] == False: ax.legend()

        # color bar
        if "bar_off" not in features or features["bar_off"] == False:
            ax_inset   = inset_axes(ax, width="2%", height="30%", loc="upper left")
            normalizer = mpl.colors.Normalize(vmin=features["zrange"][0], vmax=features["zrange"][1])
            cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=self.cmap), cax=ax_inset, orientation="vertical",
                                                      label=features["zlabel"], format="%.4g", alpha=0.7)
            cbar.ax.locator_params(nbins=3)

        return ax


    # checked
    def fit_linear(self, x, A, B):
        y = A+B*x
        return y


    # checked
    def forest_plot(self, ax, data, features):
        if len(features["xerror"]) == 2:
            xerr = np.transpose(data[features["xerror"]].to_numpy())

            # error type is absolute, values must be converted to match matplotlib format
            if features["error_type"] == "absolute":
                xerr[0] = np.subtract(data[features["xcol"]], xerr[0])
                xerr[1] = np.subtract(xerr[1], data[features["xcol"]])

        ax.errorbar(data[features["xcol"]], np.arange(data.shape[0]), xerr=xerr, fmt='s', color="dimgray",
                    capsize=4, capthick=2, elinewidth=2, markeredgecolor="black", markersize=8)
        ax.set_yticks(np.arange(data.shape[0]))
        ax.set_yticklabels(features["categories"])
        ax.invert_yaxis()
        ax.axvline(x=features["center"], color='grey', linestyle='--', linewidth=2)
        ax.set_xlabel(features["xlabel"])
        return ax


    # checked
    def get_edge_index(self, bin_edges, cutoff):
        edge_index = None
        for i, bin_edge in enumerate(bin_edges):
            if bin_edge < cutoff:
                edge_index = i

        return edge_index


    # checked
    def get_fit(self, x, y):
        y_fit = []
        parameters, _ = curve_fit(self.fit_linear, x, y)

        for k in range(len(x)):
            y_fit.append(self.fit_linear(x[k], parameters[0], parameters[1]))

        return y_fit
    

    # checked
    def _get_missing_values(self, text):
        print("< missing values @", text)
        exit()


    # checked
    def get_missing_values(self, data, target_cols, text=""):
        if type(data) == list:
            for i in range(len(data)):
                for target_col in target_cols:
                    if target_col in data[i].columns:
                        if data[i][data[i][target_col] == -1].shape[0] > 0:  self._get_missing_values(text+", "+target_col+", -1")
                        if data[i][data[i][target_col].isna()].shape[0] > 0: self._get_missing_values(text+", "+target_col+", None")

        else:
            for target_col in target_cols:
                if target_col in data.columns:
                    if data[data[target_col] == -1].shape[0] > 0:  self._get_missing_values(text+", "+target_col+", -1")
                    if data[data[target_col].isna()].shape[0] > 0: self._get_missing_values(text+", "+target_col+", None")
    

    def get_positions(self, data_3d, data_xd, threshold):
        pvalues_x = [[] for _ in range(3)]
        pvalues_y = [[] for _ in range(3)]

        for i in range(len(data_xd)):
            for j in range(len(data_xd[i])):
                if data_3d[i][j] > threshold:
                    if data_xd[i][j] < 0.001 and pd.isna(data_xd[i][j]) == False:
                        pvalues_x[0].append(j)
                        pvalues_y[0].append(i)

                    elif data_xd[i][j] < 0.01 and pd.isna(data_xd[i][j]) == False:
                        pvalues_x[1].append(j)
                        pvalues_y[1].append(i)

                    elif data_xd[i][j] < 0.05 and pd.isna(data_xd[i][j]) == False:
                        pvalues_x[2].append(j)
                        pvalues_y[2].append(i)

        return pvalues_x, pvalues_y


    # checked
    def get_ticks(self, x, steps=3):
        end       = round_digitwise(max(x))
        increment = round_digitwise(max(x)-min(x)) / (steps-1)
        return [end-(steps-1-i)*increment for i in range(steps)]
    

    def plot_3d(self, ax, data, features):
        data_3d = data.to_numpy()
        # default parameters
        if "cmap" not in features:     features["cmap"]         = plt.get_cmap('viridis_r')
        if "3d_range" not in features: features["3d_range"]     = (np.nanmin(data_3d), np.nanmax(data_3d))
        
        x = np.arange(len(data_3d[0]))
        y = np.arange(len(data_3d))

        contour_x, contour_y = np.meshgrid(x, y)
        ax.pcolor(contour_x, contour_y, data_3d, cmap=features["cmap"], shading='nearest',
                  vmin=features["3d_range"][0], vmax=features["3d_range"][1], edgecolors='none')
        
        ax.set_xticks(np.arange(data.shape[1]), data.columns, rotation="vertical", va="bottom", position=(-0.05, -0.05))

        normalizer = mpl.colors.Normalize(vmin=features["3d_range"][0], vmax=features["3d_range"][1])
        ax_inset   = inset_axes(ax, width="2%", height="30%", loc="upper left")
        cbar       = plt.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=features["cmap"]), cax=ax_inset, label=features["bar_label"],
                                  shrink=0.5, format="%.1f", alpha=0.7)
        cbar.ax.locator_params(nbins=6)

        if "yticks" in features: ax.set_yticks(features["yticks"], features["yticks"])
        if "ylabel" in features: ax.set_ylabel(features["ylabel"])
        
        return ax


    def plot_4d(self, ax, data_3d, features, pos, data_4d=np.array([]), data_5d=np.array([])):
        # default parameters
        if "cmap" not in features:          features["cmap"]         = plt.get_cmap('viridis_r')
        if "3d_range" not in features:      features["3d_range"]     = (np.nanmin(data_3d), np.nanmax(data_3d))
        if "4d_threshold" not in features:  features["4d_threshold"] = 0

        # switch axis off
        ax[pos[0]-1][pos[1]+1].axis("off")
        
        # fill aggregated observations subplots, x-axis (aggregated observations)
        y = np.nanmean(data_3d, 0)
        x = np.arange(len(y)) 
        ax[pos[0]-1][pos[1]].plot([0, max(x)], [1, 1], "-", color="grey", linewidth=2)
        ax[pos[0]-1][pos[1]].bar(x, y, color=features["cmap"](0.5))
        ax[pos[0]-1][pos[1]].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax[pos[0]-1][pos[1]].set_xmargin(0)

        # fill aggregated observations subplots, x-axis (aggregated observations)
        y = np.nanmean(data_3d, 1)
        x = np.arange(len(y))
        ax[pos[0]][pos[1]+1].plot([1, 1], [0, max(x)], "-", color="grey", linewidth=2)
        ax[pos[0]][pos[1]+1].barh(x, y, color=features["cmap"](0.5))
        ax[pos[0]][pos[1]+1].tick_params(axis='y', which='both', left=False, labelleft=False)
        ax[pos[0]][pos[1]+1].set_ymargin(0)

        # fill contour plot subplots
        x = np.arange(len(data_3d[0]))
        y = np.arange(len(data_3d))

        contour_x, contour_y = np.meshgrid(x, y)
        ax[pos[0]][pos[1]].pcolor(contour_x, contour_y, data_3d, cmap=features["cmap"], shading='nearest',
                                  vmin=features["3d_range"][0], vmax=features["3d_range"][1])

        if len(data_4d.shape) > 0:
            pvalues_4d_x, pvalues_4d_y = self.get_positions(data_3d, data_4d, features["4d_threshold"])
            ax[pos[0]][pos[1]].scatter(pvalues_4d_x[0], pvalues_4d_y[0], facecolors='#343835', edgecolor="#343835", s=30, linewidths=0.5)
            ax[pos[0]][pos[1]].scatter(pvalues_4d_x[1], pvalues_4d_y[1], facecolors='#343835', edgecolor="#343835", s=15, linewidths=0.5)
            ax[pos[0]][pos[1]].scatter(pvalues_4d_x[2], pvalues_4d_y[2], facecolors='#343835', edgecolor="#343835", s=5, linewidths=0.5)

        if len(data_5d.shape) > 0:
            pvalues_5d_x, pvalues_5d_y = self.get_positions(data_3d, data_5d, features["4d_threshold"])
            ax[pos[0]][pos[1]].scatter(pvalues_5d_x[0], pvalues_5d_y[0], facecolors='white', edgecolor="white", s=2, linewidths=0.5)
            ax[pos[0]][pos[1]].scatter(pvalues_5d_x[1], pvalues_5d_y[1], facecolors='white', edgecolor="white", s=1, linewidths=0.5)
            ax[pos[0]][pos[1]].scatter(pvalues_5d_x[2], pvalues_5d_y[2], facecolors='white', edgecolor="white", s=0.5, linewidths=0.5)
        
        ax[pos[0]][pos[1]].set_xticks(np.arange(len(features["xlabels"])), features["xlabels"], rotation="vertical", va="bottom", position=(-0.05, -0.05))
        ax[pos[0]][pos[1]].set_yticks(np.arange(len(features["ylabels"])), features["ylabels"], rotation="horizontal", va="bottom", position=(0, 0))

        if "title" in features:
            ax[pos[0]][pos[1]].set_title(features["title"])

        normalizer = mpl.colors.Normalize(vmin=features["3d_range"][0], vmax=features["3d_range"][1])
        ax_inset   = inset_axes(ax[pos[0]][pos[1]], width="2%", height="30%", loc="upper left")
        cbar       = plt.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=features["cmap"]), cax=ax_inset, #location="right",
                                  shrink=0.5, format="%.1f", alpha=0.7)
        cbar.ax.locator_params(nbins=6)
        return ax


    # checked
    def plot_auroc(self, ax, data, features):
        # skip first entry that holds average auroc and iterate only over cohort-wise aurocs
        for i in range(1, len(data)):
            x           = data[i][features["x"]].tolist()
            y           = data[i][features["y"]].tolist()
            color_scale = 0.2*float(i)/float(len(data)-1)
            ax.plot(x, y, '-', color=self.cmap(0.4+color_scale), linewidth=2, alpha=0.5)

        ax.plot(data[0][features["x"]], data[0][features["y"]], '-', color=self.cmap(0.3), linewidth=3, alpha=1)
        ax.set_xlabel(features["xlabel"])
        ax.set_ylabel(features["ylabel"])
        return ax


    # checked
    def plot_class_selection(self, ax, data, features):
        x        = json.loads(data.iloc[0].loc["x"])
        y        = json.loads(data.iloc[0].loc["y"])
        class1_x = json.loads(data.iloc[0].loc["class1_x"])
        class1_y = json.loads(data.iloc[0].loc["class1_y"])
        class2_x = json.loads(data.iloc[0].loc["class2_x"])
        class2_y = json.loads(data.iloc[0].loc["class2_y"])

        if "colors" not in features or "data" not in features["colors"]: ax.plot(x, y, linestyle="", marker="o", markeredgecolor=self.cmap(0.5),
                                                                                 markerfacecolor=self.cmap(0.7), markersize=2.5)
        else:                                                            ax.plot(x, y, linestyle="", marker="o", markeredgecolor=features["edgecolors"]["data"],
                                                                                 markerfacecolor=features["colors"]["data"], markersize=2.5)
        
        ax.plot(class1_x, class1_y, linestyle="", marker="o", markeredgecolor=features["colors"]["class1"], markerfacecolor=features["colors"]["class1"], markersize=2.5)
        ax.plot(class2_x, class2_y, linestyle="", marker="o", markeredgecolor=features["colors"]["class2"], markerfacecolor=features["colors"]["class2"], markersize=2.5)

        y_fit1 = self.get_fit(class1_x, class1_y)
        if "colors" not in features: ax.plot(class1_x, y_fit1, linestyle="--", linewidth=3, color=self.cmap(0.3))
        else:                        ax.plot(class1_x, y_fit1, linestyle="--", linewidth=3, color=features["colors"]["class1"], alpha=0.5)

        y_fit2 = self.get_fit(class2_x, class2_y)
        if "colors" not in features: ax.plot(class2_x, y_fit2, linestyle="--", linewidth=3, color=self.cmap(0.3))
        else:                        ax.plot(class2_x, y_fit2, linestyle="--", linewidth=3, color=features["colors"]["class2"], alpha=0.5)

        if "x_mute" not in features or ("x_mute" in features and features["x_mute"] == False):
            ax.set_xlabel(features["xlabel"])

        if "yticks" not in features or features["yticks"] != None: ax.set_yticks(self.get_ticks(y))
        ax.set_ylabel(features["ylabel"])

        if "x_mute" in features and features["x_mute"] == True:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        return ax


    # checked
    def plot_coefficients(self, ax, data, features, max_coefficients=10):
        # default settings
        if "max_coefficients" in features: max_coefficients = np.min(features["max_coefficients"], data.shape[0])
        else:                              max_coefficients = data.shape[0]

        data = data.sort_values(by=[features["x"]])
        data = data.iloc[0:max_coefficients]

        coefficient_names = [coefficient[8::] for coefficient in data.index]
        for i in range(data.shape[0]):
            ax.annotate(coefficient_names[i], (i, data.iloc[i].loc[features["x"]]), fontsize=rcParams['legend.fontsize'])
            color_scale = data.iloc[i].loc[features["z"]] / (features["zcutoff"][1]-features["zcutoff"][0])
            ax.plot(i, data.iloc[i].loc[features["x"]], 's', color=self.cmap(color_scale*0.8), markerfacecolor=self.cmap(color_scale), markersize=5, alpha=0.7)
        
        # color bar
        ax_inset   = inset_axes(ax, width="2%", height="30%", loc="upper left")
        normalizer = mpl.colors.Normalize(vmin=features["zcutoff"][0], vmax=features["zcutoff"][1])
        cbar       = plt.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=self.cmap), cax=ax_inset, orientation="vertical",
                                                        label=features["zlabel"], format="%.3f", alpha=0.7)
        cbar.ax.locator_params(nbins=3)
        
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.set_ylabel(features["ylabel"])
        ax.grid(axis='y', color='0.8')
        return ax
    

    # checked
    def plot_correlation(self, ax, data, features):
        # default settings
        if "color" not in features:     features["color"]     = self.cmap(0.7)
        if "edgecolor" not in features: features["edgecolor"] = self.cmap(0.5)
        if "linecolor" not in features: features["linecolor"] = self.cmap(0.3)


        # check for missing values
        self.get_missing_values(data, [features["x"], features["y"]], "plot_correlation")

        x = data[features["x"]]
        y = data[features["y"]]
        ax.plot(x, y, linestyle="", marker="o", markeredgecolor=features["edgecolor"], markerfacecolor=features["color"], markersize=5, alpha=0.3)

        if "stats" in features and "pearson" in features["stats"]:
            y_fit = self.get_fit(x.tolist(), y.tolist())
            ax.plot(x, y_fit, linestyle="--", linewidth=2, color=features["linecolor"])

        if "stats" in features and features["stats"] == "pearson2":
            r = stats.pearsonr(x, y)

            # <- if...else statement added on 250624
            if r.pvalue > 0.000001: 
                anchored_text = AnchoredText("r² " + str(round(math.pow(r.statistic, 2), 4)) + " (" + str('{:.2e}'.format(r.pvalue)) + ")", loc=2, prop=dict(size=rcParams['legend.fontsize']))

            else:
                anchored_text = AnchoredText("r² " + str(round(math.pow(r.statistic, 2), 4)) + " (< " + str('{:.2e}'.format(0.000001)) + ")", loc=2, prop=dict(size=rcParams['legend.fontsize']))

            ax.add_artist(anchored_text)

        if "stats" in features and features["stats"] == "spearman":
            s = stats.spearmanr(x, y)
            # <- if...else statement added on 250624
            if s.pvalue > 0.000001: 
                anchored_text = AnchoredText("s " + str(round(s.statistic, 4)) + "  (" + str('{:.2e}'.format(s.pvalue)) + ")", loc=2, prop=dict(size=rcParams['legend.fontsize']))

            else:
                anchored_text = AnchoredText("s " + str(round(s.statistic, 4)) + "  (< " + str('{:.2e}'.format(0.000001)) + ")", loc=2, prop=dict(size=rcParams['legend.fontsize']))

            ax.add_artist(anchored_text)

        if "x_mute" not in features or ("x_mute" in features and features["x_mute"] == False):
            ax.set_xlabel(features["xlabel"])

        #ax.set_yticks(self.get_ticks(y))
        ax.set_ylabel(features["ylabel"])

        if "x_mute" in features and features["x_mute"] == True:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        return ax
    

    # checked
    def plot_coordinates(self, ax, data, features):
        # marked (<-) added / removed on 250529
        model_positions   = np.arange(len(data["model predictions"])) # <- added
        model_predictions = data["model predictions"] # <- added

        if "remove_placeholders" in features and features["remove_placeholders"] == True: # <- added
            model_positions   = [i for i, pred in enumerate(data["model predictions"]) if pd.isna(pred) == False] # <- added
            model_predictions = [pred for pred in data["model predictions"] if pd.isna(pred) == False] # <- added

        # default values
        if "color" not in features: features["color"] = self.cmap(0.55)

        if "stats" in features and features["stats"] == True:
            if len(model_predictions) > 0: # <- added
                label = (features["label"] + "\n"
                        + "mean: " + str(round(np.mean(model_predictions), 4)) + "\n"
                        + "std:  " + str(round(np.std(model_predictions), 4)) + "\n") # <- added
            
        else:
            label = features["label"]

        # ax.plot(np.arange(len(data["model predictions"])), data["model predictions"], color=features["color"], alpha=0.3, label=label, zorder=0) # <- removed
        ax.plot(model_positions, model_predictions, color=features["color"], alpha=0.3, label=label, zorder=0) # <- added
        ejc_pos = 0

        for exon in data["exon predictions"]:
            if "yrange" not in features:
                # ax.plot([ejc_pos, ejc_pos], [min(data["model predictions"]), max(data["model predictions"])], "--", color="lightgray", linewidth=2) # <- removed
                ax.plot([ejc_pos, ejc_pos], [min(model_predictions), max(model_predictions)], "--", color="lightgray", linewidth=0.5) # <- added

            else:
                ax.plot([ejc_pos, ejc_pos], [features["yrange"][0], features["yrange"][1]], "--", color="lightgray", linewidth=0.5)

            ejc_pos += len(data["exon predictions"][exon])

        if pd.isna(data["real positions"]) == False:
            for exon in data["real positions"]:
                for i in range(len(data["real positions"][exon])):
                    ax.scatter(data["real positions"][exon][i], data["real predictions"][exon][i], color=features["color"], marker="o", s=3, zorder=1)
                    
        ax.legend()
        if "xlabel" in features: ax.set_xlabel(features["xlabel"])
        if "ylabel" in features: ax.set_ylabel(features["ylabel"])
        if "yrange" in features: ax.set_ylim(features["yrange"])
        if "yticks" in features: ax.set_yticks(features["yticks"])
        return ax


    # checked
    def plot_gene_histogram(self, ax, data, features):
        # default values
        if "color" not in features: features["color"] = self.cmap(0.55)
        if "bins" not in features:  features["bins"]  = 20
        if "label" not in features: features["label"] = ""
        
        # marked (<-) added / replaced on 250529
        model_predictions = data["model predictions"] # <- added
        if "remove_placeholders" in features and features["remove_placeholders"] == True: # <- added
            model_predictions = [pred for pred in data["model predictions"] if pd.isna(pred) == False] # <- added

        # hist, bin_edges = np.histogram(data["model predictions"], bins=features["bins"], density=True) # <- replaced
        hist, bin_edges = np.histogram(model_predictions, bins=features["bins"], density=True) # <- added
        ax.plot(bin_edges[0:bin_edges.shape[0]-1], hist, alpha=0.3, color=features["color"], linewidth=2, drawstyle="steps-post", label="expected")

        predictions = []
        # unclear whether None needs to be expressed by np.isnan
        if pd.isna(data["real predictions"]) == False:
            for exon in data["real predictions"]:
                for i in range(len(data["real predictions"][exon])):
                    predictions.append(data["real predictions"][exon][i])

        if len(predictions) > 0:
            #pvalue, _, _ = bootstrapping_test(data["model predictions"], np.mean(predictions), len(predictions), simulation_steps=100) # <- replaced
            pvalue, _, _ = bootstrapping_test(model_predictions, np.mean(predictions), len(predictions), simulation_steps=100) # <- added

        if "prediction_mode" in features and features["prediction_mode"] == "histogram":
            hist, bin_edges = np.histogram(predictions, bins=features["bins"], density=True)
            
            if "stats" in features and features["stats"] == True:
                ax.plot(bin_edges[0:bin_edges.shape[0]-1], hist, color=features["color"], linewidth=2, drawstyle="steps-post", label="pvalue "+str(round(pvalue, 4)))

            else:
                ax.plot(bin_edges[0:bin_edges.shape[0]-1], hist, color=features["color"], linewidth=2, drawstyle="steps-post", label=features["label"])

        else:
            for i in range(len(predictions)):
                if "stats" in features and features["stats"] == True:
                    if i == 0: ax.scatter(predictions[i], np.max(hist)/2, color=features["color"], marker="o", edgecolor="lightgray", s=20, alpha=0.3, label="pvalue "+str(round(pvalue, 4)))
                    else:      ax.scatter(predictions[i], np.max(hist)/2, color=features["color"], marker="o", edgecolor="lightgray", s=20, alpha=0.3, label="pvalue "+str(round(pvalue, 4)))

                else:
                    if i == 0: ax.scatter(predictions[i], np.max(hist)/2, color=features["color"], marker="o", edgecolor="lightgray", s=20, alpha=0.3, label=features["label"])
                    else:      ax.scatter(predictions[i], np.max(hist)/2, color=features["color"], marker="o", edgecolor="lightgray", s=20, alpha=0.3, label=features["label"])

        if "xlabel" in features: ax.set_xlabel(features["xlabel"])
        if "ylabel" in features: ax.set_ylabel(features["ylabel"])
        if "xrange" in features: ax.set_xlim(features["xrange"])
        if "xticks" in features: ax.set_xticks(features["xticks"]) # <- added on 250527
        ax.legend()
        return ax


    # checked
    def plot_dendrogram(self, ax, data, features):
        hierarchy.set_link_color_palette(['darkcyan', 'cornflowerblue', 'mediumslateblue', 'darkmagenta'])
        plot_features = hierarchy.dendrogram(data, ax=ax, color_threshold=0.1, p=20, leaf_rotation=90., leaf_font_size=7,#rcParams["xtick.labelsize"],
                                             show_contracted=True, labels=features["projects"])
        ax.tick_params(axis='both', which='both', bottom=False, left=False, right=False, top=False, labelleft=False)
        return ax, plot_features["ivl"]
    

    # marked (<-) added on 250705 (modified from survival_analysis_utils)
    def plot_kaplan_meier(self, ax, data, features, pvalue=None):
        kmf1 = KaplanMeierFitter()
        kmf1.fit(data[0][data[0].columns[0]], data[0][data[0].columns[1]])
        sf1  = kmf1.survival_function_
        ci1  = kmf1.confidence_interval_survival_function_
        ax.plot(sf1.index, sf1[sf1.columns[0]], color=features["colors"][0], label=features["labels"][0])
        ax.fill_between(ci1.index, ci1[ci1.columns[0]], ci1[ci1.columns[1]], alpha=.15)

        kmf2 = KaplanMeierFitter()
        kmf2.fit(data[1][data[1].columns[0]], data[1][data[1].columns[1]])
        sf2  = kmf2.survival_function_
        ci2  = kmf2.confidence_interval_survival_function_
        ax.plot(sf2.index, sf2[sf2.columns[0]], color=features["colors"][1], label=features["labels"][1])
        ax.fill_between(ci2.index, ci2[ci2.columns[0]], ci2[ci2.columns[1]], alpha=.15)

        ax.legend()
        ax.set_xlabel(features["xlabel"])
        ax.set_ylabel(features["ylabel"])

        if pvalue == None:
            pvalue = logrank_test(data[0][data[0].columns[0]], data[1][data[1].columns[0]],
                                event_observed_A=data[0][data[0].columns[1]], event_observed_B=data[1][data[1].columns[1]]).p_value

        if pvalue > 0.000001: anchored_text = AnchoredText("p-value " + str('{:.2e}'.format(pvalue)), loc=2)
        else:                 anchored_text = AnchoredText("p-value < " + str('{:.2e}'.format(0.000001)), loc=2)

        ax.add_artist(anchored_text)
        return ax


    # checked
    def plot_matrix(self, ax, data, features):
        if "bar_dimension" not in features: features["bar_dimension"] = ("2%", "30%")
        if "cmap" not in features:          features["cmap"]  = plt.get_cmap('viridis_r')
        if "scale" not in features:         features["scale"] = (0, 0.5)
        
        ax.imshow(data.to_numpy(), aspect='auto', cmap=features["cmap"], vmin=features["scale"][0], vmax=features["scale"][1])

        # color bar
        ax_inset   = inset_axes(ax, width=features["bar_dimension"][0], height=features["bar_dimension"][1], loc="upper left")
        normalizer = mpl.colors.Normalize(vmin=features["scale"][0], vmax=features["scale"][1])

        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=features["cmap"]), cax=ax_inset, orientation="vertical",
                                                  label=features["bar_label"], format="%.2f", alpha=0.5)
        cbar.ax.locator_params(nbins=3)

        # Adding labels to the matrix
        ax.set_yticks([i for i in range(data.shape[0])], data.index, rotation="horizontal", va="center")

        if "xmute" in features and features["xmute"] == True: ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        else:                                                 ax.set_xticks([i for i in range(data.shape[1])], data.columns, rotation="vertical", va="center")
        if "ymute" in features and features["ymute"] == True: ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        return ax


    def plot_pie(self, ax, data, features):
        ax.pie(data[features["size_col"]], colors=features["colors"], labels=["n="+str(size) for size in data[features["size_col"]]], shadow=True)
        ax.set_xlabel(features["xlabel"])
        ax.set_alpha(0.5)
        return ax


    # checked (with exception of fit)
    def plot_step_histogram(self, ax, data, features):
        # determine xcolumns
        xcols = []
        for i in range(len(features["xcol"])):
            for j in range(len(features["xcol"][i])):
                xcols.append(features["xcol"][i][j])

        # check for missing values
        self.get_missing_values(data, xcols, "plot_step_histogram")

        # default parameters
        if "bins" not in features:                bins    = 100
        else:                                     bins    = features["bins"]
        if "colors" not in features:              colors  = [self.cmap(0.3+0.25*i) for i in range(len(data))]
        else:                                     colors  = features["colors"]
        if "density" not in features:             density = False
        else:                                     density = features["density"]
        if "xrange" not in list(features.keys()): xrange  = (0, 1)
        else:                                     xrange  = features["xrange"]
        
        if "split_histogram" in features and features["split_histogram"] == True:
            all_values = []
            [all_values.extend(data[i][xcols[i]].tolist()) for i in range(len(data))]
            if "split_limits" not in features: mean = np.mean(all_values)

        all_bin_edges = []; all_hists = []; values = []
        for i in range(len(data)):
            hist, bin_edges = np.histogram(data[i][xcols[i]], bins=bins, density=density)
            values.append(data[i][xcols[i]])

            # append pseudo-element for better visualization
            bin_edges = np.append(bin_edges, 1+1/bins)
            hist      = np.append(hist, hist[-1])

            if "fill_area" in features:
                all_bin_edges.append(bin_edges)
                all_hists.append(hist)

            if "fill_area" in features and "full" in features["fill_area"] and i in features["fill_area"]["full"]:
                ax.fill_between(bin_edges[0:bin_edges.shape[0]-1], hist, alpha=0.3, color=colors[i], step="post", zorder=0.5)
  
            ax.plot(bin_edges[0:bin_edges.shape[0]-1], hist, color=colors[i], drawstyle="steps-post", label=features["labels"][i])

            if "fit" in features and features["fit"][i] == True:
                params, fit_y = fit(bin_edges[0:bin_edges.shape[0]-1], hist, function="gaussian", components=3)
                print("< params")
                print(params)

                for j in range(len(fit_y)):
                    ax.plot(bin_edges[0:bin_edges.shape[0]-1], fit_y[j], "--", color="dimgray")

            # marked (<-) added on 250522
            if "line" in features:
                for key in features["line"][i]:
                    pos = [j for j, bin_edge in enumerate(bin_edges) if bin_edge >= key][0]
                    ax.plot([bin_edges[pos], bin_edges[pos]], [0, hist[pos]], "--", color=features["line"][i][key])

        # calculates area within specified range between data 1 and 2
        if "fill_area" in features:
            start_index = [i for i, bin_edge in enumerate(all_bin_edges[0]) if bin_edge >= features["fill_area"]["range"][0]][0]
            end_index   = [i for i, bin_edge in enumerate(all_bin_edges[0]) if bin_edge >= features["fill_area"]["range"][1]][0]
            index1 = []; index2 = []

            if "auto-range" in features["fill_area"] and features["fill_area"]["auto-range"] == True:
                for i in range(start_index, end_index):
                    if all_hists[0][i] > all_hists[1][i] and (len(index1) == 0 or index1[-1] == i-1): index1.append(i)
                    elif all_hists[0][i] > all_hists[1][i]:                                           index1 = [i]
                    if all_hists[0][i] < all_hists[1][i] and (len(index2) == 0 or index2[-1] == i-1): index2.append(i)
                    elif all_hists[0][i] < all_hists[1][i]:                                           index2 = [i]
                
                if len(index1) > len(index2): index = index1
                if len(index1) < len(index2): index = index2

            else:
                index = [i for i in range(start_index, end_index)]

            ax.fill_between(bin_edges[index[0]:index[-1]], all_hists[0][index[0]:index[-1]], all_hists[1][index[0]:index[-1]], alpha=0.5,
                            color=features["fill_area"]["color"], step="post")

        if len(values) == 2:
            print("< Mann-Whitney", np.mean(values[1]) / np.mean(values[0]), stats.mannwhitneyu(values[0], values[1]).pvalue)

        if "x_mute" not in features or ("x_mute" in features and features["x_mute"] == False): ax.set_xlabel(features["xlabel"])
        ax.set_xlim(xrange)
        ax.set_ylabel(features["ylabel"])
        ax.legend()

        if "x_mute" in features and features["x_mute"] == True:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        return ax


    # <- subfunction added on 250524
    def _plot_projections(self, ax, data, features):
        if "xrange" not in list(features.keys()): xrange = (0, 100)
        else:                                     xrange = features["xrange"]

        if "bar_items" in features:
            item_count = 0
            for i in range(len(data)):
                bar_item = features["bar_items"][i]

                if bar_item in data[i]:
                    max = float(np.max([len(data[i][bar_item][j]) for j in range(len(data[i][bar_item]))]))

                    for j in range(len(data[i][bar_item])):
                        if len(data[i][bar_item][j]) > 0 and max > 0:
                            ax.bar(j, 0.6*item_count+float(len(data[i][bar_item][j]))/((len(features["bar_items"])-1)*max),
                                   alpha=0.5, color=self.cmap(0.3+0.25*i), width=0.5)
                    
                    item_count += 1

        ax.set_xlim(xrange) 
        ax.set_ylim((0, 1))
        return ax


    # checked
    def plot_projections(self, ax, data, features):
        # default parameters
        if "xrange" not in list(features.keys()): xrange = (0, 100)
        else:                                     xrange = features["xrange"]

        if "yrange" not in list(features.keys()): yrange = (0, 1)
        else:                                     yrange = features["yrange"]

        if "colors" not in features:              colors = [self.cmap(0.3+0.25*i) for i in range(len(data))]
        else:                                     colors = features["colors"]

        xcols = []
        for i in range(len(features["xcol"])):
            for j in range(len(features["xcol"][i])):
                xcols.append(features["xcol"][i][j])

        offset = np.mean(data[0][features["offset_items"][0]][xrange[0]])-np.mean(data[1][features["offset_items"][1]][xrange[0]])
        #print("offset", offset, np.mean(data[0][features["offset_items"][0]][xrange[0]]), np.mean(data[1][features["offset_items"][1]][xrange[0]]))
        exon_means1 = []
        for j in range(len(data[0][features["offset_items"][0]])):
            if "offset_correction" in features and features["offset_correction"] == True:
                exon_means1.append(np.mean(data[0][features["offset_items"][0]][j])-offset)

            else:
                exon_means1.append(np.mean(data[0][features["offset_items"][0]][j]))
        
        ax.plot(np.arange(len(exon_means1)), exon_means1, alpha=0.3, color=colors[0], linestyle="-", linewidth=3, label=features["labels"][0])

        exon_means2 = []
        for j in range(len(data[1][xcols[0]])):
            exon_means2.append(np.mean(data[1][xcols[0]][j]))

        ax.plot(np.arange(len(exon_means2)), exon_means2, color=colors[1], linestyle="-", linewidth=3, label=features["labels"][1])
        
        # marked (<-) added on 250524, content of if statement previously un-indented
        if "duplicate_axis" in features and features["duplicate_axis"] == True: # <- added
            ax2 = ax.twinx() # <- indention changed from here
            ax.set_zorder(ax2.get_zorder()+1)
            
            if "bar_items" in features:
                item_count = 0
                for i in range(len(data)):
                    bar_item = features["bar_items"][i]

                    if bar_item in data[i]:
                        max = float(np.max([len(data[i][bar_item][j]) for j in range(len(data[i][bar_item]))]))

                        for j in range(len(data[i][bar_item])):
                            if len(data[i][bar_item][j]) > 0 and max > 0:
                                #print(bar_item, "i", i, "j", j, "max", max, len(data[i][bar_item][j]))
                                ax2.plot([j, j], [0.6*item_count, 0.6*item_count+float(len(data[i][bar_item][j]))/((len(features["bar_items"])-1)*max)], "-",
                                color=self.cmap(0.3+0.25*i), linewidth=1.5, label=features["labels"][i])
                        
                        item_count += 1 

            ax2.spines["right"].set_color(self.cmap(0.3))
            ax2.set_ylabel(features["ylabel2"]) # <- until here


        if "x_mute" not in features or ("x_mute" in features and features["x_mute"] == False): ax.set_xlabel(features["xlabel"])
        ax.patch.set_visible(False) # <- shifted from previous location
        ax.set_xlim(xrange) # <- shifted from previous location
        ax.set_ylim(yrange) # <- shifted from previous location
        ax.set_ylabel(features["ylabel"])
        ax.legend()

        if "x_mute" in features and features["x_mute"] == True:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        return ax
    
   
    def xy_plot(self, ax, data, features):
        if "colors" not in features: features["colors"]         = [self.cmap(0.25+0.25*i) for i in range(len(features["x"]))]
        if "edgecolors" not in features: features["edgecolors"] = [self.cmap(0.15+0.25*i) for i in range(len(features["x"]))]

        for i in range(len(features["x"])):  
            ax.plot(data[features["x"][i]], data[features["y"][i]], label=features["labels"][i], linestyle='None',
                    marker="o", markeredgecolor=features["edgecolors"][i], markerfacecolor=features["colors"][i])
        
        if "xrange" in features: ax.set_ylim(features["xrange"])
        if "yrange" in features: ax.set_ylim(features["yrange"])
        ax.set_xlabel(features["xlabel"])
        ax.set_ylabel(features["ylabel"])

        ax.legend()
        return ax