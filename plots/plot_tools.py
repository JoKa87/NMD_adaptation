import os
from scipy.cluster.hierarchy import linkage
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


def apply_grid(fig, grid, dims, g, margin):
    subplots = []
    h_margin = int(dims[1]*g[1]*margin[1])
    v_margin = int(dims[0]*g[0]*margin[0])

    gs  = fig.add_gridspec(g[0]*dims[0]+v_margin, g[1]*dims[1]+h_margin)
    row = 0
    for slots in grid:
        col = 0
        for i, slot in enumerate(slots):
            if type(slot[0]) == tuple:
                subplots.append(fig.add_subplot(gs[int(g[0]*(slot[0][0]+row*margin[0])):int(g[0]*(slot[0][1]+row*margin[0])),
                                                   int(g[1]*(slot[1][0]+col*margin[1])):int(g[1]*(slot[1][1]+col*margin[1]))]))

                if i < len(slots)-1 and type(slots[i+1][0]) != tuple:
                    col += 1

                elif i < len(slots)-1 and (slots[i+1][1][0] != slot[1][0] or slots[i+1][1][1] != slot[1][1]):
                    col += 1

            elif pd.isna(slot[0]) == False:
                subplots.append(fig.add_subplot(gs[int(g[0]*(slot[0]+row*margin[0])):int(g[0]*(slot[0]+row*margin[0]+1)),
                                                   int(g[1]*(slot[1]+col*margin[1])):int(g[1]*(slot[1]+col*margin[1]+1))]))
                col += 1

            else:
                col += 1

        row += 1

    return fig, subplots


def calculate_linkage(data, features):
    projects = []
    values   = np.array([])

    for feature in features["features"]:
        temp_values = []
        for i in range(data.shape[0]):
            temp_values.append(data.iloc[i].loc[feature+"_mean"])
    
        if values.shape[0] > 0: values = np.column_stack((values, temp_values))
        else:                   values = np.array(temp_values)

    projects = np.array([data.iloc[i].loc["project"].replace("TCGA-", "") for i in range(data.shape[0])])

    # sum normalization
    if features["normalization_mode"] == "sum":
        values /= values.sum(axis=0)

    # min-max normalization
    if features["normalization_mode"] == "min-max":
        for i in range(values.shape[1]):
            max_val = np.max(values[i])
            min_val = np.min(values[i])
            [(values[j][i]-min_val)/(max_val-min_val) for j in range(values.shape[0])]

    # Perform hierarchical clustering with complete linkage
    linkage_matrix = linkage(values, method='complete', metric='euclidean')

    return linkage_matrix, projects


def normalize(data, target_cols):
    values = []
    for i in range(data.shape[0]):
        for target_col in target_cols:
            if type(data.iloc[i].loc[target_col]) == str:
                data.at[data.index[i], target_col] = json.loads(data.iloc[i].loc[target_col])

            values.extend(data.iloc[i].loc[target_col])

    max_value = np.max(values)
    min_value = np.min(values)

    for i in range(data.shape[0]):
        for target_col in target_cols:
            data.at[data.index[i], target_col] = [(data.iloc[i].loc[target_col][j]-min_value)/(max_value-min_value) for j in range(len(data.iloc[i].loc[target_col]))]
    
    return data


# prepare feature data for boxplot of escape/target values, yselector is the NMD prediction score column
def prepare_boxplot(data, features):
    projects     = data.drop_duplicates(subset="ID:project")["ID:project"].tolist()
    updated_data = {"project": [], "escape_values": [], "target_values": []}

    # changed on 250224
    if data[data[features["yselector"]] == "-"].shape[0] > 0:
        print("< old placeholder '-' detected.")
        exit()
    
    data = data[~data[features["yselector"]].isna()]
    data = data[~data[features["ycol"]].isna()]

    # fill new datafarame with features seelected according to yselector
    for project in projects:
        project_data = data[data["ID:project"] == project]
        updated_data["project"].append(project)
        updated_data["escape_values"].append(project_data[project_data[features["yselector"]] <= features["score_threshold"]["FEATURE:escape"][1]][features["ycol"]].tolist())
        updated_data["target_values"].append(project_data[project_data[features["yselector"]] >= features["score_threshold"]["FEATURE:target"][0]][features["ycol"]].tolist())

    updated_data["project"].append("all")
    updated_data["escape_values"].append(data[data[features["yselector"]] <= features["score_threshold"]["FEATURE:escape"][1]][features["ycol"]].tolist())
    updated_data["target_values"].append(data[data[features["yselector"]] >= features["score_threshold"]["FEATURE:target"][0]][features["ycol"]].tolist())
    return pd.DataFrame(updated_data)