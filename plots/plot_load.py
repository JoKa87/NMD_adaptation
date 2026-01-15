import json
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


def check_extensions(extensions, it):
    if type(extensions) != list:
        return extensions
    
    if type(extensions) == list and it < len(extensions)-1:
        return extensions[it]

    else:
        print("< dimension error @plot_load, check_extensions for extensions:", extensions, "iterator:", it)
        exit()


def check_item(item, it):
    if type(item) != list:
        return item
    
    if type(item) == list:
        if it < len(item)-1:
            return item[it]
        
        else:
            print("< dimension error @plot_load, check_item for item:", item, "iterator:", it)
            exit()


def create_load_map(current_data_map, load_map, path, datatype, separator, apply_all_xcols):
    if path not in load_map:
        load_map[path] = {"data": None, "datatype": datatype, "separator": separator}

    current_data_map.append({"apply_all_xcols": apply_all_xcols, "path": path, "datatype": datatype})
    return current_data_map, load_map


def create_maps(data):
    data_map = []
    load_map = {}

    for i in range(len(data)):
        current_data_map = []
        apply_all_xcols  = False
        updated_xcol     = data[i]["features"]["xcol"]

        if "apply_all_xcols" in data[i]["features"] and data[i]["features"]["apply_all_xcols"] == True: apply_all_xcols  = True
        if type(data[i]["paths"][0]) != list:                                                           data[i]["paths"] = [data[i]["paths"]]

        for j in range(len(data[i]["paths"])):
            for k in range(len(data[i]["paths"][j])):
                datatype  = check_item(data[i]["datatype"][j], k)
                separator = check_item(data[i]["separators"][j], k)
                extension = check_extensions(data[i]["extensions"][j], k)

                if extension == None:
                    current_data_map, load_map = create_load_map(current_data_map, load_map, data[i]["paths"][j][k], datatype, separator, apply_all_xcols)

                else:
                    paths      = os.listdir(data[i]["paths"][j][k])
                    path_index = 0
                    for path in paths:
                        if extension in path:
                            # added so that "[None]"" tag that activates loading of the full dataset can be appended to match the number of input files
                            if data[i]["features"]["xcol"][j][k] == None and path_index > 0: updated_xcol.append([None])
                            current_data_map, load_map = create_load_map(current_data_map, load_map, data[i]["paths"][j][k]+"\\"+path, datatype, separator, apply_all_xcols)
                            path_index                += 1

        data[i]["features"]["xcol"] = updated_xcol
        data_map.append(current_data_map)
    
    #print(json.dumps(load_map, indent=4))
    
    #for i in range(len(data_map)):
    #    print(i, data_map[i])
        
    return data_map, load_map


def extract(load_map, path, cols):
    if cols[0] != None: data = {col: [] for col in cols} # copy specified columns only
    else:               data = {col: [] for col in load_map[path]["data"].columns} # copy all columns

    for col in data:
        if col in load_map[path]["data"]: data[col] = load_map[path]["data"][col] # if type(load_map[data_map[i]]["data"]) == dict:
        else:                             print("<", col, "not found in", path)

    return pd.DataFrame(data)


def _find_misses(test_df, col):
    if test_df[test_df[col] == -1].shape[0] > 0:
        print("<", test_df[test_df[col] == -1].shape[0], "times '-1' was found in", col+".")
        exit()

    if test_df[test_df[col] == "-"].shape[0] > 0:
        print("<", test_df[test_df[col] == "-"].shape[0], "times '-' was found in", col+".")
        exit()

    if test_df[test_df[col].isna()].shape[0] > 0:
        print("<", test_df[test_df[col].isna()].shape[0], "times 'None' was found in", col+".")


def find_misses(data, features):
    # get all value column names
    cols = []
    if "xcol" in features: cols.extend(np.ravel(features["xcol"]).tolist())
    if "ycol" in features: cols.extend(np.ravel(features["ycol"]).tolist())
    cols = np.unique([col for col in cols if pd.isna(col) == False])

    for col in cols:
        if col in data:
            test_df = pd.DataFrame({col: data[col]})

            if type(test_df.iloc[0].loc[col]) == list:
                for i in range(test_df.shape[0]):
                    _find_misses(pd.DataFrame({col: test_df.iloc[i].loc[col]}), col)

            if type(test_df.iloc[0].loc[col]) != list:
                _find_misses(pd.DataFrame({col: test_df[col]}), col)


# loading section
def load(data):
    data_map, load_map = create_maps(data)

    bar = IncrementalBar(set_bar("loading data"), max=len(load_map))
    for path in load_map:
        if load_map[path]["datatype"] == "pandas":
            load_map[path]["data"] = pd.read_csv(path, sep=load_map[path]["separator"])

        if load_map[path]["datatype"] == "json":
            with open(path, "r") as file:
                load_map[path]["data"] = json.load(file)

        bar.next()
    bar.finish()

    for i in range(len(data)):
        for j in range(len(data_map[i])):
            # keep dataframe or convert json to dataframe
            #print("i", i, "j", j, data_map[i][j]["path"])
            if "to_json" not in data[i] or data[i]["to_json"] == False:
                if data_map[i][j]["datatype"] == "pandas":
                    data[i]["data"].append(extract(load_map, data_map[i][j]["path"], data[i]["features"]["xcol"][j]))

                elif len(data[i]["features"]["xcol"]) != len(data_map[i]) or data_map[i][j]["apply_all_xcols"] == True:
                    for k in range(len(data[i]["features"]["xcol"])):
                        data[i]["data"].append(extract(load_map, data_map[i][j]["path"], data[i]["features"]["xcol"][k]))

                else:
                    #print("i", i, "j", j, data_map[i][j]["path"], data[i]["features"]["xcol"][j])
                    data[i]["data"].append(extract(load_map, data_map[i][j]["path"], data[i]["features"]["xcol"][j]))

                #find_misses(data[i]["data"][-1], data[i]["features"])

            # keep json
            if "to_json" in data[i] and data[i]["to_json"] == True:
                data[i]["data"].append(load_map[data_map[i][j]["path"]]["data"])
                find_misses(data[i]["data"][-1], data[i]["features"])

                
        if "xlabels" in data[i]["features"] and len(data[i]["features"]["xlabels"]) == 0:
            data[i]["features"]["xlabels"] = [data_map[i][j]["path"] for j in range(len(data_map[i]))] #get_label_names(data_map[i])

    return data
