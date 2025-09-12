from datetime import datetime
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import platform
import random
import scipy
from sklearn import metrics
from sklearn.tree import export_graphviz


def __balance__(df, params, phase):
    lower_label_index = [i for i in range(df.shape[0]) if df.iloc[i].loc[params["label_id"]] < params["decision_threshold"]] # df[df[params["label_id"]] < params["decision_threshold"]].index.tolist()
    upper_label_index = [i for i in range(df.shape[0]) if df.iloc[i].loc[params["label_id"]] >= params["decision_threshold"]] # df[df[params["label_id"]] >= params["decision_threshold"]].index.tolist()
    print("lower_label_index", len(lower_label_index), "upper_label_index", len(upper_label_index))

    if params["randomize_balancing"] == True:
        random.shuffle(lower_label_index)
        random.shuffle(upper_label_index)

    min_size    = min(len(lower_label_index), len(upper_label_index))
    label_index = lower_label_index[0:min_size] + upper_label_index[0:min_size]
    updated_df  = df.iloc[label_index]

    print("< balancing of", phase, "data reduced size from", df.shape[0], "to", len(label_index))
    return updated_df


def _find_missing_values(df, feature):
    # check for outdated placeholders
    if df[df[feature] == -1].shape[0] > 0:
        print("< potential outdated placeholder '-1' detected @", feature, "("+str(df[df[feature] == -1].shape[0])+")")
        # marked (<-) removed on 250414 because causing error with expected -1 value
        # exit() # <-
    
    if df[df[feature] == "-"].shape[0] > 0:
        print("< potential outdated placeholder '-' detected @", feature, "("+str(df[df[feature] == "-"].shape[0])+")")
        exit()
        
    return [i for i in range(df.shape[0]) if pd.isna(df.iloc[i].loc[feature]) == True]


class Forest_utils:
    def __init__(self, params, test_metrics, train_metrics):
        self.best_cohort        = None
        self.best_metric        = None
        self.full_value_cols    = None
        self.full_value_index   = None
        self.init_df            = pd.DataFrame()
        self.params             = params
        self.test_labels        = []
        self.test_metrics       = test_metrics
        self.test_preds         = []
        self.train_metrics      = train_metrics
        
        # marked (<-) added on 250411 to allow printing of gap values
        self.cohort_values      = {} # <-
        self.values             = {} # <-

        os_name = platform.system()
        self.os_sep = "//"
        if os_name == "Windows":
            self.os_sep = "\\"

        self.newdir = self.create_newdir()


    def append_preds(self, all_preds, all_index):
        if len(all_preds.shape) > 1: all_preds = np.average(all_preds, axis=1)
        self.init_df["FEATURE:prediction"] = [all_preds[i] for i in all_index]
        self.init_df.to_csv(path_or_buf=self.params["data_dir"]+self.os_sep+self.params["fname"], sep=",", index=False)
        return


    def convert_labels(self, df):
        if self.params["label_input"] == "ASE":
            if self.params["label_output"] == "nonASE":
                df["LABEL:NMD score"] = [-math.log2(1/(2*df.iloc[i].loc["LABEL:NMD score"])) if df.iloc[i].loc["LABEL:NMD score"] > 0 
                                         else -1000 for i in range(df.shape[0])]
                
                df = df[df["LABEL:NMD score"] != -1000]

            # marked (<-) was added on 2050410 to allow transformation according to Kim et al. (r = -ln2(VAF(RNA) / VAF(DNA)))
            # only True for normal data (with VAF(DNA) = 0.5)
            if self.params["label_output"] == "nonASE2":
                df["LABEL:NMD score"] = [-math.log2((1-df.iloc[i].loc["LABEL:NMD score"])/0.5) if df.iloc[i].loc["LABEL:NMD score"] > 0 
                                         else None for i in range(df.shape[0])]
                df = df[~df["LABEL:NMD score"].isna()]
            # <- ends here 
        
        # transformation of predictions by Lindeboom et al. 2019 (r = (A+B) / 2A)
        if self.params["lindeboom_output"] == "ASE":
            # marked (<-) edited on 250409 to avoid string conversion that was due to old placeholder usage ("-" instead of None)
            #df["FEATURE:lindeboom prediction"] = [str(1/(2*math.exp(-math.log(2)*float(df.iloc[i].loc["FEATURE:lindeboom prediction"])))) <- removed
            #                                      if pd.isna(df.iloc[i].loc["FEATURE:lindeboom prediction"]) == False else None for i in range(df.shape[0])] <- removed

            df["FEATURE:lindeboom prediction"] = [1/(2*math.exp(-math.log(2)*float(df.iloc[i].loc["FEATURE:lindeboom prediction"]))) # <- added
                                                  if pd.isna(df.iloc[i].loc["FEATURE:lindeboom prediction"]) == False else None for i in range(df.shape[0])] # <- added

        return df


    def create_newdir(self):
        datestr = str(datetime.now())
        index   = datestr.index(".")
        dirname = datestr[:index].replace(" ", "_").replace(":", "-")
        newdir  = self.params["target_dir"] + self.os_sep + dirname + "_" + self.params["tag"]
        if not os.path.exists(newdir): os.mkdir(newdir)
        return newdir


    def evaluate(self, labels, preds, cohort, phase):
        cat_labels = np.where(labels.iloc[:, 0] < self.params["decision_threshold"], 0, 1).tolist()
        cat_preds  = np.where(preds < self.params["decision_threshold"], 0, 1).tolist()

        accuracy   = metrics.accuracy_score(cat_labels, cat_preds)
        auroc      = metrics.roc_auc_score(cat_labels, preds)
        r          = scipy.stats.pearsonr(labels.iloc[:,0].tolist(), preds)
        r2         = r.statistic**2
        rmse       = metrics.mean_squared_error(labels, preds)

        if phase == "training":
            self.train_metrics["accuracy"].append(accuracy)
            self.train_metrics["auroc"].append(auroc)
            self.train_metrics["rmse"].append(rmse)
            self.train_metrics["r2"].append(r2)

        elif phase == "test":
            self.test_metrics["accuracy"].append(accuracy)
            self.test_metrics["auroc"].append(auroc)
            self.test_metrics["rmse"].append(rmse)
            self.test_metrics["r2"].append(r2)
            
            if len(self.test_labels) != 0:
                self.test_labels += labels.iloc[:,0].tolist()
                self.test_preds  += preds.tolist()

            else:
                self.test_labels = labels.iloc[:,0].tolist()
                self.test_preds  = preds.tolist()

            if self.best_cohort == None or self.best_metric < self.test_metrics[self.params["opt_metric"]][-1]:
                self.best_cohort = cohort
                self.best_metric = self.test_metrics[self.params["opt_metric"]][-1]

        print("phase:", phase, "cohort:", cohort, "AUROC:", auroc, "accuracy:", accuracy, "RMSE:", rmse, "R2", r2)
        return
    

    def evaluate_metrics(self, metrics, phase):
        full_r2 = None
        stats   = {key: [0 for _ in range(3)] for key in metrics}
        print(phase, "results:")

        for key in metrics:
            if key != "labels":
                stats[key][0] = np.average(np.array(metrics[key]))
                stats[key][1] = scipy.stats.sem(np.array(metrics[key]))
                stats[key][2] = np.argmax(np.array(metrics[key]))
                if phase == "training" or (phase == "test" and key != "r2"): print("  ", key, ":", round(stats[key][0], 3), "+/-", round(stats[key][1], 3))

                if key == "r2" and phase == "test":
                    full_r             = scipy.stats.pearsonr(self.test_labels, self.test_preds)
                    full_r2            = full_r.statistic**2
                    print("  ", key, ":", round(stats[key][0], 3), "+/-", round(stats[key][1], 3), "(", round(full_r2, 3), ")")
                
        return stats, full_r2
    

    def evaluate_model(self, model):
        std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)
        forest_importances = pd.Series(model.feature_importances_, index=model.feature_names_in_)
        forest_importances = forest_importances.sort_values(ascending=False)

        fig, ax = plt.subplots()
        forest_importances.plot.bar(yerr=std, ax=ax)
        ax.set_title("Feature importances using MDI")
        ax.set_ylabel("Mean decrease in impurity")
        fig.tight_layout()
        plt.show()
        return
    

    def evaluate_models(self, models):
        importances = np.zeros((0))
        for model in models:
            if importances.shape[0] > 0: importances = np.column_stack((importances, model.feature_importances_))
            else:                        importances = model.feature_importances_
            
        std              = np.std(importances, axis=1)
        mean_importances = pd.Series(np.average(importances, axis=1), index=models[0].feature_names_in_)
        mean_importances = mean_importances.sort_values(ascending=False)

        fig, ax = plt.subplots()
        mean_importances.plot.bar(yerr=std, ax=ax)
        ax.set_title("Feature importances using MDI")
        ax.set_ylabel("Mean decrease in impurity")
        fig.tight_layout()
        plt.show()

        # print to file
        forest_importances = pd.DataFrame({"average": mean_importances, "abs": [abs(mean_importance) for mean_importance in mean_importances],
                                           "sum": np.sum(importances, axis=1), "sem": [std[i]/math.sqrt(len(models)) for i in range(std.shape[0])]},
                                           index=models[0].feature_names_in_)
        forest_importances.to_csv(path_or_buf=self.newdir+self.params["os_sep"]+"forest_importances.txt", sep=",")
        return


    def extract(self, df, phase="training"):
        # balance datasets if selected
        if self.params["balance"] == True:
            df = __balance__(df, self.params, phase)

        # split dataframe into labels and features
        labels = []
        for column in df.columns:
            if self.params["label_id"] in column: labels.append(column)
        
        features = []
        for column in df.columns:
            if "FEATURE" in column: features.append(column)
    
        label_df   = df.loc[:, labels]
        feature_df = df.loc[:, features]
        return label_df, feature_df
    
       
    def full(self, df):
        labels = []
        for column in df.columns:
            if "LABEL" in column: labels.append(column)
        
        labels_df  = df.loc[:,labels]
            
        features = []
        for column in df.columns:
            if "FEATURE" in column: features.append(column)

        features_df = df.loc[:,features]

        return labels_df, features_df
    

    def load(self):
        paths = os.listdir(self.params["data_dir"])

        for path in paths:
            if self.params["fname"] == path:
                with open(self.params["data_dir"]+self.os_sep+path, 'r') as _:
                    df = pd.read_csv(self.params["data_dir"]+self.os_sep+path, delimiter=",", index_col=False)
                    # marked (<-) added / removed on 250510
                    # df = df.dropna(subset=self.params["label_id"]) # <- removed
                    if self.params["mode"] != "predict": df = df.dropna(subset=self.params["label_id"]) # <- added

                    # marked (<-) added on 250415 to remove gaps from features
                    if self.params["remove_gaps"] == True: # <-
                        df = df.dropna(subset=[col for col in df.columns if "FEATURE" in col and col != "FEATURE:rna half-life"]) # <-
                        df = df.reset_index() # <-
                
                # check if label is present
                if self.params["label_id"] not in df.columns:
                    print("< error. label", self.params["label_id"], "not found.")
                    exit()
                    
                # convert labels to specified format
                df = self.convert_labels(df)

                #if self.params["lindeboom_output"] != None: df = self.convert_labels(df)
                df           = df.reset_index(drop=True)
                self.init_df = df

                if self.params["decision_threshold"] == None: self.params["decision_threshold"] = df[self.params["label_id"]].median()
                    
                # initial modifications: replacement of missing values (with the exception of the skip feature)
                for feature in df.columns:
                    if "FEATURE" in feature:
                        missing_index = _find_missing_values(df, feature)

                        if len(missing_index) > 0:
                            # marked (<-) added/removed on 250414 to allow hyper-parameter und and model selection
                            #if self.params["mode"] == "fit": # <- removed
                            if self.params["mode"] == "fit" or self.params["mode"] == "tuning": # <- added
                                if self.params["cohortwise_misses"] == False:
                                    # marked (<-) added/removed on 250428 to allow usage of default values (specified in params)
                                    # mean = df[~df[feature].isna()][feature].mean() # <- removed
                                    if feature in self.params["default_values"]: mean = self.params["default_values"][feature] # <- added
                                    else:                                        mean = df[~df[feature].isna()][feature].mean() # <- added
                                    
                                    for missing_index_ in missing_index:
                                        df.at[df.index[missing_index_], feature] = mean

                                    print("< missing values for feature", feature, "replaced by", mean)
                                    # marked (<-) added on 250411 to allow printing of updated values
                                    self.values[feature] = mean # <-


                                if self.params["cohortwise_misses"] == True:
                                    testsize = 0
                                    # determine mean values to replace missing values for each cohort individually (241102)
                                    for cohort in range(self.params["cohorts"]):
                                        cohort_df       = df[df["ID:cohort"] == cohort+1]
                                        non_cohort_df   = df[df["ID:cohort"] != cohort+1]
                                        # marked (<-) added/removed on 250428 to allow usage of default values (specified in params)
                                        # mean            = non_cohort_df[~non_cohort_df[feature].isna()][feature].mean() # <-
                                        if feature in self.params["default_values"]: mean = self.params["default_values"][feature] # <- added
                                        else:                                        mean = non_cohort_df[~non_cohort_df[feature].isna()][feature].mean() # <- added

                                        for missing_index_ in missing_index:
                                            if missing_index_ in cohort_df.index:
                                                df.at[df.index[missing_index_], feature] = mean
                                                testsize                                += 1
                                        
                                        print("< missing values for feature", feature, "replaced by", mean, "in cohort", cohort)
                                        # marked (<-) added on 250411 to allow printing of updated values
                                        if feature in self.cohort_values: self.cohort_values[feature][cohort] = mean # <-
                                        else:                             self.cohort_values[feature] = {cohort: mean} # <-

                                    if testsize != len(missing_index):
                                        print("< error occurred during replacement of missing values ("+str(testsize)+"/"+str(len(missing_index))+")")
                                        exit()


                            else:
                                if self.params["cohortwise_misses"] == False:
                                    for missing_index_ in missing_index:
                                        # marked (<-) added/removed on 250428 to allow usage of default values (specified in params)
                                        # df.at[df.index[missing_index_], feature] = self.params["values"][feature] # <- removed
                                        if feature in self.params["default_values"]: df.at[df.index[missing_index_], feature] = self.params["default_values"][feature] # <- added
                                        else:                                        df.at[df.index[missing_index_], feature] = self.params["values"][feature] # <- added

                                # apply cohort-specific mean values to replace missing values (241102)
                                if self.params["cohortwise_misses"] == True:
                                    testsize = 0

                                    for cohort in range(self.params["cohorts"]):
                                        cohort_df       = df[df["ID:cohort"] == cohort+1]

                                        for missing_index_ in missing_index:
                                            if missing_index_ in cohort_df.index:
                                                # marked (<-) added/removed on 250428 to allow usage of default values (specified in params)
                                                # df.at[df.index[missing_index_], feature] = self.params["cohort_values"][feature][cohort] # <- removed
                                                if feature in self.params["default_values"]: df.at[df.index[missing_index_], feature] = self.params["default_values"][feature] # <- added
                                                else:                                        df.at[df.index[missing_index_], feature] = self.params["cohort_values"][feature][cohort] # <- added
                                                testsize += 1

                                    if testsize != len(missing_index):
                                        print("< error occurred during replacement of missing values ("+str(testsize)+"/"+str(len(missing_index))+")")
                                        exit()


                # print features selected for fitting or prediction
                if self.params["verbosity"] == 1:
                    print("< selected features:")
                    for column in df.columns:
                        if "FEATURE" in column: print(" ", column)

        if self.params["labels_to_class"]:
            df[self.params["label_id"]] = np.where(df[self.params["label_id"]] < self.params["decision_threshold"], 0, 1).tolist()

        return df        
    

    def print_evaluation(self, feature_names):
        file = open(self.newdir+self.os_sep+"results.txt", 'w')
        
        if self.params["mode"] == "fit":
            training_stats, _ = self.evaluate_metrics(self.train_metrics, "training")
            file.write("training results:" + "\n")
            for key in training_stats:
                file.write(key + ": " + str(training_stats[key][0]) + "+/-" + str(training_stats[key][1]) + "\n")

            file.write("\n")
        
        testing_stats, full_r2  = self.evaluate_metrics(self.test_metrics, "test")

        file.write("test results:" + "\n")
        for key in testing_stats:
            if key == "r2": file.write(key + ": " + str(testing_stats[key][0]) + "+/-" + str(testing_stats[key][1]) + " (" + str(full_r2) + ")" + "\n")
            if key != "r2": file.write(key + ": " + str(testing_stats[key][0]) + "+/-" + str(testing_stats[key][1]) + "\n")

        # marked (<-) added on 250411 to allow printing of updated values
        if len(self.cohort_values) > 0: self.params["cohort_values"] = self.cohort_values # <-
        if len(self.values) > 0:        self.params["values"]        = self.values # <-

        file.write("\n")
        file.write("hyperparameters:" + "\n")
        for key in self.params:
            file.write(key + ": " + str(self.params[key])+ "\n")

        file.write("used features: ")
        for i in range(len(feature_names)):
            file.write(feature_names[i])
            if i < len(feature_names[i])-1: file.write(",")

        file.write("\n")
        file.close()       
    

    def print_model(self, model, feature_names, cohort_id):
        pickle.dump(model, open(self.newdir+self.os_sep+"model_" + str(cohort_id) + ".pickle", "wb"))
        
        tree = model[0]
        export_graphviz(tree, out_file=self.newdir+self.os_sep+"tree.dot", feature_names=feature_names,
                        rounded=True, proportion=False, precision=2, filled=True)
    
    
    def print_preds(self, df, all_preds, all_index):
        if len(all_preds.shape) > 1: all_preds = np.average(all_preds, axis=1)

        # marked (<-) added / removed on 250630
        # df["FEATURE:prediction"] = [all_preds[i] for i in all_index] # <- removed
        # df                       = df.iloc[all_index] # <- removed

        df["FEATURE:prediction"] = [None for _ in all_index] # <- removed
        
        for i, i_ in enumerate(all_index): # <- added
            df.at[i_, "FEATURE:prediction"] = float(all_preds[i]) # <- added

        df = df[[col for col in df.columns if "ID" in col or "prediction" in col or "LABEL" in col]]
        if self.params["balance"] == False: df.to_csv(path_or_buf=self.newdir+self.params["os_sep"]+"all_preds.txt", sep=",", index=False)
        if self.params["balance"] == True:  df.to_csv(path_or_buf=self.newdir+self.params["os_sep"]+"all_balanced_preds.txt", sep=",", index=False)


    def print_roc(self, all_labels, all_preds, cohort_labels, cohort_preds):
        if self.params["balance"] == False: fname = "roc"
        if self.params["balance"] == True:  fname = "balanced_roc"

        for i in range(len(cohort_labels)):
            rocs = pd.DataFrame({"fpr": [], "tpr": []})
            cat_labels  = np.where(cohort_labels[i].iloc[:, 0] < self.params["decision_threshold"], 0, 1).tolist()
            fpr, tpr, _ = metrics.roc_curve(cat_labels, cohort_preds[i])
            rocs["fpr"] = fpr
            rocs["tpr"] = tpr
            rocs.to_csv(path_or_buf=self.newdir+self.params["os_sep"]+fname+"_"+str(i)+".txt", sep=",", index=False)

        rocs        = pd.DataFrame({"fpr": [], "tpr": []})
        cat_labels  = np.where(all_labels[:, 0] < self.params["decision_threshold"], 0, 1).tolist()

        fpr, tpr, _ = metrics.roc_curve(cat_labels, all_preds)
        rocs["fpr"] = fpr
        rocs["tpr"] = tpr
        rocs.to_csv(path_or_buf=self.newdir+self.params["os_sep"]+"all_"+fname+".txt", sep=",", index=False)
        return


    def split(self, df, cohort):
        test_df  = df[df["ID:cohort"] == cohort]
        train_df = df[df["ID:cohort"] != cohort]

        # marked (<-) added/removed on 250426 to account for missing gene ids in data from Kim et al.
        # check that identical gene ids are not present in both train_df and test_df
        if "ID:gene id" in test_df.columns:
            selected = train_df[train_df["ID:gene id"].isin(test_df["ID:gene id"])] # <- added
            
            # if train_df[train_df["ID:gene id"].isin(test_df["ID:gene id"])].shape[0] > 0: # <- removed
            if selected.shape[0] > 0 and selected[~selected["ID:gene id"].isna()].shape[0] > 0: # <- added
                print("< error. ambiguous gene ids detected.")
                # print(train_df[train_df["ID:gene id"].isin(test_df["ID:gene id"])]["ID:gene id"].tolist()) # <- removed
                print(selected["ID:gene id"].tolist()) # <- added
                exit()

        if "ID:gene symbol" in test_df.columns: # <- added
            selected = train_df[train_df["ID:gene symbol"].isin(test_df["ID:gene symbol"])] # <- added
            
            if selected.shape[0] > 0 and selected[~selected["ID:gene symbol"].isna()].shape[0] > 0: # <- added
                print("< error. ambiguous gene ids detected.") # <- added
                print(selected["ID:gene symbol"].tolist()) # <- added
                exit() # <- added
        
        return train_df, test_df