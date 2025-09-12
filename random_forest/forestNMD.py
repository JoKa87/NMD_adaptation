# Importing the libraries
import numpy as np # for array operations
# scikit-learn modules
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
# marked (<-) added on 250414 to facilitate feature selection
from sklearn.inspection import permutation_importance # <-
# marked (<-) newly added on 250414 for hyper-parameter tuning
from sklearn.model_selection import GridSearchCV # <-
from sklearn import tree
from sklearn.tree import DecisionTreeRegressor

from create_genome_predictions import *
from forestNMD_utils import *


parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def __init_model__(mode="fast"):
    if mode == "fast_class":
        model = RandomForestClassifier(
                                       criterion="gini",
                                       n_estimators=100,
                                       bootstrap=True,
                                       max_features="sqrt",
                                       min_samples_leaf=2,
                                       min_samples_split=2,
                                       n_jobs=8,
                                       random_state=0
                                      )
        
    if mode == "slow_class":
        model = RandomForestClassifier(
                                       criterion="gini",
                                       n_estimators=10000,
                                       bootstrap=True,
                                       max_features="sqrt",
                                       min_samples_leaf=50,
                                       min_samples_split=50,
                                       n_jobs=9,
                                       random_state=42
                                      )
        
    if mode == "fast":
        model = RandomForestRegressor(
                                      criterion="squared_error",
                                      n_estimators=100,
                                      bootstrap=True,
                                      max_samples=500,
                                      n_jobs=8,
                                      )
        
    if mode == "slow":        
        model = RandomForestRegressor(
                                      n_estimators=10000,
                                      max_features=3,
                                      min_samples_leaf=10,
                                      min_samples_split=30,
                                      n_jobs=8,
                                     )

        
    if mode == "decision_tree":
        model = DecisionTreeRegressor(
                                      criterion="squared_error",
                                      max_features="sqrt",
                                      min_samples_leaf=50,
                                      min_samples_split=50,
                                     )
    return model


def __run__(fu, params):
    # create genome predictions
    if params["mode"] == "create_genome_predictions":
        models = []
        for cohort in range(params["cohorts"]):
            models.append(pickle.load(open(params["model_dir"]+params["os_sep"]+"model_"+str(cohort+1)+".pickle", "rb")))

        cgp = Create_genome_predictions(models, params)
        cgp.initialize()
        cgp.predict()

    else:
        df = fu.load()


    # evaluation mode
    if params["mode"] == "evaluate_fit":
        models = []

        for cohort in range(params["cohorts"]):
            train_df, test_df          = fu.split(df, cohort+1)
            test_labels, test_features = fu.extract(test_df, "test")

            model = pickle.load(open(params["model_dir"]+params["os_sep"]+"model_"+str(cohort+1)+".pickle", "rb"))
            models.append(model)

        fu.evaluate_models(models)


    # fitting mode
    if params["mode"] == "fit":
        pred_count = 0
        
        # marked (<-) added on 250414 to facilitate feature selection
        if params["permutation_importance"] == True:
            all_importances_mean = np.zeros((0)) # <-

        for cohort in range(params["cohorts"]):
            train_df, test_df            = fu.split(df, cohort+1)
            model                        = __init_model__(params["tree_mode"])
            train_labels, train_features = fu.extract(train_df, "training")
            model.fit(train_features, np.ravel(train_labels))

            train_preds                  = model.predict(train_features)
            fu.evaluate(train_labels, train_preds, cohort, "training")

            test_labels, test_features   = fu.extract(test_df, "test")
            test_preds                   = model.predict(test_features)

            # marked (<-) added on 250414 to facilitate feature selection
            if params["permutation_importance"] == True:
                pis = permutation_importance(model, test_features, test_labels, n_repeats=1, random_state=None)
                if all_importances_mean.shape[0] > 0:
                    all_importances_mean = np.column_stack((all_importances_mean, pis.importances_mean))

                else:
                    all_importances_mean = pis.importances_mean
            # <- ends here

            fu.evaluate(test_labels, test_preds, cohort, "test")
            pred_count += len(test_preds)
            
            if params["tree_mode"] != "decision_tree" and params["print_model"] == True:
                fu.print_model(model, train_features.columns, cohort+1)

        fu.print_evaluation(train_features.columns)
        print("< total predictions:", pred_count)

        # marked (<-) added on 250414 to facilitate feature selection
        if params["permutation_importance"] == True:
            all_importances_std  = np.std(all_importances_mean, axis=1)
            all_importances_mean = np.average(all_importances_mean, axis=1)

            for i in all_importances_mean.argsort()[::-1]:
                #if all_importances_mean[i] - 2 * all_importances_std[i] > 0:
                print(f"{model.feature_names_in_[i]:<8}", f"{all_importances_mean[i]:.3f}", f" +/- {all_importances_std[i]:.3f}")
        
        # <- ends here


    # perform simple prediction
    if params["mode"] == "predict" or params["mode"] == "predict_and_evaluate":
        all_preds = np.zeros((0))
        for cohort in range(params["cohorts"]):
            labels, features = fu.extract(df, "prediction")
            model = pickle.load(open(params["model_dir"]+params["os_sep"]+"model_"+str(cohort+1)+".pickle", "rb"))
            preds = model.predict(features)

            if all_preds.shape[0] > 0: all_preds = np.column_stack((all_preds, preds))
            else:                      all_preds = preds
            if params["mode"] == "predict_and_evaluate": fu.evaluate(labels, preds, cohort, "test")
        
        # modified on 250312 to use averaged predictions for calculation of total r²
        if params["mode"] == "predict_and_evaluate":
            if len(all_preds.shape) > 1: all_preds = np.average(all_preds, axis=1)
            fu.test_preds  = all_preds.tolist()
            fu.test_labels = labels.iloc[:,0].tolist()
            fu.print_evaluation(features.columns)
        
        if params["balance"] == False:
            fu.print_preds(df, all_preds, np.arange(all_preds.shape[0]))


    # conducting predictions based on cohort identity (used to compare prediction power in identical datasets with different labels)
    if params["mode"] == "predict_by_cohorts":
        all_labels = np.zeros((0))
        all_preds  = np.zeros((0))
        all_index  = np.zeros((0))
        cohort_index  = []
        cohort_labels = []
        cohort_preds  = []

        for cohort in range(params["cohorts"]):
            train_df, test_df          = fu.split(df, cohort+1)
            test_labels, test_features = fu.extract(test_df, "test")

            model      = pickle.load(open(params["model_dir"]+params["os_sep"]+"model_"+str(cohort+1)+".pickle", "rb"))
            test_preds = model.predict(test_features)

            cohort_index.append(test_labels.index.tolist())
            cohort_labels.append(test_labels)
            cohort_preds.append(test_preds)

            if all_preds.shape[0] > 0:
                all_labels = np.concatenate((all_labels, test_labels), axis=0)
                all_preds  = np.concatenate((all_preds, test_preds), axis=0)
                all_index  = np.concatenate((all_index, test_labels.index.tolist()), axis=0)

            else:
                all_labels = test_labels
                all_preds  = test_preds
                all_index  = test_labels.index.tolist()

            fu.evaluate(test_labels, test_preds, cohort, "test")

        fu.print_evaluation(test_features.columns)
        fu.print_roc(all_labels, all_preds, cohort_labels, cohort_preds)

        # tentitaviley not working with balanced data
        if params["balance"] == False:       
            fu.print_preds(df, all_preds, all_index)


    # marked (<-) added on 250414 to facilitate hyper-parameter tuning and feature selection
    if params["mode"] == "tuning":
        model  = __init_model__(params["tree_mode"])
        index = []
        for i in range(params["cohorts"]):
            index.append((np.array(df[df["ID:cohort"] != i+1].index), np.array(df[df["ID:cohort"] == i+1].index)))

        param_grid = {'max_features': ['sqrt', 'log2', None],#"max_features": [3, 20],
                      "max_samples": [500, 2000],
                      "min_samples_leaf": [2, 250],
                      "min_samples_split": [2, 250]}

        grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=index, return_train_score=False)
        grid_search.fit(df[[col for col in df.columns if "FEATURE" in col]], df[params["label_id"]])
        print(grid_search.best_params_, grid_search.best_score_, grid_search.best_index_)
    # ends here <-


def main():
    params = {
            "balance":             False,
            "check_features":      ["FEATURE:dist. from last EJC", "FEATURE:lindeboom prediction"], # contains the FEATURES that should checked for missing values, not identical with the SKIP FEATURE
            "cohorts":             5,
            "cohort_values":       {'FEATURE:dist. from last EJC': {0: 649.0836768342951, 1: 637.8458871772308, 2: 651.7632813161315, 3: 651.6981959770337, 4: 626.4562481066013}, 'FEATURE:lindeboom prediction': {0: 0.8053842179778423, 1: 0.8052303500747969, 2: 0.8058681263352281, 3: 0.8056492095653642, 4: 0.8052111893707272}, 'FEATURE:mean expression': {0: 12.346146500210075, 1: 13.072785410672724, 2: 13.05020724168873, 3: 12.745097587544205, 4: 12.626723208675605}, 'FEATURE:median expression': {0: 10.741738350519437, 1: 11.389557189873653, 2: 11.365893584100705, 3: 11.118814414946037, 4: 10.967476155489756}},
            "cohortwise_misses":   False, # if True, missing values are replaced cohort-wise (cohort_values else values)
            "data_dir":            parent_dir+r"\data",
            "decision_threshold":  0.65, # if None, median is used
            # <- specified values are prioritized default values for given features (added on 250428)
            # mean expression: median of mean expressions of full dataset (tcga_expression_info) after filtering for mean_fpkm_unstranded >= 1
            # median expression: median of median expressions of full dataset (tcga_expression_info) after filtering for median_fpkm_unstranded >= 1
            # lindeboom prediction: mean of all lindeboom predictions (full genome build 38), converted mean: 0.6859838726827834, unconverted mean: 0.4419001872943089
            "default_values":      {"FEATURE:mean expression": 4.819906570601459, "FEATURE:median expression": 4.68850178030303, "FEATURE:lindeboom prediction": 0.6859838726827834}, # <-
            "expressions_fname":   "tcga_avg_expressions.txt", # <- added on 250427 to integrate novel model features
            "fname":               "model_training_variants.txt",
            "genome_index":        [0, 10000],
            "hg_build":            "hg38",
            "label_id":            "LABEL:NMD score",
            "label_input":         "ASE", # "ASE", "nonASE" expected input for Lindeboom predictions is nonASE
            "label_output":        "ASE", # "ASE", "nonASE", "nonASE2" # "nonASE2" added on 250410
            "labels_to_class":     False,
            "lindeboom_output":    "ASE", # "ASE", "nonASE", None
            "mode":                "fit", # "create_genome_predictions", "evaluate_fit", "fit", "predict", "predict_and_evaluate", "predict_by_cohorts" "tuning" # 'tuning' added on 250414
            "model_dir":            parent_dir+r"\data\nmdelphi",
            "opt_metric":          "auroc",
            "os_sep":              "\\",
            "permutation_importance": True, # newly added on 250414
            "print_model":         True, # newly added on 250410
            "randomize_balancing": True,
            "remove_gaps":         False, # newly added on 250414
            "session_dir":         None,
            "tag":                 "nmdelphi",
            "target_dir":          parent_dir+r"\data",
            "test_fraction":       1,
            "threads":             1,
            "tree_mode":           "fast", # "fast" "slow" "decision_tree" "fast_class" "slow_class"
            "values":              {'FEATURE:dist. from last EJC': 627.5703261140154, 'FEATURE:lindeboom prediction': 0.8081107710910813, 'FEATURE:mean expression': 16.03474136632121, 'FEATURE:median expression': 14.210174319898858}, # values for kim_reduced_avg
            "verbosity":           1
            }
    
    if params["print_model"] == False:
        print("< warning. model will not be printed.")
    
    # initialize class with container for selected evaluation metrics
    fu = Forest_utils(params, {"auroc": [], "accuracy": [], "r2": [], "rmse": []}, {"auroc": [], "accuracy": [], "r2": [], "rmse": []})
    __run__(fu, params)


if __name__ == '__main__':
    main()