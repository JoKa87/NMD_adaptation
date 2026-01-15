import numpy as np
import scipy

from analyze_predictions_utils import *


class Read_predictions_utils(Analyze_predictions_utils):
    def __init__(self, params):
        super().__init__(params, {})

        self.params                       = params
        self.params["appris_selection"]   = False
        self.params["use_variant_filter"] = False

            
    def compare_scores(self, scores, targets, predictions):
        matches    = 0
        mismatches = 0
        report     = {}

        for i in range(len(scores)):
            if pd.isna(scores[i]) == False and scores[i] != round(targets.iloc[i].loc[self.params["target"]["prediction_identifier"]], 4):
                mismatches += 1
                transcript  = targets.iloc[i].loc[self.params["target"]["transcript_identifier"]]
                print("< score mismatch @", transcript, scores[i], "/", targets.iloc[i].loc[self.params["target"]["prediction_identifier"]])

                report[transcript] = targets.iloc[i].loc[self.params["target"]["position_identifier"]]

                selected_predictions = predictions[predictions[self.params["prediction_transcript_identifier"]] == transcript]
                for j in range(selected_predictions.shape[0]):
                    print("  transcript", j, ", cds position:", targets.iloc[i].loc[self.params["target"]["position_identifier"]],
                          ", cds length:", len(selected_predictions.iloc[j].loc["predictions"]))
                
            elif pd.isna(scores[i]) == False:
                matches += 1

        print("<", mismatches, "mismatches and", matches, "matches",
            ". pearsonr", scipy.stats.pearsonr([score for score in scores if pd.isna(score) == False], [targets.iloc[i].loc[self.params["target"]["prediction_identifier"]] for i in range(targets.shape[0]) if scores[i] != None]))

        print(report)
        print(len(report))


    def print_scores(self, scores, targets):
        targets["FEATURE:prediction"] = scores
        targets.to_csv(path_or_buf=self.params["target_dir"]+self.params["os_sep"]+self.params["target_fname"].split(".")[0]+"_predictions.txt", sep=",", index=False)


    def read_scores(self, targets, predictions):
        scores             = []
        failed_identifiers = []

        bar = IncrementalBar(set_bar("reading scores"), max=targets.shape[0])
        for i in range(targets.shape[0]):
            selected_predictions = predictions[predictions[self.params["prediction_transcript_identifier"]] == targets.iloc[i].loc[self.params["target"]["transcript_identifier"]]]

            if selected_predictions.shape[0] > 0:
                selected_scores = [selected_predictions.iloc[j].loc["predictions"][targets.iloc[i].loc[self.params["target"]["position_identifier"]]]
                                   for j in range(selected_predictions.shape[0])
                                   if targets.iloc[i].loc[self.params["target"]["position_identifier"]] < len(selected_predictions.iloc[j].loc["predictions"])]
                
                if len(selected_scores) > 0:
                    scores.append(np.mean(selected_scores))

                if len(selected_scores) == 0:
                    scores.append(None)
                    print("< requested position exceeded for", targets.iloc[i].loc[self.params["target"]["transcript_identifier"]], "position:", targets.iloc[i].loc[self.params["target"]["position_identifier"]],
                          "compared sizes:", [len(selected_predictions.iloc[j].loc["predictions"]) for j in range(selected_predictions.shape[0])])
                 
            else:
                scores.append(None)
                if targets.iloc[i].loc[self.params["target"]["transcript_identifier"]] not in failed_identifiers:
                    failed_identifiers.append(targets.iloc[i].loc[self.params["target"]["transcript_identifier"]])

            bar.next()
        bar.finish()
        
        print("<", len(failed_identifiers), "identifiers could not be found.")
        print(failed_identifiers)
        return scores