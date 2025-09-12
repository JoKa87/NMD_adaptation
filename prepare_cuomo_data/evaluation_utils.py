import json
import numpy as np
import os
import pandas as pd
from progress.bar import IncrementalBar
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")
from shared_utils import *


class Evaluation_utils():
    def __init__(self, params):
        self.params = params

    
    def _check_annotation(self, appris_annotation):
        if self.params["filter"]["appris_mode"] == "inclusive":
            annotation_found = False
            
            for i in range(len(self.params["filter"]["appris"])):
                if self.params["filter"]["appris"][i] in appris_annotation:
                    annotation_found = True

        if self.params["filter"]["appris_mode"] == "exclusive":
            annotation_found = True
            
            for i in range(len(self.params["filter"]["appris"])):
                if self.params["filter"]["appris"][i] == appris_annotation:
                    annotation_found = False

        return annotation_found


    def _apply_filter(self, blocks, stats, mutants, temp):
        block_mutants = mutants.iloc[temp]

        if len(self.params["filter"]["duplicate_targets"]) > 0:
            for duplicate_target in self.params["filter"]["duplicate_targets"]:
                block_mutants = block_mutants.drop_duplicates(subset=duplicate_target)
        
        if len(self.params["filter"]["appris_priorities"]) > 0:
            appris_priorities = {key: [] for key in self.params["filter"]["appris_priorities"]}
            [appris_priorities[block_mutants.iloc[i].loc["ID:appris annotation"]].append(i) for i in range(block_mutants.shape[0])
             if block_mutants.iloc[i].loc["ID:appris annotation"] in appris_priorities]

            temp = []
            i = 0
            while i < len(list(appris_priorities.keys())) and len(temp) == 0:
                if len(appris_priorities[list(appris_priorities.keys())[i]]) > 0:
                    temp                                      = appris_priorities[list(appris_priorities.keys())[i]]
                    stats[list(appris_priorities.keys())[i]] += len(temp)

                i += 1

            # edited 250320 to allow selection of only a single appris term even if multiple terms are present with identical appris priority
            if self.params["filter"]["appris_redundant"] == True:  block_mutants = block_mutants.iloc[temp] # old version without if statement
            if self.params["filter"]["appris_redundant"] == False: block_mutants = block_mutants.iloc[[temp[0]]] # newly added

        for col in block_mutants.columns:
            blocks[col].extend(block_mutants[col].tolist())

        return blocks, stats


    def apply_filter(self, mutants):
        # filter using read-cutoff
        mutants = mutants[mutants["LABEL:total counts"] >= self.params["filter"]["min_reads"]]

        # filter for entries present in the hg38 database as well
        '''
        filtered transcripts checked on ensembl.org on 250319
        deprecated without novel assignment:
        ENST00000278360, ENST00000378880, ENST00000532848, ENST00000431837, ENST00000379218, ENST00000444838, ENST00000393992, ENST00000534850, ENST00000331336, ENST00000590902,
        ENST00000328550, ENST00000291187, ENST00000591604, ENST00000420102, ENST00000339613, ENST00000359009, ENST00000301439, ENST00000544146, ENST00000426889, ENST00000438534,
        ENST00000258324, ENST00000373023, ENST00000422516, ENST00000215780, ENST00000398141, ENST00000542275, ENST00000260630, ENST00000454352, ENST00000264122, ENST00000393363,
        ENST00000398965, ENST00000438599, ENST00000543135, ENST00000541730, ENST00000392738, ENST00000525043, ENST00000434271, ENST00000311322, ENST00000427548

        deprecated with novel assignment:
        ENST00000357647 (several novel assignments), 
        ENST00000451963 (novel assignment with missing 5'CDS)
        ENST00000540964 (1 novel assignment)
        '''
        
        if self.params["filter"]["hg38"] == True:
            # load knownGene database
            init_shape                  = mutants.shape[0]
            init_variants               = mutants.drop_duplicates(subset="ID:transcript id").shape[0]
            known_gene                  = pd.read_csv(self.params["data_dir1"]+self.params["os_sep"]+"hg38_knownGene_appended.txt", delimiter=",", index_col=False)
            known_gene["transcript id"] = [known_gene.iloc[i].loc["transcript id"].split(".")[0] for i in range(known_gene.shape[0])]
            # marked (<-) added on 250428
            mutants["ID:transcript id"] = [mutants.iloc[i].loc["ID:transcript id"].split(".")[0] for i in range(mutants.shape[0])] # <- added
            test                        = mutants[~mutants["ID:transcript id"].isin(known_gene["transcript id"])]
            mutants                     = mutants[mutants["ID:transcript id"].isin(known_gene["transcript id"])]
            print("< hg38 filtering reduced data from", init_shape, "to", mutants.shape[0], "and from", init_variants,
                  "to", mutants.drop_duplicates(subset="ID:transcript id").shape[0], "variants")

            transcripts = test.drop_duplicates(subset="ID:transcript id")["ID:transcript id"].tolist()
            print(transcripts)

        # filter to remove or keep specified elements of appris annotation
        if self.params["filter"]["appris"] != None:
            selected_index = [self._check_annotation(mutants.iloc[i].loc["ID:appris annotation"]) for i in range(mutants.shape[0])]
            mutants        = mutants[selected_index]

        # filter according to specified appris priorities
        if self.params["filter"]["block_targets"] != None:
            mutants  = mutants.sort_values(by=self.params["filter"]["block_targets"])
            last_ids = {block_target: None for block_target in self.params["filter"]["block_targets"]}
            temp     = []
            
            blocks   = {col: [] for col in mutants.columns}
            stats    = {key: 0 for key in self.params["filter"]["appris_priorities"]}
            bar      = IncrementalBar("filtering blocks", max=mutants.shape[0])

            for i in range(mutants.shape[0]):
                check = [block_target for block_target in last_ids if last_ids[block_target] == mutants.iloc[i].loc[block_target]]

                # append existing block if elements are identical
                if len(temp) == 0 or len(check) == len(self.params["filter"]["block_targets"]):
                    temp.append(i)

                # create new block if elements diverge
                elif len(temp) > 0:
                    blocks, stats = self._apply_filter(blocks, stats, mutants, temp)
                    temp          = [i]

                last_ids = {block_target: mutants.iloc[i].loc[block_target] for block_target in self.params["filter"]["block_targets"]}

                bar.next()
            bar.finish()

            # create last block
            if len(temp) > 0: blocks, stats = self._apply_filter(blocks, stats, mutants, temp)

        mutants         = pd.DataFrame(blocks)
        print("< appris stats")
        print(json.dumps(stats, indent=4))
        appris_sum      = np.sum(stats[key] for key in stats)
        if self.params["filter"]["appris_redundant"] == True and appris_sum != mutants.shape[0]: print("< inconsistent sizes @apply_filter", appris_sum, "/", mutants.shape[0])
        duplicates_test = mutants.drop_duplicates(subset=self.params["filter"]["block_targets"])
        print("< total blocks:", mutants.shape[0], " duplicates test after filtering: ", mutants.shape[0]-duplicates_test.shape[0])
        return mutants
    

    def separate_counts(self, mutants):
        mutant_dict = {col: [] for col in mutants.columns}
        bar = IncrementalBar(set_bar("separating counts"), max=mutants.shape[0])

        for i in range(mutants.shape[0]):
            alt_counts            = mutants.iloc[i].loc["LABEL:alt counts"].replace("[", "").replace(" ", "").strip("]").split(",")
            total_counts          = mutants.iloc[i].loc["LABEL:total counts"].replace("[", "").replace(" ", "").strip("]").split(",")
            filtered_alt_counts   = [float(alt_counts[j]) for j in range(len(alt_counts)) if alt_counts[j] != "nan" and float(total_counts[j]) >= self.params["min_reads"]]
            filtered_total_counts = [float(total_counts[j]) for j in range(len(total_counts)) if total_counts[j] != "nan" and float(total_counts[j]) >= self.params["min_reads"]]

            sample_ids            = mutants.iloc[i].loc["ID:sample id"].replace("[", "").replace(" ", "").replace("'", "").strip("]").split(",")
            sample_ids            = [sample_ids[j] for j in range(len(sample_ids)) if float(total_counts[j]) >= self.params["min_reads"]]

            if len(filtered_alt_counts) != len(filtered_total_counts) or len(filtered_alt_counts) != len(sample_ids):
                print("< length mismatch @separate_counts")

            for j in range(len(filtered_alt_counts)):
                for col in mutants.columns:
                    if "LABEL" not in col and "ID:sample id" not in col:
                        mutant_dict[col].append(mutants.iloc[i].loc[col])

                    elif "ID:sample id" in col:
                        mutant_dict[col].append(sample_ids[j])
                    
                    elif "alt counts" in col:
                        mutant_dict[col].append(filtered_alt_counts[j])

                    elif "total counts" in col:
                        mutant_dict[col].append(filtered_total_counts[j])

            bar.next()
        bar.finish()

        mutants = pd.DataFrame(mutant_dict)
        mutants.insert(mutants.columns.get_loc("LABEL:total counts")+1, "LABEL:NMD score",
                       [(mutants.iloc[i].loc["LABEL:total counts"]-mutants.iloc[i].loc["LABEL:alt counts"]) / mutants.iloc[i].loc["LABEL:total counts"]
                        for i in range(mutants.shape[0])])
        return mutants
    
    
    def sum_counts(self, mutants):
        summed_alt_counts   = []
        summed_total_counts = []
        bar = IncrementalBar(set_bar("summing counts"), max=mutants.shape[0])


        for i in range(mutants.shape[0]):
            alt_counts   = mutants.iloc[i].loc["LABEL:alt counts"].replace("[", "").replace(" ", "").strip("]").split(",")
            total_counts = mutants.iloc[i].loc["LABEL:total counts"].replace("[", "").replace(" ", "").strip("]").split(",")
            alt_counts_sum   = 0
            total_counts_sum = 0

            for j in range(len(alt_counts)):
                alt_counts[j]   = alt_counts[j].replace("np.float64", "").replace("(", "").replace(")", "") # demanded from output created by run on MacOS
                total_counts[j] = total_counts[j].replace("np.float64", "").replace("(", "").replace(")", "") # demanded from output created by run on MacOS

                if alt_counts[j] != "nan": 
                    alt_counts_sum   += float(alt_counts[j])
                    total_counts_sum += float(total_counts[j])

            summed_alt_counts.append(alt_counts_sum)
            summed_total_counts.append(total_counts_sum)
            
            bar.next()
        bar.finish()

        mutants["LABEL:alt counts"]   = summed_alt_counts
        mutants["LABEL:total counts"] = summed_total_counts
        mutants.insert(mutants.columns.get_loc("LABEL:total counts")+1, "LABEL:NMD score",
                       [(summed_total_counts[i]-summed_alt_counts[i]) / summed_total_counts[i] for i in range(len(summed_total_counts)) if summed_total_counts[i] > 0])
        
        return mutants