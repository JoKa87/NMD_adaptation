# marked (<-) added on 250427 to allow integration of allele frequency
# from gnomad_db.database import gnomAD_DB # <- added
import math
import os
import pandas as pd
from progress.bar import IncrementalBar
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\apply_lindeboom")
sys.path.insert(0, parent_dir+"\\shared")

import threading
import time

from apply_lindeboom_utils import *
from shared_utils import *


class Create_genome_predictions():
    def __init__(self, models, params):
        self.expressions           = None
        self.gene_identifier       = None
        self.genome                = None
        #self.gnomad_db             = None
        self.knowngene             = None
        self.lindeboom_dictionary  = None
        self.models                = models
        self.params                = params
        self.pattern_cutoff        = 5
        # marked (<-) added on 250428 (collection of extracted features)
        self.selected_features = ["FEATURE:5'utr-size", "FEATURE:3'utr-size", "FEATURE:upstream cds exons", "FEATURE:downstream cds exons", "FEATURE:downstream exons",
                                  "FEATURE:downstream EJC density", "FEATURE:cds EJC density", "FEATURE:downstream cds EJC density", "FEATURE:ptc downstream distance",
                                  "FEATURE:dist. from last EJC", "FEATURE:ptc cds position", "FEATURE:ptc exon position", "FEATURE:ptc-wt stop codon distance",
                                  "FEATURE:last cds exon", "FEATURE:50 nt to last cds EJC", "FEATURE:ptc exon size", "FEATURE:total GC count exons",
                                  "FEATURE:downstream GC ptc count ptc exon", "FEATURE:lindeboom prediction", "FEATURE:mean expression", "FEATURE:median expression"] # <- added

        # marked (<-) added on 250428 to allow flexible usage of params->default_values or params->values (see initialize)
        self.values                = {} # <- added


    # marked (<-) added/replaced on 250427 to adapt to novel model with new features
    def extract_features(self, features, fails, it):
        #transient_features   = {**{feature: 0 for feature in ["3'utr-size", "current exon size", "last ejc", "regular lindeboom", "total cds size", "total exon size"]},
        #                        **{"lindeboom prediction": [], "ptc cds position": []}} # replaced
        transient_features   = {**{feature: 0 for feature in ["5'utr-size", "3'utr-size", "cds exon", "cds exons", "current exon size", "last ejc",
                                                              "regular lindeboom", "total cds size", "total exon size", "total lindeboom size", "last cds correction"]},
                                **{"downstream cds EJC density": [], "lindeboom prediction": [], "ptc cds position": []}, # , "ptc abs. position": []
                                **{"all exons": "", "all cds": ""}} # <- added
        init_size            = len(features["FEATURE:dist. from last EJC"])

        transcript           = self.knowngene.iloc[it]
        #transcript_id        = transcript.loc["transcript id"].split(".")[0]
        transcript_id        = transcript.loc[self.gene_identifier].split(".")[0]
        fails[transcript_id] = {"cds_mismatch": 0, "chromosome_not_found": 0, "last_ejc_unknown": 0, "lindeboom_cds_size_error": 0, "lindeboom_exon_index_error": 0, "lindeboom_exon_size_error": 0,
                                "lindeboom_not_found": 0, "lindeboom_size_test_error": 0, "total_exon_size_unknown": 0}

        exonsstart = [int(i) for i in transcript.loc["exonsstart"].split(",") if len(i) > 0] # condition required because last letter is a comma
        exonsend   = [int(i) for i in transcript.loc["exonsend"].split(",") if len(i) > 0]   # condition required because last letter is a comma

        # check whether chromosome is present
        chr = transcript.loc["chr"]
        if chr not in self.genome: fails[transcript_id]["chromosome_not_found"] += 1
        
        total_exon_size_counter = [exonsend[i]-exonsstart[i] for i in range(transcript.loc["exons"])]
        if len(total_exon_size_counter) > 0: transient_features["total exon size"] = np.sum(total_exon_size_counter)
        else:                                fails[transcript_id]["total_exon_size_unknown"] += 1

        if transcript.loc["strand"] == "+":
            last_ejc_counter = [exonsend[i]-exonsstart[i] for i in range(transcript.loc["exons"]-1)]
            if len(last_ejc_counter) > 0:    transient_features["last ejc"] = np.sum(last_ejc_counter)
            else:                            fails[transcript_id]["last_ejc_unknown"] += 1
            i = 0

        if transcript.loc["strand"] == "-":
            last_ejc_counter = [exonsend[transcript.loc["exons"]-i-1]-exonsstart[transcript.loc["exons"]-i-1] for i in range(transcript.loc["exons"]-1)]
            if len(last_ejc_counter) > 0:    transient_features["last ejc"] = np.sum(last_ejc_counter)
            else:                            fails[transcript_id]["last_ejc_unknown"] += 1
            i = transcript.loc["exons"]-1

        # determine number of exons (<- added)
        transient_features["cds exons"] = len([i for i in range(len(exonsstart))
                                               if not (exonsend[i] < transcript.loc["cdsstart"] or exonsstart[i] > transcript.loc["cdsend"])]) # <- added

        # determine lindeboom size for later testing
        if transcript.loc["uc id"].split(".")[0] in self.lindeboom_dictionary.keys(): # <- added
	        transient_features["total lindeboom size"] = np.sum([len(self.lindeboom_dictionary[transcript.loc["uc id"].split(".")[0]][i]["scores"])
                                                                 for i in range(len(self.lindeboom_dictionary[transcript.loc["uc id"].split(".")[0]]))]) # <- added
        
        cds_index = None
        cds_start = False
        show = False
        while ((transcript.loc["strand"] == "+" and i < transcript.loc["exons"]) or (transcript.loc["strand"] == "-" and i >= 0)):
            if show == True: print("i", i, "exons", exonsstart[i], "/", exonsend[i], "cds", transcript.loc["cdsstart"], "/", transcript.loc["cdsend"])

            if transcript.loc["strand"] == "+":   transient_features["all exons"] += self.genome[chr][exonsstart[i]:exonsend[i]] # <- added
            elif transcript.loc["strand"] == "-": transient_features["all exons"] += self.invert(self.genome[chr][exonsstart[i]:exonsend[i]]) # <- added

            # calculations for exons with no CDS (UTR-only)
            if exonsend[i] < transcript.loc["cdsstart"] or exonsstart[i] > transcript.loc["cdsend"]:
                if cds_start == True:
                    if transcript.loc["strand"] == "+":   transient_features["3'utr-size"] += exonsend[i]-exonsstart[i]
                    elif transcript.loc["strand"] == "-": transient_features["3'utr-size"] += exonsend[i]-exonsstart[i]
                    if show == True: print("3'UTR", transient_features["3'utr-size"])

                else: # <- added
                    transient_features["5'utr-size"] += exonsend[i]-exonsstart[i] # <- added
                    if show == True: print("5'UTR", transient_features["5'utr-size"]) # <- added

            # calculations for exons with CDS
            else:
                cds_start = True
                if cds_index == None:
                    if transcript.loc["strand"] == "+":
                        cds_index = 0

                    if transcript.loc["strand"] == "-":
                        if transcript.loc["uc id"].split(".")[0] in self.lindeboom_dictionary.keys():
                            cds_index = len(self.lindeboom_dictionary[transcript.loc["uc id"].split(".")[0]])-1

                # check if exon contains non-cds
                if exonsstart[i] < transcript.loc["cdsstart"]: exon_cds_start = transcript.loc["cdsstart"]
                else:                                          exon_cds_start = exonsstart[i]
                if exonsend[i] > transcript.loc["cdsend"]:     exon_cds_end   = transcript.loc["cdsend"]
                else:                                          exon_cds_end   = exonsend[i]

                if transcript.loc["strand"] == "+":   transient_features["all cds"] += self.genome[chr][exon_cds_start:exon_cds_end] # <- added
                elif transcript.loc["strand"] == "-": transient_features["all cds"] += self.invert(self.genome[chr][exon_cds_start:exon_cds_end]) # <- added
                if show == True: print("CDS", exon_cds_start, exon_cds_end, exon_cds_end-exon_cds_start)

                is_last_cds_exon = False
                if transcript.loc["strand"] == "+":
                    if transcript.loc["cdsend"] <= exonsend[i]:     is_last_cds_exon = True

                if transcript.loc["strand"] == "-":
                    if transcript.loc["cdsstart"] >= exonsstart[i]: is_last_cds_exon = True

                # introduced in order to take care of a few genes in which stop codon is distributed over two exons
                if is_last_cds_exon == True and exon_cds_end-exon_cds_start < 3: # <- added
                    transient_features["last cds correction"] = 3-(exon_cds_end-exon_cds_start) # <- added

                
                if show == True: print("cds", exon_cds_start, "/", exon_cds_end)

                if transcript.loc["strand"] == "+": # <- added
                    transient_features["5'utr-size"] += exon_cds_start-exonsstart[i] # <- added

                elif transcript.loc["strand"] == "-": # <- added
                    transient_features["5'utr-size"] += exonsend[i]-exon_cds_end # <- added

                if transcript.loc["strand"] == "+":
                    transient_features["3'utr-size"] += exonsend[i]-exon_cds_end

                elif transcript.loc["strand"] == "-":
                    transient_features["3'utr-size"] += exon_cds_start-exonsstart[i]                
                    
                # iterate over exon
                if transcript.loc["strand"] == "+":
                    j = 0
                    if is_last_cds_exon == False: max_it = exon_cds_end-exon_cds_start
                    if is_last_cds_exon == True:  max_it = exon_cds_end-exon_cds_start-3

                if transcript.loc["strand"] == "-":
                    j = exon_cds_end-exon_cds_start-1
                    if is_last_cds_exon == False: max_it = 0
                    if is_last_cds_exon == True:  max_it = 3
                
                while ((transcript.loc["strand"] == "+" and j < max_it) or (transcript.loc["strand"] == "-" and j >= max_it)):                    
                    if transcript.loc["strand"] == "+":
                        # <- added from here
                        features["FEATURE:ptc cds position"].append(len(transient_features["downstream cds EJC density"]))
                        features["FEATURE:ptc exon size"].append(exonsend[i]-exonsstart[i]-1) # '-1' is one less than actual exon size (as in prepare_data), ptc side excluded!
                        features["FEATURE:upstream cds exons"].append(transient_features["cds exon"])
                        features["FEATURE:downstream cds exons"].append(transient_features["cds exons"]-transient_features["cds exon"])
                        transient_features["downstream cds EJC density"].append(transient_features["cds exons"]-transient_features["cds exon"])
                        #transient_features["ptc abs. position"].append(exon_cds_start+j)

                        if transient_features["cds exons"]-transient_features["cds exon"] != 1:
                            features["FEATURE:last cds exon"].append(0)

                            if transient_features["cds exons"]-transient_features["cds exon"] == 2 and exonsend[i]-1-exon_cds_start-j <= 50: 
                                features["FEATURE:50 nt to last cds EJC"].append(1)

                            else:
                                features["FEATURE:50 nt to last cds EJC"].append(0)
                        
                        else:
                            features["FEATURE:last cds exon"].append(1)
                            features["FEATURE:50 nt to last cds EJC"].append(1)
                        # <- until here
                            
                        transient_features["ptc cds position"].append(transient_features["total cds size"]+j)
                        features["FEATURE:downstream exons"].append(transcript.loc["exons"]-i)
                        features["FEATURE:downstream EJC density"].append((transcript.loc["exons"]-i-1) / transient_features["total exon size"])
                        #features["FEATURE:ptc downstream distance"].append(exonsend[i]-exon_cds_start-j)
                        features["FEATURE:ptc downstream distance"].append(exonsend[i]-1-exon_cds_start-j)
                        features["FEATURE:ptc exon position"].append(transient_features["current exon size"]+exon_cds_start-exonsstart[i]+j)

                        if transient_features["last ejc"] > 0: features["FEATURE:dist. from last EJC"].append(transient_features["last ejc"]-features["FEATURE:ptc exon position"][-1]-1)
                        # marked (<-) added / removed on 250428 to instead use upstream distance as measure
                        # else:                                  features["FEATURE:dist. from last EJC"].append(self.params["values"]["FEATURE:dist. from last EJC"]) # <- removed
                        else:                                  features["FEATURE:dist. from last EJC"].append(-((exonsend[i]-exonsstart[i]-1)-(exonsend[i]-1-exon_cds_start-j))) # <- added
                        
                        gc_counter = [1 for k in range(exon_cds_start+j+1, exonsend[i], 1) if self.genome[chr][k] in ["G", "C"]]
                        #if exonsend[i]-exon_cds_start-j >= self.pattern_cutoff: features["FEATURE:downstream GC ptc count ptc exon"].append(np.sum(gc_counter)/(exonsend[i]-exon_cds_start-j))
                        if exonsend[i]-exon_cds_start-j-1 >= self.pattern_cutoff: features["FEATURE:downstream GC ptc count ptc exon"].append(np.sum(gc_counter)/(exonsend[i]-exon_cds_start-j-1))
                        else:                                                     features["FEATURE:downstream GC ptc count ptc exon"].append(0.5)

                        features["exon_cds_start"].append(exon_cds_start)

                    if transcript.loc["strand"] == "-":
                        # <- added from here
                        features["FEATURE:ptc cds position"].append(len(transient_features["downstream cds EJC density"]))
                        features["FEATURE:ptc exon size"].append(exonsend[i]-exonsstart[i]-1) # '-1' is one less than actual exon size (as in prepare_data), ptc side excluded!
                        features["FEATURE:upstream cds exons"].append(transient_features["cds exon"])
                        features["FEATURE:downstream cds exons"].append(transient_features["cds exons"]-transient_features["cds exon"])
                        transient_features["downstream cds EJC density"].append(transient_features["cds exons"]-transient_features["cds exon"])
                        #transient_features["ptc abs. position"].append(exon_cds_start+j)

                        if transient_features["cds exons"]-transient_features["cds exon"] != 1:
                            features["FEATURE:last cds exon"].append(0)

                            if transient_features["cds exons"]-transient_features["cds exon"] == 2 and exon_cds_start+j-exonsstart[i] <= 50: 
                                features["FEATURE:50 nt to last cds EJC"].append(1)

                            else:
                                features["FEATURE:50 nt to last cds EJC"].append(0)
                        
                        else:
                            features["FEATURE:last cds exon"].append(1)
                            features["FEATURE:50 nt to last cds EJC"].append(1)
                        # <- until here

                        transient_features["ptc cds position"].append(transient_features["total cds size"]+exon_cds_end-exon_cds_start-1-j)
                        features["FEATURE:downstream exons"].append(i+1)
                        features["FEATURE:downstream EJC density"].append(i / transient_features["total exon size"])
                        features["FEATURE:ptc downstream distance"].append(exon_cds_start+j-exonsstart[i])
                        features["FEATURE:ptc exon position"].append(transient_features["current exon size"]+exonsend[i]-exon_cds_end+exon_cds_end-exon_cds_start-1-j)
                        
                        if transient_features["last ejc"] > 0: features["FEATURE:dist. from last EJC"].append(transient_features["last ejc"]-features["FEATURE:ptc exon position"][-1]-1)
                        # marked (<-) added / removed on 250428 to instead use upstream distance as measure
                        # else:                                  features["FEATURE:dist. from last EJC"].append(self.params["values"]["FEATURE:dist. from last EJC"]) # <- removed
                        else:                                  features["FEATURE:dist. from last EJC"].append(-((exonsend[i]-exonsstart[i]-1)-(exon_cds_start+j-exonsstart[i]))) # <- added

                        gc_counter = [1 for k in range(exonsstart[i], exon_cds_start+j, 1) if self.genome[chr][k] in ["G", "C"]]

                        if exon_cds_start+j-exonsstart[i] >= self.pattern_cutoff: features["FEATURE:downstream GC ptc count ptc exon"].append(np.sum(gc_counter)/(exon_cds_start+j-exonsstart[i]))
                        else:                                                     features["FEATURE:downstream GC ptc count ptc exon"].append(0.5)

                        if is_last_cds_exon == False: features["exon_cds_start"].append(exon_cds_start)
                        if is_last_cds_exon == True:  features["exon_cds_start"].append(exon_cds_start+3)

                    if transcript.loc["strand"] == "+": lindeboom_prediction, fails, failed = self.get_lindeboom_predictions(fails, transcript_id, transcript.loc["uc id"], cds_index, j,
                                                                                                                             exon_cds_end-exon_cds_start, is_last_cds_exon)
                    if transcript.loc["strand"] == "-": lindeboom_prediction, fails, failed = self.get_lindeboom_predictions(fails, transcript_id, transcript.loc["uc id"], cds_index, j-max_it,
                                                                                                                             exon_cds_end-exon_cds_start, is_last_cds_exon)
                    
                    # count the no. of regularly read Lindeboom predictions for later filtering
                    if failed == False: transient_features["regular lindeboom"] += 1
                    
                    #features["FEATURE:lindeboom prediction"].append(lindeboom_prediction)
                    transient_features["lindeboom prediction"].append(lindeboom_prediction)

                    if transcript.loc["strand"] == "+": j += 1
                    if transcript.loc["strand"] == "-": j -= 1

                    if show == True and j < 0:
                        print("j", j, "cds size", exon_cds_end-exon_cds_start, "dist. from last EJC", features["FEATURE:dist. from last EJC"][-1], "downstream exons", features["FEATURE:downstream exons"][-1],
                              "downstream EJC density", round(features["FEATURE:downstream EJC density"][-1], 4), "ptc downstream distance", features["FEATURE:ptc downstream distance"][-1],
                              "ptc exon position", features["FEATURE:ptc exon position"][-1])
                
                transient_features["cds exon"]       += 1 # <- added
                transient_features["total cds size"] += exon_cds_end-exon_cds_start
                if transcript.loc["strand"] == "+" and cds_index != None: cds_index += 1
                if transcript.loc["strand"] == "-" and cds_index != None: cds_index -= 1

            transient_features["current exon size"] += exonsend[i]-exonsstart[i]
            if transcript.loc["strand"] == "+": i += 1
            if transcript.loc["strand"] == "-": i -= 1

        #features["transcript id"].extend([transcript_id for _ in range(len(features["FEATURE:ptc downstream distance"])-init_size)])
        features[self.gene_identifier].extend([transcript_id for _ in range(len(features["FEATURE:ptc downstream distance"])-init_size)])
        features["FEATURE:5'utr-size"].extend([transient_features["5'utr-size"] for _ in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <- added
        features["FEATURE:3'utr-size"].extend([transient_features["3'utr-size"] for _ in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) 

        if transcript.loc["strand"] == "+": # +1 is empirical to match prepare data output
            #features["FEATURE:ptc-wt stop codon distance"].extend([transient_features["total cds size"]-transient_features["ptc cds position"][i]+1 for i in range(len(transient_features["ptc cds position"]))])
            features["FEATURE:ptc-wt stop codon distance"].extend([transient_features["total cds size"]-transient_features["ptc cds position"][i] for i in range(len(transient_features["ptc cds position"]))])

        if transcript.loc["strand"] == "-":
            features["FEATURE:ptc-wt stop codon distance"].extend([transient_features["total cds size"]-transient_features["ptc cds position"][i] for i in range(len(transient_features["ptc cds position"]))])

        # determine cds EJC density (<- added)
        features["FEATURE:cds EJC density"].extend([(transient_features["cds exons"]-1) / transient_features["total cds size"]
                                                    for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-

        # determine downstream cds EJC density (<- added)
        features["FEATURE:downstream cds EJC density"].extend([(transient_features["downstream cds EJC density"][i]-1) / transient_features["total cds size"]
                                                               for i in range(len(transient_features["downstream cds EJC density"]))]) # <-

        # determine all exon GC count (<- added)
        gc_ratio = len([1 for base in transient_features["all exons"] if base in ["G", "C"]])/len(transient_features["all exons"]) # <-
        features["FEATURE:total GC count exons"].extend([gc_ratio for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-

        # append expression features (<- added)
        selected_expressions = pd.DataFrame() # <-

        if type(self.knowngene.iloc[it].loc["gene id"]) == str: # <-
            selected_expressions = self.expressions[self.expressions["gene_id"] == self.knowngene.iloc[it].loc["gene id"].split(".")[0]] # <-

        if selected_expressions.shape[0] > 0: # <-
            features["FEATURE:mean expression"].extend([selected_expressions.iloc[0].loc["mean_fpkm_unstranded"]
                                                        if pd.isna(selected_expressions.iloc[0].loc["mean_fpkm_unstranded"]) == False
                                                        else self.values["FEATURE:mean expression"]
                                                        for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-
            features["FEATURE:median expression"].extend([selected_expressions.iloc[0].loc["median_fpkm_unstranded"]
                                                          if pd.isna(selected_expressions.iloc[0].loc["median_fpkm_unstranded"]) == False
                                                          else self.values["FEATURE:median expression"]
                                                          for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-

        else: # <-
            features["FEATURE:mean expression"].extend([self.values["FEATURE:mean expression"]
                                                        for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-
            features["FEATURE:median expression"].extend([self.values["FEATURE:median expression"]
                                                          for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-

        # append allele frequencies (currently, only estimates as in Lindeboom et al. 2019) <-
        # features["FEATURE:AF"].extend([self.values["FEATURE:AF"] for i in range(len(features["FEATURE:ptc downstream distance"])-init_size)]) # <-

        # check if total cds size equals lindeboom size, if not replace all selected Lindeboom predictions by default (to synchronize behaviour from apply_lindeboom)
        #if transient_features["total cds size"]-3 == transient_features["regular lindeboom"]: # <- removed
        if transient_features["total cds size"]-3 == transient_features["total lindeboom size"]: # <- added
            features["FEATURE:lindeboom prediction"].extend(transient_features["lindeboom prediction"])

        else:
            # marked (<-) added / removed on 250428 to allow usage of default values and training values
            #features["FEATURE:lindeboom prediction"].extend([self.params["values"]["FEATURE:lindeboom prediction"]
            #                                                for _ in range(len(transient_features["lindeboom prediction"]))]) # <- removed
            features["FEATURE:lindeboom prediction"].extend([self.values["FEATURE:lindeboom prediction"]
                                                             for _ in range(len(transient_features["lindeboom prediction"]))]) # <- added
            fails[transcript_id]["cds_mismatch"] += 1

        # correct for small error due to stop codon distributed over two exons
        if transient_features["last cds correction"] > 0: # <- added
            for col in features: # <- added
                for _ in range(transient_features["last cds correction"]): # <- added
                    features[col].pop() # <- added


        if len(features["FEATURE:lindeboom prediction"])-init_size != len(transient_features["all cds"])-3: # <- added
            print("< inconsistent size for "+transcript.loc["transcript id"] # <- added
                  +" ("+str(len(features["FEATURE:lindeboom prediction"])-init_size)+"/"+str(len(transient_features["all cds"])-3)+")") # <- added
            exit() # <- added

        return features, fails


    def get_lindeboom_predictions(self, fails, transcript_id, ucid, exon_index, cds_position, cds_size, is_last_cds_exon):
        failed = False
        # marked (<-) added / removed on 250428 to allow usage of default values and training values
        # lindeboom_prediction = self.params["values"]["FEATURE:lindeboom prediction"] # <- removed
        lindeboom_prediction = self.values["FEATURE:lindeboom prediction"] # <- added

        if is_last_cds_exon == True: cds_size -= 3

        if ucid.split(".")[0] in self.lindeboom_dictionary.keys():
            if exon_index != None and exon_index >= 0 and exon_index < len(self.lindeboom_dictionary[ucid.split(".")[0]]):
                # marked (<-) added/removed on 250430 to avoid error due to start codon split over two exons
                #if cds_size != len(self.lindeboom_dictionary[ucid.split(".")[0]][exon_index]["scores"]): # <- removed
                #    fails[transcript_id]["lindeboom_size_test_error"] += 1; failed = True # <- removed
                
                #elif cds_position < len(self.lindeboom_dictionary[ucid.split(".")[0]][exon_index]["scores"]): # <- removed
                if cds_position < len(self.lindeboom_dictionary[ucid.split(".")[0]][exon_index]["scores"]): # <- added
                    if self.params["lindeboom_output"] == "nonASE":
                        lindeboom_prediction = self.lindeboom_dictionary[ucid.split(".")[0]][exon_index]["scores"][cds_position]
                    
                    if self.params["lindeboom_output"] == "ASE":
                        lindeboom_prediction = 1/(2*math.exp(-math.log(2)*self.lindeboom_dictionary[ucid.split(".")[0]][exon_index]["scores"][cds_position]))
                
                else:
                    fails[transcript_id]["lindeboom_cds_size_error"] += 1; failed = True; #print("2", exon_index, cds_position)
            
            elif exon_index < 0 or exon_index == None:
                    fails[transcript_id]["lindeboom_exon_index_error"] += 1; failed = True

            else:
                fails[transcript_id]["lindeboom_exon_size_error"] += 1; failed = True

        else:
            fails[transcript_id]["lindeboom_not_found"] += 1; failed = True
            
        return lindeboom_prediction, fails, failed


    def initialize(self):       
        # load genome builds 
        if self.params["hg_build"] == "hg19":
            self.gene_identifier = "uc id"
            hg_build_dir         = "hg19"
            knowngene_fname      = "hg19_knownGene.txt"
            lindeboom_fname      = "hg19_NMDetectiveA_Lindeboom_et_al.v2.gtf"

        if self.params["hg_build"] == "hg38":
            self.gene_identifier = "transcript id"
            hg_build_dir         = "hg38.p14"
            knowngene_fname      = "hg38_knownGene.txt"
            lindeboom_fname      = "hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf"


        # load lindeboom file
        with open(self.params["data_dir"]+self.params["os_sep"]+lindeboom_fname, "r") as f:
            lines = f.readlines()

        # create gene dictionary
        self.lindeboom_dictionary = create_lindeboom_dictionary(lines)
        print("< gene dictionary created.")


        # load knowGene dictionary
        with open(self.params["data_dir"]+self.params["os_sep"]+knowngene_fname, 'r') as _:
            self.knowngene = pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+knowngene_fname, delimiter=",", index_col=False).sort_values(by=["chr", "cdsstart"])


        if len(self.params["genome_index"]) != 0:
            self.knowngene = self.knowngene.iloc[self.params["genome_index"][0]:min(self.params["genome_index"][1], self.knowngene.shape[0])]
            print("< knowngene database shrinked to", self.knowngene.shape[0])

        # load genome
        self.genome = load_split_genome(self.params["data_dir"]+self.params["os_sep"]+hg_build_dir, os_sep=self.params["os_sep"])

        # marked (<-) added on 250427 to integrate expression as features
        self.expressions            = pd.read_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["expressions_fname"], delimiter=",", index_col=False) # <-
        self.expressions["gene_id"] = [self.expressions.iloc[i].loc["gene_id"].split(".")[0] for i in range(self.expressions.shape[0])]

        # marked (<-) added on 250427 to integrate allele frequency as feature
        # self.gnomad_db = gnomAD_DB(self.params["data_dir"], gnomad_version="v4") # <- added

        # marked (<-) added on 250428 to allow flexible usage of params->default_values or params->values (see initialize)
        # default values are prioritized over values
        self.values = {**self.params["default_values"], **{key: self.params["values"][key] for key in self.params["values"] if key not in self.params["default_values"]}} # <-


    def invert(self, seq):
        inverse_seq = ""
        #bases         = ["A", "C", "G", "T"]
        #inverse_bases = ["T", "G", "C", "A"]
        bases         = ["A", "C", "G", "T", "N", "-"]
        inverse_bases = ["T", "G", "C", "A", "N", "-"]
        
        for i in range(len(seq)):
            j = 0
            match = False
            while j < len(bases) and match == False:
                if seq[len(seq)-1-i] == bases[j]:
                    inverse_seq += inverse_bases[j]
                    match = True
                j += 1

        return inverse_seq


    def _predict(self, predictions, thread_index, thread_id):
        print("< thread id:", thread_id, "-", thread_index[0], "/", thread_index[1])

        # iterate over entries in genome dictionary and conduct predictions for coding sequences
        #features  = {**{"transcript id": []}, **{"exon_cds_start": []}, **{feature: [] for feature in self.models[0].feature_names_in_}}

        # marked (<-) added / removed on 250428
        # features  = {**{self.gene_identifier: []}, **{"exon_cds_start": []}, **{feature: [] for feature in self.models[0].feature_names_in_}} # <- removed
        features  = {**{self.gene_identifier: []}, **{"exon_cds_start": []}, **{feature: [] for feature in self.selected_features}} # <- added
        fails     = {}
        
        # feature extraction
        feature_start_time = time.time()
        steps = 0
        for i in range(thread_index[0], thread_index[1]+1, 1):
            if self.knowngene.iloc[i].loc["cdsstart"] < self.knowngene.iloc[i].loc["cdsend"] and self.knowngene.iloc[i].loc["chr"] in self.genome: 
                steps += 1
                features, fails = self.extract_features(features, fails, i)

            if i % 100 == 0: print("< thread id", thread_id, i, thread_index[1])
            
        feature_end_time = time.time()
        size_test        = pd.DataFrame({"feature": [col for col in features], "sizes": [len(features[col]) for col in features]})
        if size_test.drop_duplicates(subset=["sizes"]).shape[0] != 1:
            print("< error occurred. inconsisent sizes in feature container. thread id:", thread_id); print(size_test)
        
        features = pd.DataFrame(features)

        # prediction based on extracted features       
        all_preds        = np.zeros((0))

        for i in range(len(self.models)):
            # marked (<-) added / removed on 250428 to select only features required by the model
            # preds = self.models[i].predict(features.loc[:,features.columns.tolist()[2::]]) # <- removed
            preds = self.models[i].predict(features.loc[:, self.models[0].feature_names_in_]) # <- added
            if all_preds.shape[0] > 0: all_preds = np.column_stack((all_preds, preds))
            else:                      all_preds = preds

        if len(all_preds.shape) > 1: all_preds = np.average(all_preds, axis=1)

        last_transcript_id  = None
        predictions_by_exon = {}
        test_size           = 0
        for i in range(len(all_preds)):
            # marked (<-) added / removed on 250430 to avoid explitit output of number format on MacOS (e.g. np.int64(number))
            # predictions.at[features.iloc[i].loc[self.gene_identifier], "predictions"].append(round(all_preds[i], 4)) # <- removed
            predictions.at[features.iloc[i].loc[self.gene_identifier], "predictions"].append(float(round(all_preds[i], 4))) # <- added

            #if last_transcript_id != None and last_transcript_id != features.iloc[i].loc["transcript id"]:
            if last_transcript_id != None and last_transcript_id != features.iloc[i].loc[self.gene_identifier]:
                predictions.at[last_transcript_id, "fails"]               = fails[last_transcript_id]
                predictions.at[last_transcript_id, "predictions_by_exon"] = predictions_by_exon
                # marked (<-) added / removed on 250430 to avoid explitit output of number format on MacOS (e.g. np.int64(number))
                # predictions_by_exon                                       = {features.iloc[i].loc["exon_cds_start"]: [round(all_preds[i], 4)]} # <- removed
                predictions_by_exon                                       = {int(features.iloc[i].loc["exon_cds_start"]): [float(round(all_preds[i], 4))]} # <- added

                if len(predictions.loc[last_transcript_id].loc["predictions"]) != test_size:
                    print("< inconsistent no. of predictions for", last_transcript_id, "thread id", thread_id, len(predictions.loc[last_transcript_id].loc["predictions"]), "/", test_size)
                
                test_size = 1

            else:
                # marked (<-) added / removed on 250430 to avoid explitit output of number format on MacOS (e.g. np.int64(number))
                # if features.iloc[i].loc["exon_cds_start"] in predictions_by_exon: predictions_by_exon[features.iloc[i].loc["exon_cds_start"]].append(round(all_preds[i], 4)); test_size += 1 # <- removed
                # else:                                                             predictions_by_exon[features.iloc[i].loc["exon_cds_start"]] = [round(all_preds[i], 4)]; test_size += 1 # <- removed
                if features.iloc[i].loc["exon_cds_start"] in predictions_by_exon: predictions_by_exon[int(features.iloc[i].loc["exon_cds_start"])].append(float(round(all_preds[i], 4))); test_size += 1 # <- added
                else:                                                             predictions_by_exon[int(features.iloc[i].loc["exon_cds_start"])] = [float(round(all_preds[i], 4))]; test_size += 1 # <- added

            #last_transcript_id = features.iloc[i].loc["transcript id"]
            last_transcript_id = features.iloc[i].loc[self.gene_identifier]

        # register last block
        predictions.at[last_transcript_id, "fails"]               = fails[last_transcript_id]
        predictions.at[last_transcript_id, "predictions_by_exon"] = predictions_by_exon

        if len(predictions.loc[last_transcript_id].loc["predictions"]) != test_size:
            print("< inconsistent no. of predictions for", last_transcript_id, len(predictions.loc[last_transcript_id].loc["predictions"]), "/", test_size)
        
        model_end_time = time.time()
        print("< thread id:", thread_id, ", feature time", feature_end_time-feature_start_time, ", steps", steps, ", model time", model_end_time-feature_start_time)
        
        # back testing for consistency of exon-resolved predictions and aggregated predictions
        return


    def predict(self):
        if self.knowngene.shape[0]-self.knowngene.drop_duplicates(subset=[self.gene_identifier]).shape[0] != 0:
            print("< error. gene identifier", self.gene_identifier, "is not unique.")
            exit()

        target_index = [i for i in range(self.knowngene.shape[0]) if self.knowngene.iloc[i].loc["cdsstart"] < self.knowngene.iloc[i].loc["cdsend"] and self.knowngene.iloc[i].loc["chr"] in self.genome]
        predictions  = pd.DataFrame({"gene id":               self.knowngene.iloc[target_index]["gene id"].tolist(),
                                     "transcript id":         self.knowngene.iloc[target_index]["transcript id"].tolist(),
                                     "uc id":                 self.knowngene.iloc[target_index]["uc id"].tolist(),
                                     "chr":                   self.knowngene.iloc[target_index]["chr"].tolist(),
                                     "strand":                self.knowngene.iloc[target_index]["strand"].tolist(),
                                     "predictions":           [[] for _ in target_index],
                                     "predictions_by_exon":   [{} for _ in target_index],
                                     "fails":                 [[] for _ in target_index]})
        
        predictions.index = [transcript_id.split(".")[0] for transcript_id in predictions[self.gene_identifier]]

        thread_index = split_index(self.knowngene.shape[0], self.params["threads"])
        threads      = []

        for i in range(len(thread_index)):
            thread = threading.Thread(target=self._predict, args=(predictions, thread_index[i], i))
            threads.append(thread)
            thread.start()

        for i, thread in enumerate(threads):
            thread.join()

        if len(self.params["genome_index"]) > 0:
            predictions.to_csv(path_or_buf=self.params["data_dir"]+self.params["os_sep"]+self.params["hg_build"]+"_"+str(self.params["genome_index"][0])+"_"+str(self.params["genome_index"][1])+"_NMD_scores.txt", sep=",", index=False)
        
        else:
            predictions.to_csv(path_or_buf=self.params["data_dir"]+self.params["os_sep"]+self.params["hg_build"]+"_NMD_scores.txt", sep=",", index=False)
