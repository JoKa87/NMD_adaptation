import copy
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\shared")

from shared_utils import *


class Mask_predictions():
    def __init__(self, params, mutation_stats):
        self.params              = params

        self.errors              = {"no_gene_factor_found": [], "position_not_covered": []}
        self.masking_stats       = pd.DataFrame()
        self.mutation_stats      = mutation_stats

        if "apply_mutation_stats" in self.params:
            if self.params["apply_mutation_stats"] == True and self.params["apply_3mers"] == False:
                self.probability_trajectory = {
                                            **{"pairs":       {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["pairs"][key1]}
                                                                for key1 in self.mutation_stats["pairs"]}},
                                            **{"del-1":       {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["del-1"][key1]}
                                                                for key1 in self.mutation_stats["del-1"]}},
                                            **{"ins+1":       {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["ins+1"][key1]}
                                                                for key1 in self.mutation_stats["ins+1"]}},
                                            **{"topology":    {key: [0 for _ in range(100)] for key in ["del-1", "frameshift", "ins+1", "missense", "nonsense", "total"]}}
                                            }
                
            elif self.params["apply_mutation_stats"] == True and self.params["apply_3mers"] == True:
                self.probability_trajectory = {
                                            **{"pairs_3mers": {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["pairs_3mers"][key1]}
                                                                for key1 in self.mutation_stats["pairs_3mers"]}},
                                            **{"del-1_3mers": {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["del-1_3mers"][key1]}
                                                                for key1 in self.mutation_stats["del-1_3mers"]}},
                                            **{"ins+1_3mers": {key1: {key2: [0 for _ in range(100)] for key2 in self.mutation_stats["ins+1_3mers"][key1]}
                                                                for key1 in self.mutation_stats["ins+1_3mers"]}},
                                            **{"topology":    {key: [0 for _ in range(100)] for key in ["del-1", "frameshift", "ins+1", "missense", "nonsense", "total"]}}
                                            }
                
            else:
                self.probability_trajectory = {
                                            **{"topology":    {key: [0 for _ in range(100)] for key in ["del-1", "frameshift", "ins+1", "missense", "nonsense", "total"]}}
                                            }
            

    def _analyze_probability_trajectories(self, accumulated_probability_trajectory, mutation_type):
        topology = [key for key in self.probability_trajectory if mutation_type in key][0]

        accumulated_trajectory_correlation = pd.DataFrame({"pearson-r": [None for _ in self.probability_trajectory[topology]["total"]],
                                                           "pearson-p": [None for _ in self.probability_trajectory[topology]["total"]]},
                                                           index=self.probability_trajectory[topology]["total"])
        
        for row in accumulated_trajectory_correlation.index:
            pearson = scipy.stats.pearsonr(self.probability_trajectory[topology]["total"][row], accumulated_probability_trajectory[topology])
            accumulated_trajectory_correlation.at[row, "pearson-r"] = pearson.statistic
            accumulated_trajectory_correlation.at[row, "pearson-p"] = pearson.pvalue

        trajectory_correlation = pd.DataFrame({col: [None for _ in self.probability_trajectory[topology]["total"]]
                                               for col in self.probability_trajectory[topology]["total"]},
                                               index=self.probability_trajectory[topology]["total"])
        
        for col in trajectory_correlation.columns:
            for row in trajectory_correlation.index:
                trajectory_correlation.at[row, col] = scipy.stats.pearsonr(self.probability_trajectory[topology]["total"][row],
                                                                           self.probability_trajectory[topology]["total"][col]).pvalue
        
        return accumulated_trajectory_correlation, trajectory_correlation
    

    def analyze_probability_trajectories(self):
        # calculate accumulated probability trajectories
        accumulated_probability_trajectory = {key: [0 for _ in range(100)] for key in self.probability_trajectory if key != "topology"}

        for key1 in self.probability_trajectory:
            if "topology" not in key1:
                for key2 in self.probability_trajectory[key1]["total"]:
                    for i, value in enumerate(self.probability_trajectory[key1]["total"][key2]):
                        accumulated_probability_trajectory[key1][i] += value

        for key in accumulated_probability_trajectory:
            max_value = np.max(accumulated_probability_trajectory[key])
            accumulated_probability_trajectory[key] = [value/max_value for value in accumulated_probability_trajectory[key]]
            plt.plot(np.arange(len(accumulated_probability_trajectory[key])), accumulated_probability_trajectory[key], label=key)
        
        plt.title("accumulated probability trajectory")
        plt.show()


        # calculate trajectory correlation matrix
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)

        accumulated_trajectory_correlation, trajectory_correlation = self._analyze_probability_trajectories(accumulated_probability_trajectory, "del")
        print("< deletion correlation")
        print(accumulated_trajectory_correlation)

        accumulated_trajectory_correlation, trajectory_correlation = self._analyze_probability_trajectories(accumulated_probability_trajectory, "ins")
        print("< insertion correlation")
        print(accumulated_trajectory_correlation)

        # plot individual base/3mer trajectories for each given category 
        for key1 in self.probability_trajectory:
            if "topology" not in key1:

                for key2 in self.probability_trajectory[key1]["total"]:
                    plt.plot(np.arange(len(self.probability_trajectory[key1]["total"][key2])), self.probability_trajectory[key1]["total"][key2], label=key2)

                plt.title(key1)
                plt.legend()
                plt.show()


    def calculate_driver_genes(self):
        # select entries in masking stats relevant for downstream analysis
        self.masking_stats = self.masking_stats[[col for col in self.masking_stats.columns
                                                 if col in ["block", "transcript id", "gene symbol"] or len([mask for mask in self.params["masks"] if mask in col]) > 0]]
        self.masking_stats = self.masking_stats[self.masking_stats["transcript id"] != "total"]
        
        # calculating projections vs. observations, observations determined from mutation stats
        for mask in self.params["masks"]:
            self.masking_stats["observed "+mask] = [self.mutation_stats["observed"]["total"][mask][transcript_id] for transcript_id in self.masking_stats["transcript id"]]        
            self.masking_stats[mask+" pvalue"]   = [None for _ in self.masking_stats["transcript id"]]
            #observed_mutations[mask]             = self.masking_stats["observed "+mask].sum() # np.sum([self.mutation_stats["observed"]["total"][mask][transcript_id] for transcript_id in self.mutation_stats["observed"]["total"][mask]])

        #self.masking_stats["gene symbol"]    = [seqs[seqs["transcript id"] == gene].iloc[0].loc["gene symbol"] for gene in mp.masking_stats["transcript id"]]
        self.masking_stats["observed total"] = self.masking_stats[["observed "+mask for mask in self.params["masks"]]].sum(axis=1)
        self.masking_stats["proj. total"]    = self.masking_stats[["proj. "+mask for mask in self.params["masks"]]].sum(axis=1)
        self.masking_stats["total"]          = self.masking_stats[[mask for mask in self.params["masks"]]].sum(axis=1)
        self.masking_stats["total pvalue"]   = [None for _ in self.masking_stats["transcript id"]]
        self.masking_stats                   = self.masking_stats.reset_index()

        self.evaluate_masking_stats()
        #mp.masking_stats = mp.masking_stats[mp.masking_stats["transcript id"] != "total"]

        # calculating driver gene probabilities
        bar = IncrementalBar(set_bar("calculating driver gene probabilities"), max=self.masking_stats.shape[0])
        for i in range(self.masking_stats.shape[0]):
            for mask in [*self.params["masks"], "total"]:
                expected_probability = self.masking_stats.iloc[i].loc["proj. "+mask] / len(self.mutation_stats["cases"]["total"])
                observed_probability = self.masking_stats.iloc[i].loc["observed "+mask] / len(self.mutation_stats["cases"]["total"])

                if observed_probability >= expected_probability:
                    self.masking_stats.at[self.masking_stats.index[i], mask+" pvalue"] = 1-scipy.stats.binom.cdf(self.masking_stats.iloc[i].loc["observed "+mask],
                                                                                                                 len(self.mutation_stats["cases"]["total"]), expected_probability)

                if observed_probability < expected_probability:
                    self.masking_stats.at[self.masking_stats.index[i], mask+" pvalue"] = scipy.stats.binom.cdf(self.masking_stats.iloc[i].loc["observed "+mask],
                                                                                                               len(self.mutation_stats["cases"]["total"]), expected_probability)

            bar.next()
        bar.finish()

        for mask in [*self.params["masks"], "total"]:
            pvalues                            = self.masking_stats[mask+" pvalue"].tolist()
            self.masking_stats[mask+" pvalue"] = [*sm.stats.fdrcorrection(pvalues[0:len(pvalues)-1], alpha=0.05)[1], pvalues[0]]

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        self.masking_stats = self.masking_stats.sort_values(by="total pvalue")
        self.masking_stats = self.masking_stats[self.masking_stats["total pvalue"] < 0.01]
        print(self.masking_stats)
        self.masking_stats.to_csv(self.params["data_dir"]+self.params["os_sep"]+self.params["fname"].replace(".txt", "_masking_stats.txt"), index=False, sep=",")
                    
    
    def evaluate_masking_stats(self):
        for mask in [*self.params["masks"], "total"]:
            self.masking_stats[mask+" distance"] = [self.masking_stats.iloc[i].loc["observed "+mask]-self.masking_stats.iloc[i].loc["proj. "+mask]
                                                    for i in range(self.masking_stats.shape[0])]
        
        masking_stats = pd.DataFrame({
                                      **{"block"               : ["total"],
                                         "transcript id"       : ["total"],
                                         "gene symbol"         : ["total"]},
                                      **{mask                  : [self.masking_stats[mask].mean()] for mask in [*self.params["masks"], "total"]},
                                      **{"proj. "+mask         : [self.masking_stats["proj. "+mask].sum()] for mask in [*self.params["masks"], "total"]},
                                      **{"observed "+mask      : [self.masking_stats["observed "+mask].sum()]}, # self.masking_stats["observed "+mask][observed_variants[mask].sum()] for mask in [*self.params["masks"], "total"]},
                                      **{mask+" distance"      : [self.masking_stats[mask+" distance"].mean()] for mask in [*self.params["masks"], "total"]},
                                     })
        
        self.masking_stats = pd.concat([self.masking_stats, masking_stats]).sort_values(by="observed total", ascending=False)
        pd.set_option('display.max_columns', None)
        #print(self.masking_stats)


    def __calculate_probabilities(self, cds, it, gene_factor, pair_target, mode="-", limits=None):
        is_stop = False
        prob    = 0

        hist_index = int(100*(it-3)/(len(cds)-3))

        if limits == None:
            frames = [0, 1, 2]

        else:
            frames = []

            if it-limits[0] < 0:
                for i in range(limits[0]-it, 3):
                    frames.append(i)

            elif limits[1]-it < 3:
                for i in range(0, limits[1]-it+1):
                    frames.append(i)

            else:
                frames = [0, 1, 2]

        # test whether 'it' is in-frame
        if it-3*int(it/3) != 0:
            print("< error @__calculate_probabilities. iterator is not in-frame")
            exit()

        if mode == "-":
            for i in range(3):
                quadruplett = cds[it:it+4]
                wt_base     = quadruplett[i]
                wt_context  = cds[it+i-1:it+i+2]
                triplett    = "".join(b for j, b in enumerate(quadruplett) if j != i)

                if i in frames and triplett in ["TAA", "TAG", "TGA"]:
                    if self.params["apply_mutation_stats"] == True:
                        if self.params["apply_3mers"] == False:
                            prob += self.mutation_stats["del-1"][pair_target][wt_base]
                            self.probability_trajectory["del-1"][pair_target][wt_base][hist_index] += gene_factor*self.mutation_stats["del-1"][pair_target][wt_base]
                        
                        if self.params["apply_3mers"] == True:
                            prob += self.mutation_stats["del-1_3mers"][pair_target][wt_context]
                            self.probability_trajectory["del-1_3mers"][pair_target][wt_context][hist_index] += gene_factor*self.mutation_stats["del-1_3mers"][pair_target][wt_context]

                    is_stop = True

        if "+" in mode:
            for i in range(3):
                for base in ["A", "G", "T"]:
                    if i == 0:
                        base_context = cds[it-1]+base+cds[it]
                        triplett     = base+cds[it:it+2]
                    
                    if i == 1:
                        base_context = cds[it]+base+cds[it+1]
                        triplett     = cds[it]+base+cds[it+1]

                    if i == 2:
                        base_context = cds[it+1]+base+cds[it+2]
                        triplett     = cds[it:it+2]+base

                    if i in frames and triplett in ["TAA", "TAG", "TGA"]:
                        if self.params["apply_mutation_stats"] == True:
                            if self.params["apply_3mers"] == False:
                                prob += self.mutation_stats["ins+1"][pair_target][base]
                                self.probability_trajectory["ins+1"][pair_target][base][hist_index] += gene_factor*self.mutation_stats["ins+1"][pair_target][base]
                                
                            if self.params["apply_3mers"] == True:
                                prob += self.mutation_stats["ins+1_3mers"][pair_target][base_context]
                                self.probability_trajectory["ins+1_3mers"][pair_target][base_context][hist_index] += gene_factor*self.mutation_stats["ins+1_3mers"][pair_target][base_context]
                        #print(i, base, cds[it:it+3], triplett, frames, prob)
                        #print(it, hist_index, i, base_context, cds[it:it+3], frames, gene_factor*self.mutation_stats["ins+1_3mers"][pair_target][base_context])

                        is_stop = True

        return (prob, is_stop)


    def _calculate_probabilities(self, cds, it, last_stop_index, gene_factor, pair_target, mode="-"):
        hist_index = int(100*(it-3)/(len(cds)-3))

        if mode == "-":
            if self.params["apply_3mers"] == False: target_key = "del-1"
            if self.params["apply_3mers"] == True:  target_key = "del-1_3mers"

            pre_del1   = {key: np.sum(self.probability_trajectory[target_key][pair_target][key])
                          for key in self.probability_trajectory[target_key][pair_target]}

            local_del1 = [self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="-", limits=(last_stop_index, it))[0]
                          for i in range(3*int(last_stop_index/3), it+1, 3)]
            
            post_del1  = {key: np.sum(self.probability_trajectory[target_key][pair_target][key])-pre_del1[key]
                          for key in self.probability_trajectory[target_key][pair_target]}
            
            if self.params["apply_3mers"] == False:
                total_del1 = [self.mutation_stats["del-1"][pair_target][cds[i]] for i in range(last_stop_index, it+1)]

                for i in range(last_stop_index, it+1):
                    self.probability_trajectory["del-1"][pair_target][cds[i]][hist_index] += gene_factor*self.mutation_stats["del-1"][pair_target][cds[i]]

            if self.params["apply_3mers"] == True:
                total_del1 = [self.mutation_stats["del-1_3mers"][pair_target][cds[i-1:i+2]] for i in range(last_stop_index, it+1)]

                for i in range(last_stop_index, it+1):
                    self.probability_trajectory["del-1_3mers"][pair_target][cds[i-1:i+2]][hist_index] += gene_factor*self.mutation_stats["del-1_3mers"][pair_target][cds[i-1:i+2]]
            
            for key in post_del1:
                self.probability_trajectory[target_key][pair_target][key][hist_index] -= post_del1[key]

            prob = np.sum(total_del1)-np.sum(local_del1)


        if mode == "+":
            if self.params["apply_3mers"] == False: target_key = "ins+1"
            if self.params["apply_3mers"] == True:  target_key = "ins+1_3mers"

            pre_ins1   = {key: np.sum(self.probability_trajectory[target_key][pair_target][key])
                          for key in self.probability_trajectory[target_key][pair_target]}
            
            local_ins1 = [self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="+", limits=(last_stop_index, it))[0]
                          for i in range(3*int(last_stop_index/3), it+1, 3)]
            
            post_ins1  = {key: np.sum(self.probability_trajectory[target_key][pair_target][key])-pre_ins1[key]
                          for key in self.probability_trajectory[target_key][pair_target]}

            if self.params["apply_3mers"] == False:
                total_ins1 = [self.mutation_stats["ins+1"][pair_target]["total"] for _ in range(last_stop_index, it+1)]

                for i in range(last_stop_index, it+1):
                    for base in ["A", "C", "G", "T"]:
                        self.probability_trajectory["ins+1"][pair_target][base][hist_index] += gene_factor*self.mutation_stats["ins+1"][pair_target][base]
            
            if self.params["apply_3mers"] == True:
                total_ins1 = [np.sum([self.mutation_stats["ins+1_3mers"][pair_target][cds[i-1]+base+cds[i]] for base in ["A", "C", "G", "T"]])
                                      for i in range(last_stop_index, it+1)]
                
                for i in range(last_stop_index, it+1):
                    for base in ["A", "C", "G", "T"]:
                        self.probability_trajectory["ins+1_3mers"][pair_target][cds[i-1]+base+cds[i]][hist_index] += gene_factor*self.mutation_stats["ins+1_3mers"][pair_target][cds[i-1]+base+cds[i]]

            for key in post_ins1:
                self.probability_trajectory[target_key][pair_target][key][hist_index] -= post_ins1[key]
                #print(key, self.probability_trajectory[target_key][pair_target][key][hist_index], post_ins1[key])
            
            prob = np.sum(total_ins1)-np.sum(local_ins1)

            if len(total_ins1) % 3 != 0: # new
                print("< size error @_calculate_probabilities ("+str(len(total_ins1))+")") # new
            

        if prob < 0:
            print("< error @_calculate_probabilities. negative probability detected.")
            print(last_stop_index, it, mode)
            exit()

        return prob


    # should be reconsidered under the aspect that for variant predictions, last position of stop codon was used (i+2)
    def calculate_probabilities(self, cds, transcript_id, gene_target, pair_target, variant_targets={}, show=False):
        prob_stats      = {"del-1": 0, "frameshift": 0, "ins+1": 0, "missense": 0, "nonsense": 0, "total": 0}
        last_prob_stats = {"del-1": 0, "frameshift": 0, "ins+1": 0, "missense": 0, "nonsense": 0, "total": 0}
        probabilities   = [np.nan for _ in range(len(cds)-3)]
        
        # in MSK data, on rare occasions deviations can occur between cancer-wise and all cancer runs as for some patients,
        # mutation data are assigned to one cancer type, although multiple primary cancer types are assigned with the patient 
        if self.params["apply_mutation_stats"] == True:
            if transcript_id in self.mutation_stats["genes"][gene_target]:
                gene_factor = self.params["mutation_stats_scale"]*self.mutation_stats["genes"][gene_target][transcript_id]

            elif self.params["id_filter"] != "TCGA":
                print("< no gene factor found for", transcript_id+" in gene target "+gene_target+". gene factor set to zero")
                self.errors["no_gene_factor_found"].append(gene_target+","+transcript_id)
                gene_factor = 0

            else:
                print("< no gene factor found for", transcript_id+" in gene target "+gene_target)
                exit()

        else:
            gene_factor = 1

        # set gene factor to no. of PTC mutations if chosen (250821)
        if self.params["mutation_stats_ptc_weights"] == True and "FEATURE:ptc cds position" in variant_targets:
            gene_factor = 0
            for abs_pos in variant_targets["FEATURE:ptc cds position"]:
                for _ in variant_targets["FEATURE:ptc cds position"][abs_pos]:
                    gene_factor += 1

        elif self.params["apply_mutation_weights"] == True: # recheck!
            gene_factor = self.params["mutation_stats_scale"]*self.mutation_stats["genes"][gene_target][transcript_id]

        # test for start codon
        if cds[0:3] != "ATG":
            print("< start codon not found @mask_predictions")
            exit()

        test_index = {"global_frameshift-1": [], "local_frameshift-1": []}
        last_stop_index_1 = 3; last_stop_index_2 = 3
        for i in range(3, len(cds)-3, 3):
            current_probs = [0 for _ in range(3)]
            hist_index    = int(100*(i-3)/(len(cds)-3))

            # detect missense or nonsense candidates
            if "missense" in self.params["masks"] or "nonsense" in self.params["masks"]:
                for j in range(3):                    
                    for base in ["A", "C", "G", "T"]:
                        triplett    = cds[i:i+3]
                        wt_base     = triplett[j]
                        wt_context  = cds[i+j-1:i+j+2]
                        triplett    = "".join(b if k != j else base for k, b in enumerate(triplett))

                        if "nonsense" in self.params["masks"] and triplett in ["TAA", "TAG", "TGA"]:
                            if self.params["apply_mutation_stats"] == False:
                                current_probs[2] = gene_factor # 1
                            
                            else:
                                if self.params["apply_3mers"] == False:
                                    current_prob = gene_factor*(self.mutation_stats["pairs"][pair_target][wt_base+base])
                                    self.probability_trajectory["pairs"][pair_target][wt_base+base][hist_index] += gene_factor*self.mutation_stats["pairs"][pair_target][wt_base+base]

                                if self.params["apply_3mers"] == True:
                                    current_prob = gene_factor*self.mutation_stats["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]]
                                    self.probability_trajectory["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]][hist_index] += gene_factor*self.mutation_stats["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]]

                                current_probs[2]       += current_prob
                                prob_stats["nonsense"] += current_prob/self.params["mutation_stats_scale"]
                                if show == True: print("nonsense-"+str(j+1), i, cds[i:i+3], gene_factor, "/", current_prob/gene_factor)

                        
                        if "missense" in self.params["masks"] and triplett not in ["TAA", "TAG", "TGA"] and genetic_code[cds[i:i+3]] != genetic_code[triplett]:
                            if self.params["apply_mutation_stats"] == False:
                                current_probs[j] = gene_factor # 1
                            
                            else:
                                if self.params["apply_3mers"] == False:
                                    current_prob = gene_factor*(self.mutation_stats["pairs"][pair_target][wt_base+base])
                                    self.probability_trajectory["pairs"][pair_target][wt_base+base][hist_index] += gene_factor*self.mutation_stats["pairs"][pair_target][wt_base+base]

                                if self.params["apply_3mers"] == True:
                                    current_prob = gene_factor*self.mutation_stats["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]]
                                    self.probability_trajectory["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]][hist_index] += gene_factor*self.mutation_stats["pairs_3mers"][pair_target][wt_context[0:2]+base+wt_context[2]]

                                current_probs[j]       += current_prob
                                prob_stats["missense"] += current_prob/self.params["mutation_stats_scale"]
                                if show == True: print("missense-"+str(j+1), i, cds[i:i+3], triplett, gene_factor, "/", current_prob/gene_factor)

            if "frameshift" in self.params["masks"]:
                # detect candidates for frameshift -1
                if i+4 < len(cds) and cds[i+1:i+4] in ["TAA", "TAG", "TGA"]:
                    if self.params["apply_mutation_stats"] == False:
                        current_probs[2] = gene_factor # 1

                    else:
                        current_prob = self._calculate_probabilities(cds, i, last_stop_index_1, gene_factor, pair_target, mode="-")
                        test_index["global_frameshift-1"].append(i+2)

                        if current_prob > 0:
                            current_prob             *= gene_factor
                            current_probs[2]         += current_prob
                            prob_stats["del-1"]      += current_prob/self.params["mutation_stats_scale"]
                            prob_stats["frameshift"] += current_prob/self.params["mutation_stats_scale"]

                            if show == True: print("global frameshift-1", i, cds[i:i+5], last_stop_index_1, gene_factor, "/", current_prob/gene_factor)
                            last_stop_index_1 = i+1
                            

                # detect candidates for frameshift +1
                #if i+5 < len(cds) and cds[i+2:i+5] in ["TAA", "TAG", "TGA"]:
                if cds[i-1:i+2] in ["TAA", "TAG", "TGA"]:
                    if self.params["apply_mutation_stats"] == False:
                        current_probs[2] = gene_factor # 1

                    else:
                        #current_prob = self._calculate_probabilities(cds, i+1, last_stop_index_2, gene_factor, pair_target, mode="+")
                        current_prob = self._calculate_probabilities(cds, i-1, last_stop_index_2, gene_factor, pair_target, mode="+")

                        if current_prob > 0:
                            current_prob             *= gene_factor
                            current_probs[2]         += current_prob
                            prob_stats["ins+1"]      += current_prob/self.params["mutation_stats_scale"]
                            prob_stats["frameshift"] += current_prob/self.params["mutation_stats_scale"]
                            
                            #if show == True: print("global frameshift-2", i, cds[i:i+5], last_stop_index_2, gene_factor, "/", current_prob/gene_factor)
                            #last_stop_index_2 = i+2
                            if show == True: print("global frameshift-2", i, cds[i-3:i+3], last_stop_index_2, gene_factor, "/", current_prob/gene_factor)
                            last_stop_index_2 = i

                # deletions at positions 1, 2, or 3
                if self.params["apply_mutation_stats"] == False:
                    if self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="-")[1] == True:
                        current_probs[2] = gene_factor # 1

                else:
                    current_prob, is_stop = self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="-")
                    if is_stop == True: test_index["local_frameshift-1"].append(i+2)

                    if current_prob > 0:
                        current_prob             *= gene_factor
                        current_probs[2]         += current_prob
                        prob_stats["del-1"]      += current_prob/self.params["mutation_stats_scale"]
                        prob_stats["frameshift"] += current_prob/self.params["mutation_stats_scale"]
                        if show == True: print("local frameshift-1", i, cds[i:i+5], gene_factor, "/", current_prob/gene_factor)

                # insertions at positions 1, 2, or 3
                if self.params["apply_mutation_stats"] == False:
                    if self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="+")[1] == True:
                        current_probs[2] = gene_factor # 1

                else:
                    current_prob, _ = self.__calculate_probabilities(cds, i, gene_factor, pair_target, mode="+")

                    if current_prob > 0:
                        current_prob             *= gene_factor
                        current_probs[2]         += current_prob
                        prob_stats["ins+1"]      += current_prob/self.params["mutation_stats_scale"]
                        prob_stats["frameshift"] += current_prob/self.params["mutation_stats_scale"]
                        if show == True: print("local frameshift-2", i, cds[i:i+5], gene_factor, "/", current_prob/gene_factor)
            
            prob_stats["total"] += (current_probs[0]+current_probs[1]+current_probs[2])/self.params["mutation_stats_scale"]

            for j in range(len(current_probs)):
                if current_probs[j] > 0:
                    probabilities[i+j] = current_probs[j]

                    for key in prob_stats:
                        if prob_stats[key] > last_prob_stats[key]:
                            self.probability_trajectory["topology"][key][hist_index] += 1
                    
                    self.probability_trajectory["topology"]["total"][hist_index] += 1

            last_prob_stats = {key: copy.deepcopy(value) for key, value in prob_stats.items()}


            if self.params["apply_mutation_stats"] == True and round(prob_stats["total"], 6) != round(prob_stats["frameshift"]+prob_stats["missense"]+prob_stats["nonsense"], 6):
                print("< count error1 occurred ("+str(prob_stats["total"])+"/"+str(prob_stats["frameshift"]+prob_stats["missense"]+prob_stats["nonsense"])+")")
                print(prob_stats)

            if self.params["apply_mutation_stats"] == True and round(prob_stats["frameshift"], 6) != round(prob_stats["del-1"]+prob_stats["ins+1"], 6):
                print("< count error2 occurred ("+str(prob_stats["total"])+"/"+str(prob_stats["frameshift"]+prob_stats["missense"]+prob_stats["nonsense"])+")")
                print(prob_stats)

        
        # test for redundant calling of frameshift mutations
        for i in test_index["global_frameshift-1"]:
            if i not in test_index["local_frameshift-1"]:
                print("< global frameshift index "+str(i)+" not found in local frameshift index @mask_predictions ("+transcript_id+")")
                input("x")

        if self.params["apply_mutation_stats"] == True:
            if show == True:
                print("prob stats", prob_stats)
                print("values", np.count_nonzero(~np.isnan(probabilities)), "ratio", np.count_nonzero(~np.isnan(probabilities))/len(probabilities), len(cds), "/", len(probabilities),
                    "total score", prob_stats["total"], "frameshift score", round(prob_stats["total"]-prob_stats["nonsense"], 4), "nonsense score", round(prob_stats["nonsense"], 4),
                    "frameshift prob / TCGA", round(len(self.mutation_stats["cases"][gene_target])*(prob_stats["frameshift"])/self.mutation_stats["avg. gene factor"][gene_target], 4),
                    "nonsense prob / TCGA", round(len(self.mutation_stats["cases"][gene_target])*prob_stats["nonsense"]/self.mutation_stats["avg. gene factor"][gene_target], 4))
                input("x")


            # register results
            masking_stats = pd.DataFrame({"block"               : [gene_target],
                                          "transcript id"       : [transcript_id],
                                          "del-1"               : [round(prob_stats["del-1"], 8)],
                                          "frameshift"          : [round(prob_stats["frameshift"], 8)],
                                          "ins+1"               : [round(prob_stats["ins+1"], 8)],
                                          "missense"            : [round(prob_stats["missense"], 8)],
                                          "nonsense"            : [round(prob_stats["nonsense"], 8)],
                                          "proj. del-1"         : [(round(len(self.mutation_stats["cases"][gene_target])
                                                                    *prob_stats["del-1"]/self.mutation_stats["avg. gene factor"][gene_target], 8))],
                                          "proj. frameshift"    : [(round(len(self.mutation_stats["cases"][gene_target])
                                                                    *prob_stats["frameshift"]/self.mutation_stats["avg. gene factor"][gene_target], 8))],
                                          "proj. ins+1"         : [(round(len(self.mutation_stats["cases"][gene_target])
                                                                    *prob_stats["ins+1"]/self.mutation_stats["avg. gene factor"][gene_target], 8))],
                                          "proj. missense"      : [(round(len(self.mutation_stats["cases"][gene_target])
                                                                    *prob_stats["missense"]/self.mutation_stats["avg. gene factor"][gene_target], 8))],
                                          "proj. nonsense"      : [(round(len(self.mutation_stats["cases"][gene_target])
                                                                    *prob_stats["nonsense"]/self.mutation_stats["avg. gene factor"][gene_target], 8))],
                                          })
            
            if self.masking_stats.shape[0] > 0: self.masking_stats = pd.concat([self.masking_stats, masking_stats])
            else:                               self.masking_stats = masking_stats

        return probabilities


    def mask_predictions(self, predictions, predictions_by_exon, variant_targets, cds, transcript_id, gene_target, pair_target, show=False):
        probabilities         = self.calculate_probabilities(cds, transcript_id, gene_target, pair_target, variant_targets=variant_targets, show=show)
        probabilities_by_exon = {key: [np.nan for _ in range(len(predictions_by_exon[key]))] for key in predictions_by_exon}

        # test for recognition of observed variant (a small number of deviations can occurr as mismatching reference bases are not used
        # as filter criterion anymore during ptc preparation)
        removed_index = {abs_pos: [] for abs_pos in variant_targets["FEATURE:ptc cds position"]} # new
        if "FEATURE:ptc cds position" in variant_targets:
            for abs_pos in variant_targets["FEATURE:ptc cds position"]:
                for i, pos in enumerate(variant_targets["FEATURE:ptc cds position"][abs_pos]):
                    if pd.isna(probabilities[int(pos)]) == True:
                        if self.params["remove_inconsistencies"] == False:
                            print("< error. predicted position not recognized @mask_predictions ("+transcript_id+": "+str(pos)+", "+cds[int(pos)-2:int(pos)+1]+")")
                            exit()

                        else:
                            print("< predicted position not recognized @mask_predictions ("+transcript_id+": "+str(pos)+", "+cds[int(pos)-2:int(pos)+1]+"). position masked.")
                            self.errors["position_not_covered"].append(gene_target+","+transcript_id+","+str(pos))
                            #probabilities[int(pos)] = np.nan
                            # mark for later removal
                            removed_index[abs_pos].append(i) # new
        
        # remove by marked index
        for key in variant_targets:
            for abs_pos in removed_index:
                for i in sorted(removed_index[abs_pos], reverse=True):
                    del variant_targets[key][abs_pos][i]

        if self.params["mutation_stats_ptc_weights"] == True: # recheck!
            stop_codon_candidates = np.sum([1 for prob in probabilities if pd.isna(prob) == False])
            probabilities         = [prob/stop_codon_candidates if pd.isna(prob) == False else prob for prob in probabilities]


        step = 0
        for key in predictions_by_exon:
            for i in range(len(predictions_by_exon[key])):
                if pd.isna(probabilities[step]) == True:
                    predictions[step]           = np.nan
                    predictions_by_exon[key][i] = np.nan
                
                if pd.isna(probabilities[step]) == False: 
                    probabilities_by_exon[key][i] = probabilities[step]

                step += 1


        # test section
        if step != len(predictions):
            print("< dimension error1 occurred @mask_predictions ("+transcript_id+": "+step+" / "+str(len(predictions))+")")
            exit()

        if self.params["apply_mutation_stats"] == True:
            # test for inconsistent Nones
            test1 = [i for i in range(len(predictions)) if pd.isna(predictions[i]) == True]
            test2 = [i for i in range(len(probabilities)) if pd.isna(probabilities[i]) == True]

            if len(test1) != len(test2):
                print("< dimension error2 occurred @mask_predictions ("+transcript_id+"): "+str(len(test1))+"/"+str(len(test2))+"/"+str(len(predictions))+"/"+str(len(probabilities)))
                for i in range(min(len(predictions), len(probabilities))):
                    print(i, predictions[i], probabilities[i])

                exit()
                        
        #return predictions, predictions_by_exon, probabilities, probabilities_by_exon
        return predictions, predictions_by_exon, probabilities, probabilities_by_exon, variant_targets