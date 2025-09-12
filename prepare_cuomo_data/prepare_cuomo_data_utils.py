import os
import pandas as pd
from progress.bar import IncrementalBar
# marked (<-) added on 250427 to allow integration of hg38 in combination with liftover
from pyliftover import LiftOver # <-
import sys
import time

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\prepare_data")

from prepare_data_utils import *


chrs = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4, "6": 5, "7": 6, "8": 7, "9": 8, "10": 9, "11": 10, "12": 11,
        "13": 12, "14": 13, "15": 14, "16": 15, "17": 16, "18": 17, "19": 18, "20": 19, "21": 20, "22": 21, "X": 22, "Y": 23}


genetic_code = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                "TGT": "C", "TGC": "C",
                "GAT": "D", "GAC": "D",
                "GAA": "E", "GAG": "E",
                "TTT": "F", "TTC": "F",
                "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
                "CAT": "H", "CAC": "H",
                "ATT": "I", "ATC": "I", "ATA": "I",
                "AAA": "K", "AAG": "K",
                "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                "ATG": "M",
                "AAT": "N", "AAC": "N",
                "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                "CAA": "Q", "CAG": "Q",
                "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
                "AGT": "S", "AGC": "S", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                "TGG": "W",
                "TAT": "Y", "TAC": "Y",
                "TAA": "X", "TGA": "X", "TAG": "X"}


class Prepare_cuomo_data_utils:
    def __init__(self, params):
        self.params = params

        if self.params["hg_build"] == "hg19":
            hg_build_dir    = "hg19"
            knowngene_fname = "hg19_knownGene.txt"

        # marked (<-) added on 250427 to allow integration of hg38 in combination with liftover
        if self.params["hg_build"] == "hg38": # <-
            hg_build_dir    = "hg38.p14" # <-
            knowngene_fname = "hg38_knownGene.txt" # <-


        self.genome = load_split_genome(self.params["data_dir1"]+self.params["os_sep"]+hg_build_dir, os_sep=self.params["os_sep"])
        print("< genome loaded.")

        # load knownGene database
        self.known_gene = pd.read_csv(params["data_dir1"]+params["os_sep"]+knowngene_fname, delimiter=",", index_col=False)
        
        # create dictionary
        # marked (<-) added / removed on 250427 to speed-up gene detection
        # self.known_gene_dict = create_knowngene_dict_by_position(self.known_gene, stepsize=params["stepsize"]) # <- removed
        self.known_gene_dict = create_knowngene_dict_by_position(self.known_gene, redundant=True, stepsize=params["stepsize"]) # <- added
        print("< dictionary created")

        # marked (<-) added on 250427 to allow integration of hg38 in combination with liftover
        self.liftover = LiftOver('hg19', 'hg38') # <- added


    def assign_mutants(self, mutants, ase, lock, fname, thread_id):
        # marked (<-) added on 250616 to initialize error report
        error_report = {key: {error: [] for error in [*self.params["error_handling"], "no_error"]} for key in ["total", "+total", "-total"]} # <- 

        bar = IncrementalBar("assigning mutants", max=ase.shape[0])
        for i in range(ase.shape[0]):
            if i % 1000 == 0:
                lock.acquire()
                try:
                    print("thread", thread_id, ": ", i, "of", ase.shape[0])

                finally:
                    lock.release()

            variant_id   = ase.iloc[i].iloc[0]
            chr          = "chr" + variant_id.split("_")[0]
            chr_index    = chrs[variant_id.split("_")[0]]
            pos          = int(variant_id.split("_")[1])

            # marked (<-) added on 250427 to facilitate the check of multiple possible errors prior to extraction
            passed = True # <-

            # marked (<-) added on 250427 to allow integration of hg38 in combination with liftover
            if self.params["hg_build"] == "hg38": # <-
                liftover_results = self.liftover.convert_coordinate(chr, pos) # <-
                
                if len(liftover_results) == 1: pos = liftover_results[0][1] # <- check for valid, unequivocal results
                else:                          passed = False # <-
            
            wt_base      = variant_id.split("_")[2]
            mutated_base = variant_id.split("_")[3]

            # marked (<-) added/removed on 250427 to facilitate the check of multiple possible errors prior to extraction
            #if len(wt_base) != 1 or len(mutated_base) != 1: # <- removed
            #    print("< no SNP @", i, "variant id:", variant_id) # <- removed

            # check whether assignment is correct # <- removed
            #elif self.genome[chr][pos-1] != wt_base: # <- removed
            #    print("< wrong assignment for variant", variant_id) # <- removed

            #else: # <- removed

            if len(wt_base) != 1 or len(mutated_base) != 1: # <- added
                print("< no SNP @", i, "variant id:", variant_id) # <- added
                passed = False # <- added

            if passed == True: # <- added
                converted_pos = int(pos/self.params["stepsize"])
                if converted_pos >= len(self.known_gene_dict[chr_index]):
                    print("< error occurred @thread", thread_id, "variant id:", variant_id, ". converted position exceeds dictionary dimension")

                else:
                    uc_ids = search_knowngene_dict_by_position(self.known_gene_dict, self.known_gene, chr_index, pos, stepsize=self.params["stepsize"])
                
                if len(uc_ids) > 0:
                    # extract possible cds sequences
                    known_gene_subset = self.known_gene[self.known_gene["uc id"].isin(uc_ids)]

                    for j in range(known_gene_subset.shape[0]):
                        info, fails   = get_info(known_gene_subset, variant_id, j, pos-1, pos-1, wt_base, None, self.genome, self.params, ptc_mode="absolute")
                        rel_ptc_index = info["FEATURE:ptc cds position"]
                        strand        = known_gene_subset.iloc[j].loc["strand"]
                        #print("i1", i, "j", j, self.known_gene.iloc[j].loc["uc id"])

                        #if "mutated_wt_base_not_found" in fails: # <- removed on 250616
                        #    print("< mutated wt base not found @", i, "variant id:", variant_id) # <- removed on 250616
                        
                        # marked (<-) added / removed on 250616
                        # if len(fails) == 0: # <- removed
                        if check_fails(fails, self.params) == False: # <- added
                            if variant_id == "": print("context1", self.genome[chr][pos-5:pos+5], self.genome[chr][pos-1:pos+2], "wt", wt_base)
                            
                            # use sequence information to classify mutation
                            ptc_index = self.check_for_ptc(info, mutated_base, wt_base, strand, variant_id)

                            if self.params["target_mutations"] == "ptc" and ptc_index != None:
                                # marked (<-) added / replaced on 250425 to calculate with shifted relative PTC index (+2) to avoid assignment to wrong exons of full stop codon
                                # because of changes to relative ptc index (+2), recalculation is required regardless
                                # added on 241105: re-evaluate if relative ptc position differs from previous one after frame-correction
                                # if ptc_index != info["FEATURE:ptc cds position"]: # <- replaced
                                #    info, fails = get_info(known_gene_subset, variant_id, j, ptc_index, rel_ptc_index, wt_base,
                                #                           None, self.genome, self.params, ptc_mode="relative") # <- replaced

                                info, fails = get_info(known_gene_subset, variant_id, j, ptc_index+2, rel_ptc_index, wt_base,
                                                       None, self.genome, self.params, ptc_mode="relative") # <- added

                                # marked (<-) added / removed on 250616 to fill error report
                                error_report = fill_error_report(error_report, pd.DataFrame(), fails, strand, i) # <- added

                                # if len(fails) == 0: # <- removed
                                if check_fails(fails, self.params) == False: # <- added
                                    for key, val in info.items():
                                        if "FEATURE" in key:
                                            mutants[key].append(val)

                                    mutants["ID:variant id"].append(variant_id)
                                    mutants["ID:gene id"].append(known_gene_subset.iloc[j].loc["gene id"])
                                    mutants["ID:transcript id"].append(known_gene_subset.iloc[j].loc["transcript id"])
                                    mutants["ID:gene symbol"].append(None)
                                    mutants["ID:uc id"].append(known_gene_subset.iloc[j].loc["uc id"])
                                    mutants["ID:cell id"].append(fname)
                                    mutants["ID:sample id"].append([ase.columns[k] for k in range(1, ase.shape[1], 1)])

                                    # marked (<-) added on 250425 to allowing checking for strand bias
                                    mutants["ID:strand"].append(strand) # <- added

                                    mutants["ID:appris annotation"].append(known_gene_subset.iloc[j].loc["appris annotation"])

                                    mutants["LABEL:alt counts"].append(None)
                                    mutants["LABEL:total counts"].append([ase.iloc[i].iloc[k] for k in range(1, ase.shape[1], 1)])


                            if self.params["target_mutations"] == "non-ptc" and ptc_index == None:
                                mutation = self.check_mutation(info, mutated_base, wt_base, strand, variant_id)

                                if len([ase.iloc[i].iloc[k] for k in range(1, ase.shape[1], 1) if ase.iloc[i].iloc[k] != None and float(ase.iloc[i].iloc[k]) >= 8]) > 0 and mutation[0] != mutation[1]:
                                    for key, val in info.items():
                                        if key in mutants and "FEATURE" in key:
                                            mutants[key].append(val)

                                    mutants["ID:variant id"].append(variant_id)
                                    mutants["ID:gene id"].append(known_gene_subset.iloc[j].loc["gene id"])
                                    mutants["ID:transcript id"].append(known_gene_subset.iloc[j].loc["transcript id"])
                                    mutants["ID:gene symbol"].append(None)
                                    mutants["ID:uc id"].append(known_gene_subset.iloc[j].loc["uc id"])
                                    mutants["ID:cell id"].append(fname)
                                    mutants["ID:sample id"].append([ase.columns[k] for k in range(1, ase.shape[1], 1)])

                                    mutants["ID:appris annotation"].append(known_gene_subset.iloc[j].loc["appris annotation"])

                                    mutants["LABEL:alt counts"].append(None)
                                    mutants["LABEL:total counts"].append([ase.iloc[i].iloc[k] for k in range(1, ase.shape[1], 1)])
                                
                                    mutants["ID:mutation type"].append(mutation)
        
            #bar.next()
        bar.finish()
        # marked (<-) added / removed on 250616
        # return mutants # <- removed
        return mutants, error_report # <- added


    def check_for_ptc(self, info, mutated_base, wt_base, strand, variant_id):
        cds = self.create_mutation(info, mutated_base, wt_base, strand, variant_id)

        # calculate in-frame position that is closest to ptc cds position
        cds_ptc_index = None
        i = int(info["FEATURE:ptc cds position"]/3)*3
        #while i < len(cds)-3 and cds_ptc_index == None:
        while i < int(info["FEATURE:ptc cds position"]/3)*3+3 and i < len(cds)-3 and cds_ptc_index == None:
            if cds[i:i+3] in ["TAA", "TAG", "TGA"]:
                cds_ptc_index = i

            i += 3
        
        return cds_ptc_index
    

    def check_mutation(self, info, mutated_base, wt_base, strand, variant_id):
        pos      = max(0, int(info["FEATURE:ptc cds position"]/3)*3)
        wt_codon = info["FEATURE:all cds"][pos:pos+3]
        mutation = genetic_code[wt_codon]

        cds           = self.create_mutation(info, mutated_base, wt_base, strand, variant_id)
        mutated_codon = cds[pos:pos+3]
        mutation     += genetic_code[mutated_codon]

        if genetic_code[wt_codon] != "X" and genetic_code[mutated_codon] == "X": print("< error occurred @check_mutation. mutation leads to stop codon @", variant_id)
        #print("mutation", mutation, "wt_codon", wt_codon, "mutated_codon", mutated_codon, mutated_base, wt_base, strand, variant_id)
        return mutation

    
    def create_mutation(self, info, mutated_base, wt_base, strand, variant_id):
        if strand == "+":
            test_wt_base = wt_base
            
        if strand == "-":
            test_wt_base = self.invert(wt_base)
            mutated_base = self.invert(mutated_base)

        if info["FEATURE:all cds"][info["FEATURE:ptc cds position"]] != test_wt_base:
            print("< error occurred @check for ptc @", variant_id, "wt base not detected.")
            print("wt base:", test_wt_base, "strand", strand, "sequence context:",
                  info["FEATURE:all cds"][max(0, info["FEATURE:ptc cds position"]-5):min(len(info["FEATURE:all cds"])-1, info["FEATURE:ptc cds position"]+5)],
                  "incorrect base:", info["FEATURE:all cds"][info["FEATURE:ptc cds position"]])
            
            print("second:", info["FEATURE:all cds"][info["FEATURE:ptc cds position"]:info["FEATURE:ptc cds position"]+2])

        # insert mutated base into sequence
        cds = ''.join(info["FEATURE:all cds"][i] if i != info["FEATURE:ptc cds position"] else mutated_base for i in range(len(info["FEATURE:all cds"])))

        if cds[info["FEATURE:ptc cds position"]] != mutated_base:
            print("< error occurred @check for ptc @", variant_id, "mutated base not detected.")
            print("wt base:", test_wt_base, "strand", strand, "sequence context:",
                  cds[max(0, info["FEATURE:ptc cds position"]-5):min(len(cds)-1, info["FEATURE:ptc cds position"]+5)],
                  "incorrect base:", cds[info["FEATURE:ptc cds position"]])
            
            print("second:", cds[info["FEATURE:ptc cds position"]:info["FEATURE:ptc cds position"]+2])

        return cds
    

    def invert(self, seq):
        inverse_seq = ""
        bases         = ["A", "C", "G", "T"]
        inverse_bases = ["T", "G", "C", "A"]
        for i in range(len(seq)):
            j = 0
            match = False
            while j < len(bases) and match == False:
                if seq[len(seq)-1-i] == bases[j]:
                    inverse_seq += inverse_bases[j]
                    match = True
                j += 1

        return inverse_seq
    

    def read_counts(self, mutants, ase):
        for i in range(mutants.shape[0]):
            variant_matches = ase.iloc[:,0] == mutants.iloc[i].loc["ID:variant id"]
            index           = ase.index[variant_matches == True].tolist()
            #print("i", i, index, mutants.iloc[i].loc["ID:variant id"], ase.iloc[index[0]].iloc[0])

            if len(index) == 0:
                print("< error occurred. no entry found for variant id:", mutants.iloc[i].loc["ID:variant id"])

            if len(index) > 1:
                print("< error occurred. multiple entries found for variant id:", mutants.iloc[i].loc["ID:variant id"])
                print(ase.iloc[index])

            if len(index) == 1:
                mutants.at[mutants.index[i], "LABEL:alt counts"] = [ase.iloc[index[0]].iloc[j] for j in range(1, ase.shape[1], 1)]

        return mutants