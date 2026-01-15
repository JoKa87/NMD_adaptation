from ast import literal_eval
from datetime import datetime
import copy # <- marked added on 250627
import json
import os
import sys
import statsmodels.api as sm # <- added on 250618

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir+"\\random_forest")
sys.path.insert(0, parent_dir+"\\shared")

from mask_predictions import * # <- added on 250618
from shared_utils import *


def append_df(df, dict):
    for key in dict:
        df[key] = dict[key]

    return df


def _calculate_inner_values(values, calculation, normalization_factor=None, threshold=None, show=False):
    inner_value = None
    if len(values) > 0:
        if calculation == "all_equal":
            counts = values.count(values[0])
            if counts == len(values): inner_value = values[0]

        if normalization_factor != None:
            values = [value/normalization_factor for value in values]

        if threshold != None:
            values = [value for value in values if value > threshold]

        if len(values) > 0:
            if calculation == "avg":    inner_value = np.mean(np.array(values))
            if calculation == "max":    inner_value = max(values)
            if calculation == "median": inner_value = np.median(np.array(values)) #print("v", values, len(values), np.mean(np.array(values)))
            if calculation == "sum":    inner_value = sum(values)
    
    return inner_value


def calculate_value(values, calculation, get_stats=False, normalization_factor=None, threshold=None, show=False):
    calculated_value = None
    calculated_stats = None

    values = json.loads(values)

    middle_values = []
    for i in range(len(values)):
        inner_values = []
        for j in range(len(values[i])):
            #try:
            inner_value = _calculate_inner_values(values[i][j], calculation[2], normalization_factor, threshold, show)
            if inner_value != None: inner_values.append(inner_value)

            #except:
            #    print("< exception occurred @calculate_value", values, calculation)

        middle_value = _calculate_inner_values(inner_values, calculation[1])
        if middle_value != None: middle_values.append(middle_value)

    #print("outer", middle_values)
    calculated_value = _calculate_inner_values(middle_values, calculation[0])

    if get_stats == True:
        if calculated_value != None and len(middle_values) > 1:
            calculated_stats  = np.var(np.array(middle_values)) # WARNING: np.var calculates variance like: var = sum(i,n) [ (value[i]-mean)^2 ] / n # NOT n-1!
            calculated_counts = len(middle_values)

        else:
            calculated_stats  = None
            calculated_counts = None

        return calculated_value, calculated_stats, calculated_counts
    
    else:
        return calculated_value


def load_firehouse_(data_dir, params):
    cnf       = None
    cnf_found = False
    if os.path.isdir(data_dir):
        fnames = os.listdir(data_dir)
        for fname in fnames:
            #print(fname)
            if ".bestclus" in fname:
                with open(data_dir+params["os_sep"]+fname, 'r') as _:
                    cnf       = pd.read_csv(data_dir+params["os_sep"]+fname, skiprows=2, delimiter="\t", names=["submitter_id", "cluster", "factor"], header=None)
                    cnf_found = True

        #print(data_dir, cnf)
    return cnf, cnf_found


def aggregate(cluster, clusters, params, axis=0, delimiter="\t", folder="calculations", path=None, redundant_cols=[]):
    aggregated_data = pd.DataFrame()

    #conversion_dict = {key: literal_eval for key in params["info"]["rna"] if "RNASEQ_" in key}
    bar = IncrementalBar("aggregating", max=len(clusters))
    for i in range(len(clusters)):
        project_key = clusters[i]["project_key"]
        cluster_key = clusters[i]["cluster_key"]
        project     = cluster[project_key][cluster_key]["project"]

        #if project in params["projects"] and os.path.isdir(params["data_dir"]+params["os_sep"]+project) == True:
        if (project in params["projects"] or params["projects"] == "all") and os.path.isdir(params["data_dir"]+params["os_sep"]+project) == True:
            sub_dir = params["data_dir"]+params["os_sep"]+folder

            try:
                if os.path.isfile(sub_dir+params["os_sep"]+project_key+"_"+cluster_key+".txt") == True:
                    data = pd.read_csv(sub_dir+params["os_sep"]+project_key+"_"+cluster_key+".txt", delimiter=delimiter)
                    
                    if axis == 0:
                        if aggregated_data.shape[0] > 0: aggregated_data = pd.concat([aggregated_data, data], ignore_index=True)
                        else:                            aggregated_data = data

                    if axis == 1:
                        if aggregated_data.shape[0] > 0: data = data.drop(columns=redundant_cols); aggregated_data = pd.concat([aggregated_data, data], axis=1)
                        else:                            aggregated_data = data

                    del data
                
                else:
                    print("< loading of", sub_dir+params["os_sep"]+project_key+"_"+cluster_key+".txt", "was not successfull.")

                
            except Exception as error:
                print("< loading of", sub_dir+params["os_sep"]+project_key+"_"+cluster_key+".txt", "was not successfull.")
                print(" ", error)

        bar.next()
    bar.finish()

    if path == None: aggregated_data.to_csv(params["run_dir"]+params["os_sep"]+params["fname"], index=False, sep=delimiter)
    else:            aggregated_data.to_csv(path, index=False, sep=delimiter)
    return cluster


def calculate_expression(df, params):
    bar = IncrementalBar("expression is calculated", max=df.shape[0])

    # initialize containers to store calculation results
    rna_counts  = {rna_value + "_counts": [None for _ in range(df.shape[0])] for rna_value in params["rna_values"]}
    rna_stats   = {rna_value + "_stats":  [None for _ in range(df.shape[0])] for rna_value in params["rna_values"]}
    rna_values  = {rna_value:             [None for _ in range(df.shape[0])] for rna_value in params["rna_values"]}
    expression  = {"NMD score":           [None for _ in range(df.shape[0])],
                   "alt. NMD score":      [None for _ in range(df.shape[0])]}
    expressions = []

    # type conversion (currently, all values must be submitted in str format)
    for rna_value in params["rna_values"]:
        df[rna_value] = df[rna_value].astype(str)

    for i in range(df.shape[0]):
        if i % 1000 == 0: print("i", i, datetime.now(), df.iloc[i].loc["project"], "/", df.shape[0])
        #bar.next()

        # calculate rna values based on defined approach
        # outer dimension: expression data from different cases
        # medium dimension: expression data from different rnaseq runs of a single case (sometimes occurring)
        # inner dimension: different expression data from single rnaseq run due to overlap (unclear, whether occurring at all)
        ptc_expression = None; ref_expression = None
        for rna_value in params["rna_values"]:
            #if type(df.iloc[i].loc[rna_value]) == str or (df.iloc[i].loc[rna_value] != str and math.isnan(df.iloc[i].loc[rna_value]) == False):
            if type(df.iloc[i].loc[rna_value]) == str or (df.iloc[i].loc[rna_value] != str and pd.isna(df.iloc[i].loc[rna_value]) == False):
                # apply filter to exclude values below threshold based on category

                if "cnv" in rna_value:
                    rna_values[rna_value][i], rna_stats[rna_value + "_stats"][i], rna_counts[rna_value + "_counts"][i] = calculate_value(df.iloc[i].loc[rna_value], params["rna_calculation"],
                                                                                                                                         get_stats=True, threshold=params["cnv_threshold"])
                
                else:
                    rna_values[rna_value][i], rna_stats[rna_value + "_stats"][i], rna_counts[rna_value + "_counts"][i] = calculate_value(df.iloc[i].loc[rna_value], params["rna_calculation"],
                                                                                                                                         get_stats=True, threshold=params["rna_threshold"])
                
                #if params["rna_criterion"] in rna_value and rna_values[rna_value][i] != None:
                if params["rna_criterion"] in rna_value and pd.isna(rna_values[rna_value][i]) == False:
                    if "_noptc_" in rna_value: ref_expression = rna_values[rna_value][i]
                    elif "_ptc_" in rna_value: ptc_expression = rna_values[rna_value][i]
        
        if pd.isna(ptc_expression) == False and pd.isna(ref_expression) == False and (ref_expression + ptc_expression) > 0:
            if ptc_expression > 0 and ref_expression > 0: expression["NMD score"][i]      = -math.log2(ptc_expression/ref_expression)
            if ref_expression > 0:                        expression["alt. NMD score"][i] = (ref_expression + ptc_expression) / (2*ref_expression)

            if pd.isna(expression["NMD score"][i]) == False:
                expressions.append(expression["NMD score"][i])
                
    bar.finish()

    df = append_df(df, rna_counts)
    df = append_df(df, rna_stats)
    df = append_df(df, rna_values)
    df = append_df(df, expression)
    return df


def check_variants(wxs, ptcs, params, it):
    is_ptc = True

    # create subset of all entries with the same gene id
    # changed 241117
    #gene_variants = wxs[wxs["Gene"] == ptcs.iloc[it].loc["Gene"]] 
    gene_variants = wxs[wxs[params["wxs_identifier"]] == ptcs.iloc[it].loc[params["wxs_identifier"]]]

    if gene_variants.shape[0] > 1:    
        # ptc position of the selected variant
        ref_ptc_pos = get_ptc(ptcs.iloc[it].loc["HGVSp"])
        i = 0
        while i < gene_variants.shape[0] and is_ptc == True:
            if type(gene_variants.iloc[i].loc["HGVSp"]) == str and "ext" not in gene_variants.iloc[i].loc["HGVSp"] and "?" not in gene_variants.iloc[i].loc["HGVSp"]:
                # check whether protein position is downstream of reference (otherwise not relevant), changed 241107
                if "Ter" in gene_variants.iloc[i].loc["HGVSp"]:
                    ptc_pos = get_ptc(gene_variants.iloc[i].loc["HGVSp"])

                else:
                    ptc_pos = None
                    positions = get_numbers(gene_variants.iloc[i].loc["HGVSp"])
                    if len(positions) > 0: ptc_pos = positions[0]
                
                if ptc_pos != None and ref_ptc_pos != None and ptc_pos < ref_ptc_pos:
                    # check whether mutation contains insertions, deletions, duplications, and/or ptc (otherwise not relevant)
                    if("del" in gene_variants.iloc[i].loc["HGVSp"]
                        or "dup" in gene_variants.iloc[i].loc["HGVSp"]
                        or "ins" in gene_variants.iloc[i].loc["HGVSp"]
                        or "Ter" in gene_variants.iloc[i].loc["HGVSp"]):
                        is_ptc = False
            i += 1

    return is_ptc


# marked function (<-) added on 250627
def _create_mutation_stats(mutation_stats, category, project, seq):
    mutation_stats[category][project][seq]     += 1
    mutation_stats[category][project]["total"] += 1
    mutation_stats[category]["total"][seq]     += 1
    mutation_stats[category]["total"]["total"] += 1
    return mutation_stats
    

# marked function (<-) added on 250603
def create_mutation_stats(cluster, clusters, seqs, params):
    exon_mutations = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                      "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Site"] #  
    bases          = ["A", "C", "G", "T"]


    # added on 251027 to allow alternative project keys (required for CPTAC-3)
    if params["external_project_path"] != None:
        clinical_data       = pd.read_csv(params["external_project_path"], sep="\t")
        clinical_data.index = clinical_data["cases.case_id"]

        for i, case_id in enumerate(clinical_data["cases.case_id"]):
            if clinical_data[clinical_data["cases.case_id"] == case_id].drop_duplicates(subset="cases.disease_type").shape[0] > 1:
                print("< redundant case with deviating disease types")
                print(i, case_id, clinical_data[clinical_data["cases.case_id"] == case_id]["cases.disease_type"].tolist())
        
        init_shape          = clinical_data.shape[0]
        clinical_data       = clinical_data.drop_duplicates(subset=["cases.case_id"])
        print("< removal of duplicate cases reduced clinical dataset from", init_shape, "to", clinical_data.shape[0])
        project_keys        = [*np.unique(clinical_data["cases.disease_type"]).tolist(), "total"]

    else:
        project_keys        = [*list(cluster.keys()), "total"]


    # prepare containers
    dubletts = []; tripletts = []; quadrupletts = []
    for base1 in bases:
        for base2 in bases:
            dubletts.append(base1+base2)

            for base3 in bases:
                tripletts.append(base1+base2+base3)

                for base4 in bases:
                    quadrupletts.append(base1+base2+base3+base4)
    
    '''
    mutation_stats = {
                      **{"2mers":       {dublett: 0 for dublett in [*dubletts, "total"]}},
                      **{"3mers":       {triplett: 0 for triplett in [*tripletts, "total"]}},
                      **{"bases":       {base: 0 for base in [*bases, "total"]}},
                      **{"cases":       {key: {} for key in [*list(cluster.keys()), "total"]}},
                      **{"del-1":       {key: {base: 0 for base in [*bases, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"del-1_3mers": {key: {triplett: 0 for triplett in [*tripletts, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"del-1_full":  {key: 0 for key in [*list(cluster.keys()), "total"]}},
                      **{"genes":       {key: {} for key in [*list(cluster.keys()), "total"]}},
                      **{"ins+1":       {key: {base: 0 for base in [*bases, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"ins+1_3mers": {key: {triplett: 0 for triplett in [*tripletts, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"ins+1_full":  {key: 0 for key in [*list(cluster.keys()), "total"]}},
                      **{"observed":    {key: {**{mutation: {} for mutation in ["del-1", "frameshift", "ins+1", "missense", "nonsense"]}}
                                         for key in [*list(cluster.keys()), "total"]}},
                      **{"pairs":       {key: {dublett: 0 for dublett in [*dubletts, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"pairs_3mers": {key: {quadruplett: 0 for quadruplett in [*quadrupletts, "total"]} for key in [*list(cluster.keys()), "total"]}},
                      **{"seq_size":    0}
                     }

    report = {**{project: 0 for project in cluster}, **{"exceeding_positions": 0, "mismatching_frameshift_positions": 0, "mismatching_positions": 0,
                                                        "mismatching_positions+": 0, "mismatching_positions-": 0,
                                                        "sequence_deviations": 0, "total_mutations": 0, 
                                                        "exceeding_transcripts": [], "missing_transcripts": []}}
    '''

    mutation_stats = {
                      **{"2mers":       {dublett: 0 for dublett in [*dubletts, "total"]}},
                      **{"3mers":       {triplett: 0 for triplett in [*tripletts, "total"]}},
                      **{"bases":       {base: 0 for base in [*bases, "total"]}},
                      **{"cases":       {key: {} for key in project_keys}},
                      **{"del-1":       {key: {base: 0 for base in [*bases, "total"]} for key in project_keys}},
                      **{"del-1_3mers": {key: {triplett: 0 for triplett in [*tripletts, "total"]} for key in project_keys}},
                      **{"del-1_full":  {key: 0 for key in project_keys}},
                      **{"genes":       {key: {} for key in project_keys}},
                      **{"ins+1":       {key: {base: 0 for base in [*bases, "total"]} for key in project_keys}},
                      **{"ins+1_3mers": {key: {triplett: 0 for triplett in [*tripletts, "total"]} for key in project_keys}},
                      **{"ins+1_full":  {key: 0 for key in project_keys}},
                      **{"observed":    {key: {**{mutation: {} for mutation in ["del-1", "frameshift", "ins+1", "missense", "nonsense"]}}
                                         for key in project_keys}},
                      **{"pairs":       {key: {dublett: 0 for dublett in [*dubletts, "total"]} for key in project_keys}},
                      **{"pairs_3mers": {key: {quadruplett: 0 for quadruplett in [*quadrupletts, "total"]} for key in project_keys}},
                      **{"seq_size":    0}
                     }   
    
    report = {**{project: 0 for project in project_keys}, **{"exceeding_positions": 0, "mismatching_frameshift_positions": 0, "mismatching_positions": 0,
                                                             "mismatching_positions+": 0, "mismatching_positions-": 0,
                                                             "sequence_deviations": 0, "total_mutations": 0, 
                                                             "exceeding_transcripts": [], "missing_transcripts": []}}

    mutation_types = ["Silent", "Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "interval"]
    trajectory     = pd.DataFrame({mutation_type: [0 for _ in range(100)] for mutation_type in mutation_types})    

    bar = IncrementalBar(set_bar("creating mutation stats"), max=len(clusters))
    for project_key in cluster:
        if params["projects"] == "all" or project_key in params["projects"]:
            for cluster_key in cluster[project_key]:
                for i in range(len(cluster[project_key][cluster_key]["input"])):
                    for case_id in cluster[project_key][cluster_key]["input"][i]:
                        wxss = []
                        # added on 251027 to allow alternative project keys (required for CPTAC-3)
                        if params["external_project_path"] != None:
                            external_project_key = clinical_data.loc[case_id].loc["cases.disease_type"]

                        else:
                            external_project_key = project_key

                        # load WXS files
                        for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"])):
                            wxs = load_wxs(params["data_dir"]+params["os_sep"]+project_key+params["os_sep"]+"WXS"
                                           +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][j])

                            report[external_project_key] += 1 # report[project_key] += 1
                            if wxs.shape[0] > 0: wxss.append(wxs)

                        # merge WXS files and apply filters
                        if len(wxss) > 0:
                            # remove redundant variants
                            wxs = pd.concat(wxss).drop_duplicates(subset=["Transcript_ID", "HGVSc"])

                            # select specified mutations
                            wxs = wxs[wxs["Variant_Classification"].isin(exon_mutations)]

                            # remove irregular entries
                            wxs = wxs[~wxs["CDS_position"].isna()]
                            wxs = wxs[[False if "?" in cds_pos or "/" not in cds_pos or len(cds_pos.split("/")[0]) == 0 or len(cds_pos.split("/")[0].split("-")[0]) == 0 else True
                                       for cds_pos in wxs["CDS_position"]]]
                            
                            # filter out transcript ids with no sequence available
                            if "Transcript_ID" in wxs.columns:
                                wxs = wxs[wxs["Transcript_ID"].isin(seqs.index)]

                                if wxs[~wxs["Transcript_ID"].isin(seqs.index)].shape[0] > 0:
                                    report["Transcript_ID"].extend(wxs[~wxs["Transcript_ID"].isin(seqs.index)]["Transcript_ID"].tolist())

                            # test to exclude that the same variant is assigned multiple times to different transcripts of the same gene
                            if wxs.drop_duplicates(subset=["Gene", "Start_Position", "End_Position"]).shape[0] != wxs.shape[0]:
                                print("< redundant genes found for project", project_key, "WXS", cluster[project_key][cluster_key]["input"][i][case_id]["WXS"])
                                exit()

                        else:
                            wxs = pd.DataFrame()
                        
                        for j in range(wxs.shape[0]):
                            # register mutations
                            ref_base = wxs.iloc[j].loc["Reference_Allele"]	
                            alt_base = wxs.iloc[j].loc["Tumor_Seq_Allele2"]
                            cds      = seqs.loc[wxs.iloc[j].loc["Transcript_ID"]].loc["cds"]
                            strand   = seqs.loc[wxs.iloc[j].loc["Transcript_ID"]].loc["strand"]
                            seq_size = float(wxs.iloc[j].loc["CDS_position"].split("/")[1])

                            # determine sequence position from CDS_position
                            if "-" not in wxs.iloc[j].loc["CDS_position"].split("/")[0]:
                                seq_pos = int(wxs.iloc[j].loc["CDS_position"].split("/")[0])-1

                            else:
                                seq_pos = wxs.iloc[j].loc["CDS_position"].split("/")[0]
                                seq_pos = int(seq_pos.split("-")[0])-1

                            # filter entries (filters_passed1 for context-resolved statistics, filters_passed2 for whole-gene statistics)
                            filters_passed1 = True; filters_passed2 = True

                            # filter for transcripts exceeding CDS size (applies to sequence-specific statistics AND general statistics)
                            if seq_pos+1 >= len(cds):
                                filters_passed1 = False; filters_passed2 = False
                                report["exceeding_transcripts"].extend(wxs.iloc[j].loc["Transcript_ID"])

                            # filter to exclude mutations located in start or stop codon (applies to sequence-specific statistics AND general statistics)
                            if seq_pos < 3 or seq_pos >= seq_size-3:
                                filters_passed1 = False; filters_passed2 = False
                                report["exceeding_positions"] += 1
                            
                            # register trajectory
                            if seq_pos >= 3 and seq_pos < seq_size-3:
                                if wxs.iloc[j].loc["Variant_Classification"] in trajectory.columns:
                                    pos = int(100*seq_pos/(seq_size-3))
                                    trajectory.at[pos, wxs.iloc[j].loc["Variant_Classification"]] += 1
                                    trajectory.at[pos, "interval"]                                += ((pos+1)*(seq_size-3)/100)-((pos)*(seq_size-3)/100)

                            # filters for mutations of size 1 only (applies to sequence-specific statistics AND general statistics)
                            if len(ref_base) != 1 or len(alt_base) != 1:
                                filters_passed1 = False; filters_passed2 = False

                            # CDS_position is identical with Start_Position, HGVSc position is in agreement with HGVSp position
                            # indicated alt_base is in agreement with HGVSc position
                            # due to these ambiguities, it is probably necessary to exclude variants with ambiguous HGVSc (applies to sequence-specific statistics ONLY)                            
                            if len(get_numbers(wxs.iloc[j].loc["HGVSc"])) > 0 and get_numbers(wxs.iloc[j].loc["HGVSc"])[0]-1 != seq_pos:
                                filters_passed1 = False
                                report["mismatching_positions"]        += 1
                                report["mismatching_positions"+strand] += 1

                                if wxs.iloc[j].loc["Variant_Classification"] in ["Frame_Shift_Del", "Frame_Shift_Ins"]:
                                    report["mismatching_frameshift_positions"] += 1

                            # filter to exclude nonsense mutation in order to reduce possible bias (applies to sequence-specific statistics ONLY)
                            if wxs.iloc[j].loc["Variant_Classification"] in ["Nonsense_Mutation", "Nonstop_Mutation"]:
                                filters_passed1 = False

                            if filters_passed1 == True:
                                report["total_mutations"] += 1

                                if strand == "+":
                                    pre_base  = cds[seq_pos-1]
                                    base      = cds[seq_pos] # required for insertions
                                    post_base = cds[seq_pos+1]

                                    if cds[seq_pos] != ref_base:
                                        report["sequence_deviations"] += 1

                                if strand == "-":
                                    # test for effect of using non-inverted sequence (with inverted ref_base and/or alt_base) showed no differences for deletions or SNPs (pairs)
                                    # differences were observed for insertions due to deviating order as post_bases are not identical (only for total mutation counts)
                                    # new implementation for indels is approved for context-specific statistics
                                    pre_base  = invert_sequence(cds[seq_pos+1])
                                    base      = invert_sequence(cds[seq_pos]) # required for insertions
                                    post_base = invert_sequence(cds[seq_pos-1])
                                    '''
                                    complementary approach:
                                    # mutation_stats calculated with above code (irrelevant for full probabilities, different for triplett-resolved probabilities)
                                    pre_base  = cds[seq_pos-1]
                                    base      = cds[seq_pos] # required for insertions
                                    post_base = cds[seq_pos+1]
                                    if ref_base != "-": ref_base = invert_sequence(ref_base)
                                    if alt_base != "-": alt_base = invert_sequence(alt_base)
                                    '''

                                    # register deviations of reference base (possible due to SNPs, but should not be too many)
                                    if cds[seq_pos] != invert_sequence(ref_base):
                                        report["sequence_deviations"] += 1

                                    if invert_sequence(pre_base+base+post_base) != cds[seq_pos-1:seq_pos+2]:
                                        print("< inversion error occurred:", invert_sequence(pre_base+base+post_base), "/", cds[seq_pos-1:seq_pos+2])
                                        exit()
                                
                                '''
                                # SNP mutations
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base in bases and alt_base in bases:
                                    mutation_stats = _create_mutation_stats(mutation_stats, "pairs", project_key, ref_base+alt_base)
                                    mutation_stats = _create_mutation_stats(mutation_stats, "pairs_3mers", project_key, pre_base+ref_base+alt_base+post_base)

                                # deletion of size 1
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base in bases and alt_base == "-":
                                    mutation_stats = _create_mutation_stats(mutation_stats, "del-1", project_key, ref_base)
                                    mutation_stats = _create_mutation_stats(mutation_stats, "del-1_3mers", project_key, pre_base+ref_base+post_base)

                                # insertion of size 1
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base == "-" and alt_base in bases:
                                    mutation_stats = _create_mutation_stats(mutation_stats, "ins+1", project_key, alt_base)
                                    if strand == "+": mutation_stats = _create_mutation_stats(mutation_stats, "ins+1_3mers", project_key, base+alt_base+post_base)
                                    if strand == "-": mutation_stats = _create_mutation_stats(mutation_stats, "ins+1_3mers", project_key, pre_base+alt_base+base)
                                '''

                                # SNP mutations
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base in bases and alt_base in bases:
                                    mutation_stats = _create_mutation_stats(mutation_stats, "pairs", external_project_key, ref_base+alt_base)
                                    mutation_stats = _create_mutation_stats(mutation_stats, "pairs_3mers", external_project_key, pre_base+ref_base+alt_base+post_base)

                                # deletion of size 1
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base in bases and alt_base == "-":
                                    mutation_stats = _create_mutation_stats(mutation_stats, "del-1", external_project_key, ref_base)
                                    mutation_stats = _create_mutation_stats(mutation_stats, "del-1_3mers", external_project_key, pre_base+ref_base+post_base)

                                # insertion of size 1
                                if len(ref_base) == 1 and len(alt_base) == 1 and ref_base == "-" and alt_base in bases:
                                    mutation_stats = _create_mutation_stats(mutation_stats, "ins+1", external_project_key, alt_base)
                                    if strand == "+": mutation_stats = _create_mutation_stats(mutation_stats, "ins+1_3mers", external_project_key, base+alt_base+post_base)
                                    if strand == "-": mutation_stats = _create_mutation_stats(mutation_stats, "ins+1_3mers", external_project_key, pre_base+alt_base+base)
                                    
                                    '''
                                    complementary approach:
                                    mutation_stats = _create_mutation_stats(mutation_stats, "ins+1_3mers", project_key, base+alt_base+post_base)
                                    '''

                            # create gene-specific statistics
                            # per-gene mutation counts
                            for key in ["total", external_project_key]: # for key in ["total", project_key]:
                                if wxs.iloc[j].loc["Transcript_ID"] in mutation_stats["genes"][key]: 
                                    mutation_stats["genes"][key][wxs.iloc[j].loc["Transcript_ID"]]                  += 1 / seq_size

                                else:
                                    mutation_stats["genes"][key][wxs.iloc[j].loc["Transcript_ID"]]                   = 1 / seq_size
                                    mutation_stats["observed"][key]["nonsense"][wxs.iloc[j].loc["Transcript_ID"]]    = 0
                                    mutation_stats["observed"][key]["del-1"][wxs.iloc[j].loc["Transcript_ID"]]       = 0
                                    mutation_stats["observed"][key]["frameshift"][wxs.iloc[j].loc["Transcript_ID"]]  = 0
                                    mutation_stats["observed"][key]["ins+1"][wxs.iloc[j].loc["Transcript_ID"]]       = 0
                                    mutation_stats["observed"][key]["missense"][wxs.iloc[j].loc["Transcript_ID"]]    = 0
                            
                                # per-gene observed mutations
                                if filters_passed2 == True:
                                    if wxs.iloc[j].loc["Variant_Classification"] == "Nonsense_Mutation":
                                        mutation_stats["observed"][key]["nonsense"][wxs.iloc[j].loc["Transcript_ID"]]   += 1

                                    elif wxs.iloc[j].loc["Variant_Classification"] in ["Frame_Shift_Del", "Frame_Shift_Ins"]:
                                        mutation_stats["observed"][key]["frameshift"][wxs.iloc[j].loc["Transcript_ID"]] += 1

                                        if wxs.iloc[j].loc["Variant_Classification"]  == "Frame_Shift_Del":
                                            mutation_stats["observed"][key]["del-1"][wxs.iloc[j].loc["Transcript_ID"]]  += 1

                                        if wxs.iloc[j].loc["Variant_Classification"]  == "Frame_Shift_Ins":
                                            mutation_stats["observed"][key]["ins+1"][wxs.iloc[j].loc["Transcript_ID"]]  += 1

                                    elif wxs.iloc[j].loc["Variant_Classification"] != "Silent":
                                        mutation_stats["observed"][key]["missense"][wxs.iloc[j].loc["Transcript_ID"]]   += 1

                        if len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"]) > 0:
                            mutation_stats["cases"]["total"][case_id]              = wxs.shape[0]
                            mutation_stats["cases"][external_project_key][case_id] = wxs.shape[0] # mutation_stats["cases"][project_key][case_id] = wxs.shape[0]

                bar.next()
    bar.finish()

    # print trajectory
    trajectory.to_csv(params["run_dir"]+params["os_sep"]+params["fname"].split(".")[0]+"_trajectory.txt", index=False, sep=",")


    # print report
    for key in report:
        if type(report[key]) == list:
            report[key] = np.unique(report[key]).shape[0]

    print(json.dumps(report, indent=4))

    # code check phase 4 conducted until here (250709)
    # print observed mutations
    for key1 in mutation_stats["observed"]["total"]:
        print("<", key1, np.sum([mutation_stats["observed"]["total"][key1][key2] for key2 in mutation_stats["observed"]["total"][key1]]))


    # test whether total counts can be recalculated
    for key1 in ["cases", "genes"]:
        test_total = 0
        for key2 in mutation_stats[key1]:
            test_total += np.sum([mutation_stats[key1][key2][key3] for key3 in mutation_stats[key1][key2] if key2 != "total"])
        
        if test_total != np.sum([mutation_stats[key1]["total"][key2] for key2 in mutation_stats[key1]["total"]]):
            print("< mismatching counts1 @"+key1+"/"+key2+" ("+str(test_total)+"/"+str(np.sum([mutation_stats[key1]["total"][key2] for key2 in mutation_stats[key1]["total"]]))+")")


    for key1 in ["del-1", "del-1_3mers", "ins+1", "ins+1_3mers", "pairs", "pairs_3mers"]:
        test_total1 = {key2: 0 for key2 in mutation_stats[key1][list(mutation_stats[key1].keys())[0]]}
        for key2 in mutation_stats[key1]:
            if key2 != "total":
                for key3 in mutation_stats[key1][key2]:
                    test_total1[key3] += mutation_stats[key1][key2][key3]

            test_total2 = 0
            for key3 in mutation_stats[key1][key2]:
                if key3 != "total":
                    test_total2 += mutation_stats[key1][key2][key3]

            if test_total2 != mutation_stats[key1][key2]["total"]:
                print("< mismatching counts3 @"+key1+"/"+key2+" ("+str(test_total2)+"/"+str(mutation_stats[key1][key2]["total"])+")")

        for key2 in test_total1:
            if test_total1[key2] != mutation_stats[key1]["total"][key2]:
                print("< mismatching counts2 @"+key1+"/"+key2+" ("+str(test_total1[key2])+"/"+str(mutation_stats[key1]["total"][key2])+")")

    # get base statistics from full sequences and aggregate sequence sizes
    bar = IncrementalBar(set_bar("retrieving base statistics"), max=len(mutation_stats["genes"]["total"]))
    for transcript_id in mutation_stats["genes"]["total"]:
        if seqs.loc[transcript_id].loc["strand"] == "+":
            cds = seqs.loc[transcript_id].loc["cds"][2:len(seqs.loc[transcript_id].loc["cds"])-2]

        if seqs.loc[transcript_id].loc["strand"] == "-":
            cds = invert_sequence(seqs.loc[transcript_id].loc["cds"][2:len(seqs.loc[transcript_id].loc["cds"])-2])

        mutation_stats["seq_size"] += len(cds)-2 # only for report purpose

        total = 0
        for base in bases:
            #print(base, seqs.loc[gene].loc["cds"][0:len(seqs.loc[gene].loc["cds"])-3].count(base) )
            mutation_stats["bases"][base] += sum(cds[i] == base for i in range(1, len(cds)-1))
            total                         += sum(cds[i] == base for i in range(1, len(cds)-1))

        if total != len(cds)-2:
            print("< inconsistent dublett count ("+str(total)+"/"+str(len(cds)-2)+")")

        mutation_stats["bases"]["total"] += total

        total = 0
        for dublett in mutation_stats["2mers"]:
            mutation_stats["2mers"][dublett] += sum(cds[i:i+2] == dublett for i in range(1, len(cds)-1))
            total                            += sum(cds[i:i+2] == dublett for i in range(1, len(cds)-1))

        if total != len(cds)-2:
            print("< inconsistent dublett count ("+str(total)+"/"+str(len(cds)-2)+")")

        mutation_stats["2mers"]["total"] += total

        total = 0
        for triplett in mutation_stats["3mers"]:
            mutation_stats["3mers"][triplett] += sum(cds[i:i+3] == triplett for i in range(len(cds)-2))
            total                             += sum(cds[i:i+3] == triplett for i in range(len(cds)-2))

        if total != len(cds)-2:
            print("< inconsistent triplett count ("+str(total)+"/"+str(len(cds)-2)+")")

        mutation_stats["3mers"]["total"] += total
        bar.next()
    bar.finish()

    
    # average complementary bases
    if params["calculate_complements"] == True:
        bases = copy.deepcopy(mutation_stats["bases"]); total_2mers = copy.deepcopy(mutation_stats["2mers"]); total_3mers = copy.deepcopy(mutation_stats["3mers"])
        for base in mutation_stats["bases"]:
            if base != "total":
                inverted_base = invert_sequence(base)
                bases[base]   = (mutation_stats["bases"][base]+mutation_stats["bases"][inverted_base])/2

        for dublett in mutation_stats["2mers"]:
            if dublett != "total":
                inverted_dublett     = invert_sequence(dublett)
                total_2mers[dublett] = (mutation_stats["2mers"][dublett]+mutation_stats["2mers"][inverted_dublett])/2

        for triplett in mutation_stats["3mers"]:
            if triplett != "total":
                inverted_triplett     = invert_sequence(triplett)
                total_3mers[triplett] = (mutation_stats["3mers"][triplett]+mutation_stats["3mers"][inverted_triplett])/2

        mutation_stats["bases"] = bases
        mutation_stats["2mers"] = total_2mers
        mutation_stats["3mers"] = total_3mers

        for key in mutation_stats["del-1"]:
            # calculate complements for del-1 and ins+1
            del_1 = copy.deepcopy(mutation_stats["del-1"][key]); ins_1 = copy.deepcopy(mutation_stats["ins+1"][key])
            for base in mutation_stats["del-1"][key]:
                if base != "total":
                    inverted_base = invert_sequence(base)
                    del_1[base]   = (mutation_stats["del-1"][key][base]+mutation_stats["del-1"][key][inverted_base])/2
                    ins_1[base]   = (mutation_stats["ins+1"][key][base]+mutation_stats["ins+1"][key][inverted_base])/2
        
            mutation_stats["del-1"][key] = del_1
            mutation_stats["ins+1"][key] = ins_1

            # calculate complements for del-1_3mers and ins+1_3mers
            del_1_3mers = copy.deepcopy(mutation_stats["del-1_3mers"][key]); ins_1_3mers = copy.deepcopy(mutation_stats["ins+1_3mers"][key])
            for triplett in mutation_stats["del-1_3mers"][key]:
                if triplett != "total":
                    inverted_triplett     = invert_sequence(triplett)
                    del_1_3mers[triplett] = (mutation_stats["del-1_3mers"][key][triplett]+mutation_stats["del-1_3mers"][key][inverted_triplett])/2
                    ins_1_3mers[triplett] = (mutation_stats["ins+1_3mers"][key][triplett]+mutation_stats["ins+1_3mers"][key][inverted_triplett])/2
        
            mutation_stats["del-1_3mers"][key] = del_1_3mers
            mutation_stats["ins+1_3mers"][key] = ins_1_3mers

            # average complementary pairs
            pairs = copy.deepcopy(mutation_stats["pairs"][key])
            for pair in mutation_stats["pairs"][key]:
                if pair != "total":
                    inverted_pair = invert_sequence(pair[0])+invert_sequence(pair[1])
                    pairs[pair]   = (mutation_stats["pairs"][key][pair]+mutation_stats["pairs"][key][inverted_pair])/2

            mutation_stats["pairs"][key] = pairs

            # average complementary pairs_3mers
            pairs_3mers = copy.deepcopy(mutation_stats["pairs_3mers"][key])
            for quadruplett in mutation_stats["pairs_3mers"][key]:
                if quadruplett != "total":
                    inverted_triplett        = invert_sequence(quadruplett[0]+quadruplett[1]+quadruplett[3])
                    inverted_quadruplett     = inverted_triplett[0]+inverted_triplett[1]+invert_sequence(quadruplett[2])+inverted_triplett[2]
                    pairs_3mers[quadruplett] = (mutation_stats["pairs_3mers"][key][quadruplett]+mutation_stats["pairs_3mers"][key][inverted_quadruplett])/2
                   
            mutation_stats["pairs_3mers"][key] = pairs_3mers

    # normalize bases/pairs to base counts
    for key in mutation_stats["pairs"]:
        for pair in mutation_stats["pairs"][key]:
            if pair != "total" and mutation_stats["bases"][pair[0]]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["pairs"][key][pair] /= (mutation_stats["bases"][pair[0]]*len(mutation_stats["cases"][key]))

            elif mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["pairs"][key]["total"] /= (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]))

        for quadruplett in mutation_stats["pairs_3mers"][key]:
            if quadruplett != "total" and mutation_stats["3mers"][quadruplett[0]+quadruplett[1]+quadruplett[3]]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["pairs_3mers"][key][quadruplett] /= (mutation_stats["3mers"][quadruplett[0]+quadruplett[1]+quadruplett[3]]*len(mutation_stats["cases"][key]))

            elif mutation_stats["3mers"]["total"]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["pairs_3mers"][key]["total"] /= (mutation_stats["3mers"]["total"]*len(mutation_stats["cases"][key]))

        for base in mutation_stats["del-1"][key]:
            if base != "total" and mutation_stats["bases"][base]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["del-1"][key][base] /= (mutation_stats["bases"][base]*len(mutation_stats["cases"][key]))
                mutation_stats["ins+1"][key][base] /= (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]))

            elif mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["del-1"][key]["total"] /= (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]))
                mutation_stats["ins+1"][key]["total"] /= (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key]))

        for triplett in mutation_stats["del-1_3mers"][key]:
            if triplett != "total" and mutation_stats["3mers"][triplett]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["del-1_3mers"][key][triplett] /= (mutation_stats["3mers"][triplett]*len(mutation_stats["cases"][key]))
                mutation_stats["ins+1_3mers"][key][triplett] /= (mutation_stats["2mers"][triplett[0]+triplett[2]]*len(mutation_stats["cases"][key]))

            elif mutation_stats["3mers"][triplett]*len(mutation_stats["cases"][key]) > 0:
                mutation_stats["del-1_3mers"][key]["total"] /= (mutation_stats["3mers"]["total"]*len(mutation_stats["cases"][key]))
                mutation_stats["ins+1_3mers"][key]["total"] /= (mutation_stats["2mers"]["total"]*len(mutation_stats["cases"][key]))


    # calculate simplified mutation probabilities for indels
    for key1 in mutation_stats["del-1_full"]:
        mutation_stats["del-1_full"][key1] = (np.sum([mutation_stats["observed"][key1]["del-1"][key2] for key2 in mutation_stats["observed"][key1]["del-1"]])
                                              / (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key1])))
        mutation_stats["ins+1_full"][key1] = (np.sum([mutation_stats["observed"][key1]["ins+1"][key2] for key2 in mutation_stats["observed"][key1]["ins+1"]])
                                              / (mutation_stats["bases"]["total"]*len(mutation_stats["cases"][key1])))


    # calculate average gene factor
    mutation_stats["avg. gene factor"] = {}
    for key in mutation_stats["genes"]:
        mutation_stats["avg. gene factor"][key] = np.mean([mutation_stats["genes"][key][gene] for gene in mutation_stats["genes"][key]])

    # print to console and file
    print("<", len(mutation_stats["genes"]["total"]), "sequences of total size", mutation_stats["seq_size"])

    # test whether total probabilities can be reconstructed from single probabilities
    for key1 in ["bases", "2mers", "3mers", "del-1", "del-1_3mers", "ins+1", "ins+1_3mers", "pairs", "pairs_3mers"]:
        if type(mutation_stats[key1][list(mutation_stats[key1].keys())[0]]) != dict:
            value_sum = np.sum([mutation_stats[key1][key2] for key2 in mutation_stats[key1] if key2 != "total"])
            if value_sum != mutation_stats[key1]["total"]:
                print("< error. mismatching probabilities @"+key1+" ("+str(value_sum)+"/"+str(mutation_stats[key1]["total"])+")")

        # removed on 251120 as not mathematically sound (a/A+b/B+c/C+d/D != (a+b+c+d)/(A+B+C+D))
        #else:
        #    for key2 in mutation_stats[key1]:
        #       value_mean = np.mean([mutation_stats[key1][key2][key3] for key3 in mutation_stats[key1][key2] if key3 != "total"])
        #        if value_mean != mutation_stats[key1][key2]["total"]:
        #            print("< mismatching probabilities @"+key1+"/"+key2+" ("+str(value_mean)+"/"+str(mutation_stats[key1][key2]["total"])+")")

    # additional plausibility checks
    for key1 in ["del-1", "ins+1", "pairs"]:
        for key2 in mutation_stats[key1]:
            if round(mutation_stats[key1][key2]["total"], 8) != round(mutation_stats[key1+"_3mers"][key2]["total"], 8):
                print("< error. mismatching total probabilities @"+key1+"/"+key2+" ("+str(mutation_stats[key1][key2]["total"])+"/"+str(mutation_stats[key1+"_3mers"][key2]["total"])+")")

    with open(params["data_dir"]+params["os_sep"]+params["fname"].replace(".txt", ".json"), "w") as f:
        f.write(json.dumps(mutation_stats, indent=4))

    return mutation_stats


# marked function (<-) added on 250618
def evaluate_mutation_stats(mutation_stats, seqs, params):
    # calculate driver genes for testing purposes
    mp = Mask_predictions({**params,
                           **{"apply_3mers"                  : True,
                              "apply_mutation_stats"         : True,
                              "apply_mutation_weights"       : False,
                              "mutation_stats_scale"         : 1}},
                           mutation_stats)
    
    sorted_transcripts = sorted(mutation_stats["genes"]["total"])

    # replacement of gene factor for masking and driver gene calculation might be problematic
    for key1 in mutation_stats["genes"]:
        for key2 in mutation_stats["genes"][key1]:
            mutation_stats["genes"][key1][key2] = 1

    for key in mutation_stats["genes"]:
        mutation_stats["avg. gene factor"][key] = np.mean([mutation_stats["genes"][key][gene] for gene in mutation_stats["genes"][key]])
    
    # conduct masking
    bar = IncrementalBar(set_bar("masking"), max=len(sorted_transcripts))
    for sorted_transcript in sorted_transcripts:
        if sorted_transcript in seqs["transcript id"] and len(seqs.loc[sorted_transcript].loc["cds"]) % 3 == 0:
            mp.calculate_probabilities(seqs.loc[sorted_transcript].loc["cds"], sorted_transcript, "total", "total")

        bar.next()
    bar.finish()

    mp.analyze_probability_trajectories()
    mp.calculate_driver_genes()


def create_ptc_subset(wxs, project, wxs_fname):
    ptcs = pd.DataFrame()
    try:    ptcs = wxs.loc[wxs["HGVSp"].str.contains("Ter", na=False)]
    except: print("< error @extract_ptc_variants. creation of 'Ter'-containing mutation types failed for project:", project, wxs_fname)
    return ptcs


def extract_all_variants(cluster, project_key, cluster_key, params, lock):
    data    = {}

    project = cluster[project_key][cluster_key]["project"]
    for i in range(len(cluster[project_key][cluster_key]["input"])):
        for case_id in cluster[project_key][cluster_key]["input"][i]:
            wxss = []

            for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"])):
                wxs = load_wxs(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"WXS"+params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][j])

                if wxs.shape[0] > 0:
                    wxs["HGVSp"] = wxs["HGVSp"].astype(str)
                    if len(params["wxs_cols"]) == 0:  params["wxs_cols"] = wxs.columns.tolist() # needed if empty dataframe must be printed
                    wxss.append(wxs)

                    # define data structure for first wxs
                    if len(data) == 0:
                        cols = [col for col in wxs.columns if col in params["extraction_targets"]]
                        data = {col: [] for col in [*params["ids"], *cols, "cases"]}

            ids = []
            if len(wxss) > 0:
                for j in range(len(wxss)):
                    for k in range(wxss[j].shape[0]):
                        # register variant only once per patient
                        if (type(wxss[j].iloc[k].loc["Gene"]) == str and type(wxss[j].iloc[k].loc["HGVSp"]) == str
                            and wxss[j].iloc[k].loc["Gene"]+"_"+wxss[j].iloc[k].loc["HGVSp"] not in ids):
                            ids.append(wxss[j].iloc[k].loc["Gene"]+"_"+wxss[j].iloc[k].loc["HGVSp"])

                            for col in cols:
                                data[col].append(wxss[j].iloc[k].loc[col])

                            data["cases"].append(len(cluster[project_key][cluster_key]["input"]))
                            if "cluster_key" in data: data["cluster_key"].append(cluster_key)
                            if "project" in data:     data["project"].append(project)
                            if "variant_id" in data:  data["variant_id"].append(wxss[j].iloc[k].loc[params["wxs_identifier"]]
                                                                                +"_"+str(wxss[j].iloc[k].loc["Start_Position"])
                                                                                +"_"+str(wxss[j].iloc[k].loc["End_Position"])
                                                                                +"_"+str(wxss[j].iloc[k].loc["Variant_Type"]))

    return pd.DataFrame(data)


def extract_ptc_variants_(data, ptcs, wxs, params, project, cluster_key, ptc_position, it):
    if params["get_mutation_info"] == False: cols = [*params["ids"], *params["info"]["rna"], *wxs.columns]
    else:                                    cols = [*params["ids"], *params["info"]["rna"], *wxs.columns, *["mutation_info"]]
    ptc_dict = {col: [] for col in cols}
    
    # add RNA-seq columns for later usage
    for col in cols:
        if "RNASEQ_" in col: ptc_dict[col].append([])

    if ptc_position != None: mutation_id = ptcs.iloc[it].loc[params["wxs_identifier"]] + "_" + str(ptc_position)
    if ptc_position == None: mutation_id = ptcs.iloc[it].loc[params["wxs_identifier"]] + "_" + str(ptcs.iloc[it].loc["Start_Position"])+"_"+str(ptcs.iloc[it].loc["End_Position"])+"_"+str(ptcs.iloc[it].loc["Variant_Type"])
    
    ptc_dict["variant_id"].append(mutation_id)
    if data.shape[0] > 0:
        identical_variants = data[data["variant_id"] == mutation_id]

    # create new entry for PTC variant if not already in the dataset (for the same case)
    if data.shape[0] == 0 or (data.shape[0] > 0 and identical_variants[identical_variants["case_id"] == ptcs.iloc[it].loc["case_id"]].shape[0] == 0):
        ptc_dict["project"].append(project)
        ptc_dict["cluster_key"].append(cluster_key)
        if params["get_mutation_info"] == True: ptc_dict["mutation_info"].append(None)

        for col in ptcs.columns:
            ptc_dict[col].append(ptcs.iloc[it].loc[col])

        for i in ptc_dict:
            if len(ptc_dict[i]) != 1: print(i, len(ptc_dict[i]))

        if data.shape[0] > 0: data = pd.concat([data, pd.DataFrame(ptc_dict)], ignore_index=True)
        else:                 data = pd.DataFrame.from_dict(ptc_dict)

        # extract info on additional mutations
        if params["get_mutation_info"] == True:
            mutation_info = get_mutation_info(wxs, params, ptcs.iloc[it], ptc_position)
            if len(mutation_info["Variant_Classification"]) > 0:
                data.at[data.index[data.shape[0]-1], "mutation_info"] = mutation_info

    return data
            

def extract_ptc_variants(cluster, project_key, cluster_key, params, lock):
    data    = pd.DataFrame()

    project = cluster[project_key][cluster_key]["project"]
    for i in range(len(cluster[project_key][cluster_key]["input"])):
        for case_id in cluster[project_key][cluster_key]["input"][i]:
            wxss = []

            for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["WXS"])):
                wxs = load_wxs(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"WXS"+params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][j])

                if wxs.shape[0] > 0:
                    wxs["HGVSp"] = wxs["HGVSp"].astype(str)
                    if len(params["wxs_cols"]) == 0:  params["wxs_cols"] = wxs.columns.tolist() # needed if empty dataframe must be printed
                    if params["use_targets"] == True: wxs = wxs[wxs[params["target_identifier"]["wxs"]].isin(params["targets"])]
                    wxss.append(wxs)

            if len(wxss) > 0:
                # create sub-dataframe of ptc data using "Ter" classifier (should include non-frameshift variants)
                for j in range(len(wxss)):
                    ptcs = create_ptc_subset(wxss[j], project_key, cluster[project_key][cluster_key]["input"][i][case_id]["WXS"][j])
                    
                    for k in range(ptcs.shape[0]):
                        # check whether mutation leads to elongated protein ("ext")
                        #if (len(ptcs.iloc[k].loc["HGVSp"]) > 0 and "ext" not in ptcs.iloc[k].loc["HGVSp"] and "=" not in ptcs.iloc[k].loc["HGVSp"] and "?" not in ptcs.iloc[k].loc["HGVSp"]):
                        if (len(ptcs.iloc[k].loc["HGVSp"]) > 0 and "ext" not in ptcs.iloc[k].loc["HGVSp"] and "=" not in ptcs.iloc[k].loc["HGVSp"] and "?" not in ptcs.iloc[k].loc["HGVSp"]
                            and ptcs.iloc[k].loc["Variant_Classification"] in ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"]):
                            # check whether the gene is present multiple times, if so, only the most downstream ptc variant is selected

                            if check_variants(wxss[j], ptcs, params, k) == True:
                                # get ptc position
                                ptc_position = get_ptc(ptcs.iloc[k].loc["HGVSp"])

                                if ptc_position != None:
                                    data = extract_ptc_variants_(data, ptcs, wxss[j], params, project, cluster_key, ptc_position, k)


                if params["use_targets"] == True:
                    for j in range(len(wxss)):
                        # exclude ptc variants as they have been considered already in the previous step
                        mutants = wxss[j][~wxss[j]["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation"])]

                        # extract only selected mutations if filter is passed
                        if len(params["mutation_targets"]) > 0:
                            mutants = wxss[j][wxss[j][list(params["mutation_targets"].keys())[0]].isin(params["mutation_targets"][list(params["mutation_targets"].keys())[0]])]

                        for k in range(mutants.shape[0]):
                            data = extract_ptc_variants_(data, mutants, wxss[j], params, project, cluster_key, None, k)
                            
    return data


def extract_rnaseq(cluster, data, project_key, cluster_key, params, lock):
    rnacols = ["gene_id", "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded",
               "cnv_total", "cnv_minor", "cnv_avg"]
    project = cluster[project_key][cluster_key]["project"]

    for i in range(len(cluster[project_key][cluster_key]["input"])):
        for case_id in cluster[project_key][cluster_key]["input"][i]:
            rnas = []
            
            for j in range(len(cluster[project_key][cluster_key]["input"][i][case_id]["RNA"])):
                if params["transform_type"] == None:
                    rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+"RNA"
                                   +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)
                    
                if params["transform_type"] != None:
                    rna = load_rna(params["data_dir"]+params["os_sep"]+project+params["os_sep"]+params["transform_type"]
                                   +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j], params)

                rna["gene_id"] = [rna.iloc[k].loc["gene_id"] if "PAR_Y" in rna.iloc[k].loc["gene_id"]
                                  else rna.iloc[k].loc["gene_id"].split(".")[0] for k in range(rna.shape[0])]
                                
                # filter rna for non-negative values if transformed data are used since these can arise from transformation
                if params["transform_type"] != None and rna.shape[0] > 0:
                    try:
                        init_shape = rna.shape[0]
                        rna        = rna[rna[params["rna_criterion"]] != -1]

                        #if rna.shape[0] < init_shape:
                            #print("< removal of missing values reduced expression data from", init_shape, "to", rna.shape[0])
                            #print("  path: ", params["data_dir"]+params["os_sep"]+project+params["os_sep"]+params["transform_type"]
                            #      +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j])

                    except:
                        print("path", params["data_dir"]+params["os_sep"]+project+params["os_sep"]+params["transform_type"]
                               +params["os_sep"]+cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j])
                        print(rna)

                #pre-selection of possible candidates
                if rna.shape[0] > 0:
                    rnas.append(rna[rna["gene_id"].isin(data["Gene"])])

            if len(rnas) > 0:
                lock.acquire()
                try:     print("i", i, "/", len(cluster[project_key][cluster_key]["input"]), data.shape[0], datetime.now(), project_key, cluster_key)
                finally: lock.release()

                for k in range(data.shape[0]):
                    # search for refseq data in ptc variants
                    multiple_seqs = {rnacol: [] for rnacol in rnacols}

                    # initialize as list in order to check whether ptc variants in multiple files share ptc status
                    is_ptc = [False for _ in range(len(rnas))]

                    for l in range(len(rnas)):
                        if rnas[l].shape[0] > 0:
                            ptc_candidates = rnas[l][rnas[l]["gene_id"] == data.iloc[k].loc["Gene"]]
                            if case_id == data.iloc[k].loc["case_id"]: is_ptc[l] = True

                            # add info according to ptc status, inner list to account for overlapping ptcs (probably not occurring)
                            overlapping_seqs = {rnacol: [] for rnacol in rnacols}

                            # expression data entries with the extension "PAR_Y" are not considered, e.g.:
                            # ENSG00000124334.17	    IL9R	protein_coding	14	6	8	0.2061	0.0768	0.0791
                            # ENSG00000124334.17_PAR_Y	IL9R	protein_coding	0	0	0	0.0000	0.0000	0.0000
                            for m in range(ptc_candidates.shape[0]):
                                for rnacol in rnacols:
                                    overlapping_seqs[rnacol].append(ptc_candidates.iloc[m].loc[rnacol])

                            if ptc_candidates.shape[0] != 1:
                                print("ptc candidates", ptc_candidates.shape[0], data.iloc[k].loc["Gene"], cluster[project_key][cluster_key]["input"][i][case_id]["RNA"][j])

                            # write info to multiple seqs dictionary, outer list to account for multiple rna seq datasets
                            for rnacol in rnacols:
                                multiple_seqs[rnacol].append(overlapping_seqs[rnacol])

                    # test to check whether ptc variants in multiple files share ptc status
                    ptc_test = 0
                    for l in range(len(is_ptc)):
                        if is_ptc[l] == False: ptc_test += 1

                    if ptc_test == 0 or ptc_test == len(is_ptc):
                        for rnacol in multiple_seqs:
                            if is_ptc[0] == False: data.at[data.index[k], "RNASEQ_noptc_"+rnacol].append(multiple_seqs[rnacol])
                            else:                  data.at[data.index[k], "RNASEQ_ptc_"+rnacol].append(multiple_seqs[rnacol])

                    else:
                        print("exception. ptc status is unambiguous.")

    return data


def extract(cluster, clusters, proc_index, params, lock, thread_id):
    conversion_dict = {key: literal_eval for key in params["info"]["rna"] if "RNASEQ_" in key}
    folders         = ["ptc_variants", "rnaseq", "calculations"]

    for i in range(proc_index[0], proc_index[1]+1, 1):
        data        = pd.DataFrame()
        project_key = clusters[i]["project_key"]
        cluster_key = clusters[i]["cluster_key"]
        project     = cluster[project_key][cluster_key]["project"]
        print("< thread", thread_id, project_key, cluster_key)

        if os.path.isdir(params["data_dir"]+params["os_sep"]+project) == True:
            calculation_exists = False
            ptc_exists         = False
            rna_exists         = False

            for folder in folders:
                if os.path.isfile(params["data_dir"]+params["os_sep"]+"calculations"+params["os_sep"]+project_key+"_"+cluster_key+".txt") == False:
                    if os.path.isfile(params["data_dir"]+params["os_sep"]+"rnaseq"+params["os_sep"]+project_key+"_"+cluster_key+".txt") == False:
                        if os.path.isfile(params["data_dir"]+params["os_sep"]+"ptc_variants"+params["os_sep"]+project_key+"_"+cluster_key+".txt") == True:
                            data       = pd.read_csv(params["data_dir"]+params["os_sep"]+"ptc_variants"+params["os_sep"]+project_key+"_"+cluster_key+".txt", converters=conversion_dict, delimiter="\t")
                            ptc_exists = True

                    else:
                        data       = pd.read_csv(params["data_dir"]+params["os_sep"]+"rnaseq"+params["os_sep"]+project_key+"_"+cluster_key+".txt", converters=conversion_dict, delimiter="\t")
                        rna_exists = True

                else:
                    calculation_exists = True

                if calculation_exists == False:
                    if folder == "ptc_variants" and "ptc_variants" in params["extraction_stages"] and ptc_exists == False:
                        if params["mode"] == "full_extraction":
                            data = extract_all_variants(cluster, project_key, cluster_key, params, lock)

                        if params["mode"] == "ptc_extraction":
                            data = extract_ptc_variants(cluster, project_key, cluster_key, params, lock)

                    elif folder == "rnaseq" and "rnaseq" in params["extraction_stages"] and rna_exists == False and data.shape[0] > 0:
                        data = extract_rnaseq(cluster, data, project_key, cluster_key, params, lock)

                    elif folder == "calculations" and "calculations" in params["extraction_stages"] and data.shape[0] > 0:
                        data = calculate_expression(data, params)

                    # store data
                    store_data(data, params, project_key, cluster_key, folder)
                    
        del data


def filter_cluster(cluster, filter_mode="rna_and_wxs_present"):
    outfiltered = {project: 0 for project in cluster}
    filtered_cluster = cluster
    for project in cluster:
        for cluster_key in cluster[project]:
            temp = []
            for i in range(len(cluster[project][cluster_key]["input"])):
                for case_id in cluster[project][cluster_key]["input"][i]:
                    if filter_mode == "rna_present":
                        if len(cluster[project][cluster_key]["input"][i][case_id]["RNA"]) > 0:
                            temp.append(cluster[project][cluster_key]["input"][i])

                        else:
                            outfiltered[project] += 1

                    if filter_mode == "rna_and_wxs_present":
                        if len(cluster[project][cluster_key]["input"][i][case_id]["RNA"]) > 0 and len(cluster[project][cluster_key]["input"][i][case_id]["WXS"]) > 0:
                            temp.append(cluster[project][cluster_key]["input"][i])

                        else:
                            outfiltered[project] += 1
                
            filtered_cluster[project][cluster_key]["input"] = temp

    print("< following case ids were filtered out due to lacking files")
    print(json.dumps(outfiltered, indent=4))
    return filtered_cluster


# function to extract any non-nonsense mutations located N-terminally from PTC
def get_mutation_info(wxs, params, ptc_info, ptc_position):
    mutation_info = {"Variant_Classification": [], "Start_Position": [], "End_Position": [], "Reference_Allele": [], "HGVSc": [], "HGVSp": [],
                     "CDS_position": [], "Protein_position": [], "SIFT": [], "PolyPhen": []}
    # changed 241117
    #variants      = wxs[wxs["Gene"] == ptc_info["Gene"]]
    variants      = wxs[wxs[params["wxs_identifier"]] == ptc_info[params["wxs_identifier"]]]

    # iterate over variants and register only those N-terminal of ptc position
    for i in range(variants.shape[0]):
        if "Ter" not in variants.iloc[i].loc["HGVSp"]:
            try:
                if float(variants.iloc[i].loc["Protein_position"].split("/")[0]) <= ptc_position:
                    for key in mutation_info:
                        mutation_info[key].append(variants.iloc[i].loc[key])
            
            except:
                pass

    return mutation_info


def get_ptc(info):
    ptc = None
    # exception required because entry can be NaN
    try:  
        # last argument added to if condition to deal with cases such as Thr471_Ter474delinsArg
        if type(info) == str and len(info) >= 5 and info[0:5] != "p.Ter" and "_Ter" not in info:
            step1 = info.split("Ter") # checks for "Ter" tags
            step2 = step1[0].split("_")  # checks for insertions/ deletions
            # no insertions / deletions
            if len(step2) == 1:
                ptc = get_numbers(step2[0])[0]
            
            # insertions / deletions
            if len(step2) == 2:
                # deletion
                if "del" in step2[1] and "ins" not in step2[1]:
                    ptc = get_numbers(step2[0])[0] + 1

                # insertion or deletion / insertion
                if "ins" in step2[1] or "delins" in step2[1]:
                    ptc   = get_numbers(step2[0])[0]
                    step3 = step2[1].split("ins") # check for no. of insertions
                    if len(step3[1]) == 0: # if zero, Ter was inserted
                        ptc += 1

                    else:
                        if len(step3[1])/3 > int(len(step3[1])/3):
                            print("< error occurred @get_ptc. ptc is not an integer:", len(step3[1])/3, "from", info)

                        else:
                            ptc += int(len(step3[1])/3) + 1


            # account for frameshift variants
            if len(step1[1]) > 0:
                # marked (<-) added / removed on 250529 to correct previous error
                # ptc += get_numbers(step1[1])[0] # <- removed
                ptc += get_numbers(step1[1])[0]-1 # <- added

    except:
        print("< error @get_ptc for entry", info)

    return ptc


def init_cluster(status, params, project_key=None, cluster_key=None):
    load_firehouse = False
    if cluster_key == "Firehouse": load_firehouse = True
    cluster = None

    """
    cluster: {
              project: {
                        cluster:
                            {
                            "input":
                                [
                                    {
                                    case_id:
                                        {
                                        "CNV":        [],
                                        "RNA":        [],
                                        "WXS":        []
                                        }
                                    }
                                ],

                            "project": project,
                            "data": None # placeholder for pandas dataframe
                            }
                        }
             }
    """

    if project_key == "all": cluster = {project_key: {}}
    else:                    cluster = {project: {} for project in status["file_ids"]}

    for project in status["file_ids"]:
        cnf_found = False
        if load_firehouse == True: cnf, cnf_found = load_firehouse_(params["data_dir"]+params["os_sep"]+project, params)
        if project_key != "all":   project_key = project

        for case_id in status["file_ids"][project]:
            case_index   = status["case_ids"][project]["case_id"].index(case_id)
            submitter_id = status["case_ids"][project]["submitter_id"][case_index]

            if cnf_found == True:
                selected_cnf = cnf[cnf["submitter_id"].str.contains(submitter_id)]

                if selected_cnf.shape[0] > 0: cluster_key = "cnmf_" + str(selected_cnf.iloc[-1].loc["cluster"])
                else:                         cluster_key = "cnmf_not_assigned"


            else:
                cluster_key = "all"

            if cluster_key in cluster[project_key].keys():
                cluster[project_key][cluster_key]["input"].append({case_id: status["file_ids"][project][case_id]})

            else:
                cluster[project_key][cluster_key] = {"input": [{case_id: status["file_ids"][project][case_id]}], "project": project, "data": pd.DataFrame()}

    return cluster


# <- added on 251024 to allow random clusters (for efficient calculation)
def randomize_cluster(cluster, max_procs):
    pseudo_clusters    = ["random_cluster"+str(i) for i in range(max_procs)]
    randomized_cluster = {project_key: {pseudo_cluster: {"input": [], "project": "", "data": pd.DataFrame()}
                                        for pseudo_cluster in pseudo_clusters} for project_key in cluster}
    random.seed(int(time.time()))

    for project_key in cluster:
        for cluster_key in cluster[project_key]:
            for case in cluster[project_key][cluster_key]["input"]:
                pseudo_cluster = pseudo_clusters[int(max_procs*random.random())]
                randomized_cluster[project_key][pseudo_cluster]["input"].append(case)
                randomized_cluster[project_key][pseudo_cluster]["project"] = project_key
                #print(project_key, cluster_key, pseudo_cluster, len(cluster[project_key][cluster_key]), "/", len(randomized_cluster[project_key][pseudo_cluster]))

    for project_key in randomized_cluster:
        for cluster_key in randomized_cluster[project_key]:
            print(project_key, cluster_key, len(randomized_cluster[project_key][cluster_key]["input"]))

    return randomized_cluster


def load_df(path, delimiter="\t"):
    df = pd.DataFrame()
    if os.path.isfile(path):
        with open(path, 'r') as _:
            df = pd.read_csv(path, delimiter=delimiter)

    else:
        print("<", path, "not found.")
        exit()
    
    return df


def load_rna(path, params):
    rna = pd.DataFrame()
    if os.path.isfile(path):
        try:
            with open(path, 'r') as _:
                if params["transform_type"] == None:       rna = pd.read_csv(path, delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5])
                elif "RNA_c" in params["transform_type"]:  rna = pd.read_csv(path, delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5])
                else:                                      rna = pd.read_csv(path, delimiter="\t", skiprows=lambda x: x in [0, 2, 3, 4, 5]) #rna = pd.read_csv(path, delimiter="\t")

        except:
            pass

    return rna


def load_status(params, status=None, fname="status.json"):
    try:
        with open(params["data_dir"]+params["os_sep"]+fname, "r") as file:
            status = json.load(file)

    except:
        print("<", params["data_dir"]+params["os_sep"]+fname, "could not be loaded.")
        
    return status


def load_wxs(path):
    wxs = pd.DataFrame()
    if os.path.isfile(path):
        with open(path, 'r') as _:
            wxs = pd.read_csv(path, skiprows=7, delimiter="\t")

    else:
        print("<", path, "not found.")
    
    return wxs


# probably obsolete
def print_data(cluster, params, delimiter="\t", path=None):
    df  = pd.DataFrame()
    for project_key in cluster:
        if (project_key in params["projects"] or params["projects"] == "all"):
            for cluster_key in cluster[project_key]:
                if cluster[project_key][cluster_key]["data"].shape[0] > 0:
                    if df.shape[0] > 0: df = pd.concat([df, cluster[project_key][cluster_key]["data"]], ignore_index=True)
                    else:               df = cluster[project_key][cluster_key]["data"]

                else:
                    print("<", project_key+"_"+cluster_key, "does not contain data")

    if path == None: df.to_csv(params["run_dir"]+params["os_sep"]+params["fname"], index=False, sep=delimiter)
    else:            df.to_csv(path, index=False, sep=delimiter)


def save_status(status, params):
    status = json.dumps(status, indent=4)
    with open(params["data_dir"]+params["os_sep"]+"status.json", "w") as file:
        file.write(status)
    

def store_data(data, params, project_key, cluster_key, folder):
    sub_dir = params["data_dir"]+params["os_sep"]+folder
    if os.path.isdir(sub_dir) == False: os.mkdir(sub_dir)

    if data.shape[0] == 0:
        if params["get_mutation_info"] == False: cols = [*params["ids"], *params["info"]["rna"], *params["wxs_cols"]]
        else:                                    cols = [*params["ids"], *params["info"]["rna"], *params["wxs_cols"], *["mutation_info"]]
        data = pd.DataFrame({col: [] for col in cols})

    data.to_csv(sub_dir+params["os_sep"]+project_key+"_"+cluster_key+".txt", index=False, sep="\t")