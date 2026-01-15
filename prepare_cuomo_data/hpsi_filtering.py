import gc
import gzip
import numpy as np
import os
import pandas as pd
import time

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# assessment of files downloaded on 250407
# manual download starting from https://www.hipsci.org/#/assays/gtarray, clicking on the desired cell line
# on the appearing link, selecting the row "Exome-seq	Imputed and phased genotypes	ENA	Feeder-free	8", clicking "ENA"
# on the appearing link, selecting the files ending with "pes.vcf.gz" (and occasionally "vcf.gz.tbi") by clicking directly (Button "download selected files" did not work)
# for cell line HPSI0114i-wegi_1, the corresponding files were not found
data_dir        = parent_dir+r"\data\ase_aggregated_by_donor_open_access_lines"
genotype_dir    = parent_dir+r"\data\genotypes"
mode            = "apply_filter" # "apply_filter" "check_format" "test_downloads"
target_dir      = parent_dir+r"\data\genotype_filtering"
variants_path   = parent_dir+r"\data\ptc_variants_sums_hg19_filter8"


def _apply_filter(errors, variants, genotype_path, skiprows_info, cell_line):
    start_time = time.time()
    # load header only to infer dtypes
    genotype_header = pd.read_csv(genotype_dir+"\\"+genotype_path, skiprows=skiprows_info[cell_line], nrows=1, sep="\t")

    genotype_dtype = {
                      "#CHROM":                         str,
                      "POS":                            "int64",
                      "ID":                             str,
                      "REF":                            str,
                      "ALT":                            str,
                      "QUAL":                           str,
                      "FILTER":                         str,
                      "INFO":                           str,
                      "FORMAT":                         str,
                      genotype_header.columns[-1]:      str
                     }

    # check for presence of expected columns
    undetected_cols = [key for key in genotype_dtype if key not in genotype_header.columns]

    if len(undetected_cols) > 0:
        print("< error. several expected columns were not detected in", genotype_path)
        print(undetected_cols)
        exit()

    genotype = pd.read_csv(genotype_dir+"\\"+genotype_path, skiprows=skiprows_info[cell_line], sep="\t", dtype=genotype_dtype)
    print("< loading of", cell_line, "took", time.time()-start_time, "s.")

    # select variants with identical cell type
    selected_variants = variants[variants["ID:cell id mod."] == cell_line]

    # pre-select genotypes with identical position
    start_time = time.time()
    genotype = genotype[genotype["POS"].isin(np.unique(selected_variants["ID:pos"]))]
    print("< pre-selection of", cell_line, "took", time.time()-start_time, "s.")

    for i in range(selected_variants.shape[0]):
        # select for perfect matches
        selected_genotype = genotype[genotype["POS"] == selected_variants.iloc[i].loc["ID:pos"]]
        #print("i1", i, selected_genotype.shape[0], selected_variants.iloc[i].loc["ID:chr"], "/", selected_genotype.iloc[0].loc["#CHROM"])

        match_index = [j for j in range(selected_genotype.shape[0])
                       if str(selected_variants.iloc[i].loc["ID:chr"]) == str(selected_genotype.iloc[j].loc["#CHROM"])
                       and selected_variants.iloc[i].loc["ID:alt"] == selected_genotype.iloc[j].loc["ALT"]] #  <- check for alternative base added on 250521

        if len(match_index) == 0:
            print("< no entry found for", cell_line, ", variant", selected_variants.iloc[i].loc["ID:variant id"])

        # check for size
        if len(match_index) > 1:
            print("< genotype has multiple entries for", cell_line, ", variant", selected_variants.iloc[i].loc["ID:variant id"])
            print(selected_genotype.iloc[match_index])
            exit()

        elif len(match_index) == 1:
            # check for deviating reference
            if selected_genotype.iloc[match_index[0]].loc["REF"] != selected_variants.iloc[i].loc["ID:ref"]:
                print("< mismatching reference for", cell_line, ", variant", selected_variants.iloc[i].loc["ID:variant id"])

            # check for heterozygosity if first format entry is "GT"
            elif selected_genotype.iloc[match_index[0]].loc["FORMAT"].split(":")[0] == "GT":
                gt_info = "" # <- added on 250521

                # two alternative ways of splitting are required as two separators are used ("|" and "/")
                if "|" in selected_genotype.iloc[match_index[0]].loc[cell_line].split(":")[0]:
                    gt_info = selected_genotype.iloc[match_index[0]].loc[cell_line].split(":")[0].split("|")

                elif "/" in selected_genotype.iloc[match_index[0]].loc[cell_line].split(":")[0]:
                    gt_info = selected_genotype.iloc[match_index[0]].loc[cell_line].split(":")[0].split("/")
                
                if len(gt_info) != 2:
                    print("< inconsistent info size for gt info", cell_line, ", variant", selected_variants.iloc[i].loc["ID:variant id"])
                    errors.append(cell_line + " " +
                                  ("".join(str(selected_genotype.iloc[match_index[0]].loc[col])+"\t"
                                          if i < selected_genotype.shape[1]-1 else selected_genotype.iloc[match_index[0]].loc[col]
                                          for i, col in enumerate(selected_genotype.columns))))

                elif gt_info[0] != gt_info[1]:
                    selected_variants.at[selected_variants.index[i], "ID:gt_format"]    = selected_genotype.iloc[match_index[0]].loc["FORMAT"]
                    selected_variants.at[selected_variants.index[i], "ID:gt"]           = selected_genotype.iloc[match_index[0]].loc[cell_line]
                    selected_variants.at[selected_variants.index[i], "ID:heterozygous"] = True

                elif gt_info[0] == gt_info[1]:
                    selected_variants.at[selected_variants.index[i], "ID:gt_format"]    = selected_genotype.iloc[match_index[0]].loc["FORMAT"]
                    selected_variants.at[selected_variants.index[i], "ID:gt"]           = selected_genotype.iloc[match_index[0]].loc[cell_line]
                    selected_variants.at[selected_variants.index[i], "ID:heterozygous"] = False

            else: 
                print("< 'GT' not found in the first position of format column")
                print(selected_genotype.iloc[match_index[0]].loc["FORMAT"])
            
            #print("i2", i, selected_genotype.shape[0], selected_genotype.iloc[match_index[0]].loc["FORMAT"], selected_genotype.iloc[match_index[0]].loc[cell_line], 
            #     selected_variants.at[selected_variants.index[i], "ID:heterozygous"])

    del genotype
    gc.collect()
    return selected_variants, errors


def _check_format(genotype_dir):
    paths      = os.listdir(genotype_dir)
    genotypes  = {path.split(".")[0]: path for path in paths if "gz.tbi" not in path}
    skiprows_info = {}

    for cell_line in genotypes:
        with gzip.open(genotype_dir+"\\"+genotypes[cell_line], "rt") as f:
            lines = f.readlines(20000)

        for i, line in enumerate(lines):
            test_cols = line.replace("\t", "").replace("\n", "")
            ref_cols  = ("#CHROMPOSIDREFALTQUALFILTERINFOFORMAT"+cell_line)

            if test_cols == ref_cols:
                #print("<", cell_line, "@", i)
                skiprows_info[cell_line] = i
                #print(" ", ref_cols)
                #print(" ", test_cols)

    print("<", len(skiprows_info), "column descriptions detected.")
    return skiprows_info


if mode == "apply_filter":
    # load variant info
    variants      = pd.read_csv(variants_path, sep=",")

    skiprows_info = _check_format(genotype_dir)
    paths         = os.listdir(genotype_dir)
    genotypes     = {path.split(".")[0]: path for path in paths if "gz.tbi" not in path}

    variants["ID:chr"]          = [variants.iloc[i].loc["ID:variant id"].split("_")[0] for i in range(variants.shape[0])]
    variants["ID:cell id mod."] = [variants.iloc[i].loc["ID:cell id"].replace("_ptc_variants", "") for i in range(variants.shape[0])]
    variants["ID:pos"]          = [int(variants.iloc[i].loc["ID:variant id"].split("_")[1]) for i in range(variants.shape[0])]
    variants["ID:ref"]          = [variants.iloc[i].loc["ID:variant id"].split("_")[2] for i in range(variants.shape[0])]
    variants["ID:alt"]          = [variants.iloc[i].loc["ID:variant id"].split("_")[3] for i in range(variants.shape[0])] # <- added on 250521
    variants["ID:gt_format"]    = [None for _ in range(variants.shape[0])]
    variants["ID:gt"]           = [None for _ in range(variants.shape[0])]
    variants["ID:heterozygous"] = [None for _ in range(variants.shape[0])]

    filtered_variants = []; errors = []
    for i, cell_line in enumerate(genotypes):
        start_time = time.time()
        current_filtered_variants, errors = _apply_filter(errors, variants, genotypes[cell_line], skiprows_info, cell_line)
        filtered_variants.append(current_filtered_variants)
        print("<", cell_line, "time elapsed", time.time()-start_time, filtered_variants[-1][filtered_variants[-1]["ID:heterozygous"].isna()].shape[0], "/",
        filtered_variants[-1][filtered_variants[-1]["ID:heterozygous"] == True].shape[0], "/",
        filtered_variants[-1][filtered_variants[-1]["ID:heterozygous"] == False].shape[0])
        filtered_variants[-1].to_csv(target_dir+"\\"+variants_path.split("\\")[-1].split(".txt")[0]+"_"+cell_line+"_genotype.txt", index=False)

    if len(errors) > 0:
        print("< errors")
        for i, error in enumerate(errors):
            print(i, error)

    filtered_variants = pd.concat(filtered_variants)
    filtered_variants = filtered_variants.drop(columns=["ID:chr", "ID:cell id mod.", "ID:pos", "ID:ref", "ID:alt"])
    filtered_variants.to_csv(target_dir+"\\"+variants_path.split("\\")[-1].split(".txt")[0]+"_genotype.txt", index=False)

    filtered_variants = filtered_variants[filtered_variants["ID:heterozygous"] == True]
    filtered_variants = filtered_variants.drop(columns=["ID:heterozygous", "ID:gt_format", "ID:gt"])
    filtered_variants.to_csv(target_dir+"\\"+variants_path.split("\\")[-1].split(".txt")[0]+"_genotype_filtered.txt", index=False)


# test whether skiprows is the same for all data
if mode == "check_format":
    skiprows_info = _check_format(genotype_dir)


# test of downloaded files
if mode == "test_downloads":
    paths      = os.listdir(data_dir)
    cell_lines = [path.split(".")[0] for path in paths if "altcount" in path]
    paths      = os.listdir(genotype_dir)
    genotypes  = {path.split(".")[0]: path for path in paths if "gz.tbi" not in path}

    format_test = []
    for i, cell_line in enumerate(cell_lines):
        if cell_line in genotypes:
            entries = genotypes[cell_line].split(".")
            #print(i, cell_line, ("".join(entry for j, entry in enumerate(entries) if j > 0 and j <= 3)))
            format_test.append("".join(entry for j, entry in enumerate(entries) if j > 0 and j <= 3))
        
        else:
            print("<", cell_line, "not found.")

    if len(np.unique(format_test)) > 1:
        print("< multiple formats detected.")
        print(np.unique(format_test))