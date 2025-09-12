import os
import pandas as pd
import platform

from prepare_data_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def extract_sequence(transcript, genome):
    all_cds = ""
    fails   = []

    exonsstart = [int(i) for i in transcript.loc["exonsstart"].split(",") if len(i) > 0] # condition required because last letter is a comma
    exonsend   = [int(i) for i in transcript.loc["exonsend"].split(",") if len(i) > 0] # condition required because last letter is a comma

    # check whether chromosome is present
    chr = transcript.loc["chr"]
    if chr not in transcript.loc["chr"]: fails.append("chromosome_not_found")

    if transcript.loc["strand"] == "+": i = 0
    if transcript.loc["strand"] == "-": i = transcript.loc["exons"]-1
    
    while ((transcript.loc["strand"] == "+" and i < transcript.loc["exons"]) or (transcript.loc["strand"] == "-" and i >= 0)) and len(fails) == 0:
    	# calculations for exons with no CDS (UTR-only)
        if (transcript.loc["cdsstart"] == transcript.loc["cdsend"] # as defined in the knownGene database
            or (transcript.loc["cdsstart"] != transcript.loc["cdsend"]
            and (exonsend[i] < transcript.loc["cdsstart"] or exonsstart[i] > transcript.loc["cdsend"]))):
            pass

        # calculations for exons with CDS
        else:
            # check if exon contains non-cds
            if exonsstart[i] < transcript.loc["cdsstart"]: exon_cds_start = transcript.loc["cdsstart"]
            else:                                          exon_cds_start = exonsstart[i]
            if exonsend[i] > transcript.loc["cdsend"]:     exon_cds_end   = transcript.loc["cdsend"]
            else:                                          exon_cds_end   = exonsend[i]
        
            if transcript.loc["strand"] == "+":            all_cds       += genome[chr][exon_cds_start:exon_cds_end]
            elif transcript.loc["strand"] == "-":          all_cds       += invert(genome[chr][exon_cds_start:exon_cds_end])
               
        if transcript.loc["strand"] == "+": i += 1
        if transcript.loc["strand"] == "-": i -= 1

    if all_cds[0:3] != "ATG":                                 fails.append("no_start")
    if all_cds[len(all_cds)-3:] not in ["TAA", "TAG", "TGA"]: fails.append("no_stop")
    return all_cds, fails


# parameters (default setting: input from Teran et al. 2021 (mmc4.txt))
params = {
         "data_dir":          parent_dir+r"\data",
         "error_handling":    { # defines handling of errors during sequence information extraction (True means excluding errors)
                              "chromosome_not_found"            : True,
                              "no_start"                        : True,
                              "no_stop"                         : True
                              },
         "extract_expression": True,
         "hg_build":          ["hg38"], #["hg19", "hg38"], priority list for hg builds
         "os_sep":            "\\",
         "outfname":          "hg38_seqs"
         }


if platform.system() == "Windows": params["os_sep"] = "\\"


def main():
    for hg_build in params["hg_build"]:
        print("<", hg_build)

        if hg_build == "hg19":
            hg_build_dir    = "hg19"
            knowngene_fname = "hg19_knownGene.txt"

        if hg_build == "hg38":
            hg_build_dir    = "hg38.p14"
            knowngene_fname = "hg38_knownGene.txt"


        # load genome
        genome = load_split_genome(params["data_dir"]+params["os_sep"]+hg_build_dir, os_sep=params["os_sep"])
        
        # load knowGene dictionary
        with open(params["data_dir"]+params["os_sep"]+knowngene_fname, 'r') as _:
            knowngene = pd.read_csv(params["data_dir"]+params["os_sep"]+knowngene_fname, delimiter=",", index_col=False).sort_values(by=["chr", "cdsstart"])
        
        # initial container for errors
        error_report = {key: {error: [] for error in [*params["error_handling"], "no_error"]} 
                        for key in ["+del1", "+delins1", "+dup1", "+ins1", "+nonsense", "-del1", "-delins1", "-dup1", "-ins1", "-nonsense",
                                    "+del>1", "+delins>1", "+dup>1", "+ins>1", "-del>1", "-delins>1", "-dup>1", "-ins>1", "total", "+total", "-total"]} 
        
        extracted_sequences = {"gene id": [], "transcript id": [], "strand": [], "cds": []}

        bar = IncrementalBar(set_bar("extracting sequences"), max=knowngene.shape[0])
        for i in range(knowngene.shape[0]):
            # extract sequence if coding sequence exists
            if knowngene.iloc[i].loc["cdsstart"] < knowngene.iloc[i].loc["cdsend"]:
                cds, fails = extract_sequence(knowngene.iloc[i], genome)

                if check_fails(fails, params) == False:
                    extracted_sequences["gene id"].append(knowngene.iloc[i].loc["gene id"])
                    extracted_sequences["transcript id"].append(knowngene.iloc[i].loc["transcript id"])
                    extracted_sequences["strand"].append(knowngene.iloc[i].loc["strand"])
                    extracted_sequences["cds"].append(cds)

                #print(i, knowngene.iloc[i].loc["transcript id"], knowngene.iloc[i].loc["strand"], cds, fails)

                # fill error report
                error_report = fill_error_report(error_report, pd.DataFrame(), fails, knowngene.iloc[i].loc["strand"], i)

            bar.next()
        bar.finish()

        print("< error report")
        for key1 in error_report:
            for key2 in error_report[key1]:
                print(key1, key2, len(error_report[key1][key2]))
        
        extracted_sequences = pd.DataFrame(extracted_sequences)
        extracted_sequences.to_csv(params["data_dir"]+params["os_sep"]+params["outfname"]+".txt", index=False, sep=",")


if __name__ == '__main__':
    main()