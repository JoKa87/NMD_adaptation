This repository contains tools developed for the manuscript *Adaptation of Premature Termination Codon Mutations in Cancer to Nonsense-Mediated Decay Susceptibility*.

## Tool overview ##
| tools    | description |
| -------- | ------- |
| apply_cnv  | C++ tool for fast correction for copy number variations (CNV)    |
| apply_lindeboom    | Python tool to extract *NMDetectiveA* predictions from Lindeboom et al.    |
| NMD_activity    | Python tools to extract PTC variants from CPTAC and TCGA, calculate NMD activities, and analyze NMD susceptibility in relationship with various features (e.g. PTC burden)    |
| plots    | Python tools to create manuscript figures    |
| prepare_cuomo_data    | Python tools to extract PTC variants and NMD-related sequence features (for model training) based on the ASE calls from Cuomo et al.    |
| prepare_data    | Python tools to extract PTC variants and NMD-related sequence features for CPTAC, MSK-CHORD, TCGA, immune inhibition experiments, and the GTeX8 data from Teran et al. used for model training     |
| random_forest    | Python tools to create NMD susceptibily predictor *NMDelphi* and apply predictions to the whole human exome   |
| shared    | Python tools for shared tasks and testing   |
| survival_analysis    | Python tools for survival analysis   |

## Required files ##


### internal data ###

| filename    | description |
| -------- | ------- |
| cptac3_mutation_stats.json  | mutation statistics derived from TCGA    |
| cptac3_variants.txt    | PTC variants of the CPTAC-3 dataset    |
| hg38_NMD_susceptibilities.txt    | susceptibility predictions created with *NMDelphi*    |
| model_test_variants.txt    | PTC variants of test dataset (ASE calls from Cuomo et al.)*    |
| model_training_variants.txt    | PTC variants of training dataset (ASE calls from Teran et al.)*    |
| model_training_variants_full.txt    | PTC variants of training dataset (ASE calls from Teran et al.) with redundant PTC variants (used for estimating maximum achievable performance)*    |
| msk_chord_survival_analysis.txt    | survival data of MSK-CHORD patients alongside info on NMD susceptibility    |
| msk_chord_variants.txt    | PTC variants of MSK-CHORD dataset    |
| msk_mutation_stats.json    | mutation statistics derived from TCGA (site-specific) and MSK-CHORD (gene-specific)    |
| nmd_inhibition_variants.txt    | PTC variants from NMD inhibition experiment    |
| tcga_avg_expressions.txt    | average expressions for genes in TCGA projects    |
| tcga_data_info.json    | overview of TCGA files (CNV, expression, and WXS)    |
| tcga_avg_expressions.txt    | average expressions for genes in TCGA projects    |
| tcga_hla_mutations.txt    | HLA mutations per TCGA patient    |
| tcga_immune_editing.txt    | immune editing scores from Gong and Karchin 2022    |
| tcga_mutation_stats.json    | mutation statistics derived from TCGA only    |
| tcga_nmd_mutations.txt    | NMD mutations per TCGA patient    |
| tcga_nmd_targets_analysis.txt    | contains patient- and project-resolved expression values of bonafide NMD targets alongside mutational information    |
| tcga_variants.txt    | PTC variants of TCGA dataset    |

*values for distance to last EJC changed if genes consist of single exon using shared/tools.py, mode="change_defaults"



### external data ###
| filename    | description |
| -------- | ------- |
| ase_aggregated_by_donor_open_access_lines (folder)  | variant calls from Cuomo et al, use to extract PTC variants using prepare_cuomo_data.py and filtering for heterozygous variants (hpsi_filtering.py)*    |
| cptac3_survival_data.tsv | CPTAC-3 survival data     |
| hg19.fa	    | hg19 sequences, downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/, reformatted using shared/tools.py, mode="convert_genomes"    |
| hg19_NMDetectiveA_Lindeboom_et_al.v2.gtf    | *NMDetectiveA* predictions for hg19, Lindeboom et al. 2019, downloaded from https://figshare.com/articles/dataset/NMDetective/7803398    |
| hg19_knownGene.txt    | hg19 build version, downloaded from UCSC    |
| hg38_knownGene.txt    | hg38 build version, downloaded from UCSC    |
| hg38_knownGene_old.txt    | older hg38 build version, downloaded from UCSC, used for CNV correction    |
| hg38_lindeboom_predictions.txt    | restructured *NMDetectiveA* predictions for hg38, calculated by random_forest/analyze_predictions, mode="convert_lindeboom_predictions"    |
| hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf    | *NMDetectiveA* predictions for hg38, Lindeboom et al. 2019, downloaded from https://figshare.com/articles/dataset/NMDetective/7803398    |
| hg38.p14.fa    | hg38 sequences, downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p14/, reformatted using shared/tools.py, mode="convert_genomes"    |
| hg38_seqs.txt    | extracted coding sequences for hg38    |
| mmc4.txt    | immune score predictions from https://bioinformatics.mdanderson.org/estimate/    |
| tcga_immune_scores.txt    | TCGA survival data, Liu et al. 2018    |
| TCGA.PanCancer.onco.genes.OncoVar.tsv    | driver gene classification, Wang et al. 2021    |
| tcga_rna_example    | example data to infer RNA genes    |
| tcga_survival_data.txt    | TCGA survival data, Liu et al. 2018    |
| wt_nmd_targets.txt    | wt NMD targets and non-targets from Wang et al. 2017    |

*downloaded from https://zenodo.org/records/3625024#.Xil-0y2cZ0s