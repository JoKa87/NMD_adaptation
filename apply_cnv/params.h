#ifndef PARAMS_H
#define PARAMS_H
#include <string>
#include <vector>

using namespace std;


struct Params
{
//*basic settings
int verbosity            = 0;
int threads              = 9;

//*directory for input files
string data_dir          = ""; //* place parent directory of TCGA data here

//*tag for novel folder
string folder_tag        = "RNA_c_ptc"; //* specify subdirectory

//*name of map file
string known_gene_fname  = "hg38_knownGene_old.txt";
string type_map_fname    = "TCGA_type_map_filtered_ptc.txt"; //* type map is C++-readable file map, produced from tcga_data_info.json with NMD_analysis/tcga_tools.py, mode "map_types"

//*specific settings
bool average_mismatch    = true; //* if set to false, mismatching CNV values are treated as misses (see ignore_misses)
bool correct_rnas        = false; //* if set to true, corrections are applied to rna values
bool ignore_misses       = true; //* if set to true, values are not changed if corresponding CNV file is not found or CNV cannot be assigned. if set to false, values are set to '-1'
bool overwrite           = true;
bool print_cnv_column    = true;

vector<string> selected_projects    = {};
vector<string> type_map_selectors   = {"project", "cnv", "rna"};

vector<vector<string> > cnv_selectors        = {{"Chromosome", "Start", "End", "Copy_Number", "Minor_Copy_Number"}, {"DESC", "VAL", "VAL", "VAL", "VAL"}};
vector<vector<string> > known_gene_selectors = {{"gene id", "chr", "exonstart", "exonend"}, {"DESC", "DESC", "VAL", "VAL"}};
vector<vector<string> > rna_selectors        = {{"gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second",
                                                 "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"},
                                                 {"DESC", "DESC", "DESC", "VAL", "VAL", "VAL", "VAL", "VAL", "VAL"}};
};


#endif
