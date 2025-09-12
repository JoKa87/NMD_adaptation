#ifndef CORRECT_CNVS_H
#define CORRECT_CNVS_H
#include <fstream>
#include <iostream>
#include <windows.h>

#include "load.h"
#include "params.h"
#include "utils.h"


struct CNV
{
int total = -1;
int minor = -1;
};


class Correct_cnvs : Utils
{
protected:
Params* params;
float cnv_average;
vector<vector<double> > cnv_columns;

public:
Correct_cnvs (Params* updated_params) : Utils() {params = updated_params;}

vector<vector<vector<double> > > convert_cnv(const Data& cnv);

CNV find_cnv(const vector<vector<vector<double> > >& converted_cnv, const vector<int>& position, const string& chr);

float get_cnv_average(const Data& cnv);

void print_rna(const Data& rna, const string& project);

void run(Data& rna, const vector<Data>& cnvs, const vector<vector<int> >& positions, const vector<string>& chrs,
         const vector<string>& gene_ids, const string& project);

~Correct_cnvs() {}
};

#endif // LOAD_COHORTS_H
