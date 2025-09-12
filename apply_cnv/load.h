#ifndef LOAD_H
#define LOAD_H
#include <fstream>
#include <iostream>

#include "params.h"
#include "utils.h"


struct Data
{
string fname;
vector<string> header;
vector<vector<string> > descriptors;
vector<vector<double> > values;
};


struct Columns
{
vector<string> descriptor;
vector<unsigned int> index;
};


class Load : Utils
{
protected:
Params* params;


public:
Load (Params* updated_params) : Utils() {params = updated_params;}

vector<unsigned int> parse_columns(const vector<string>& entries, const vector<string>& selectors, const string& path);
vector<vector<vector<string> > > type_map(const vector<string>& selectors, const string& path);
Data data(const vector<vector<string> >& selectors, const string& dir, const string& fname, const int& skiprows, const int& classifier_row, const char& separator=',');

~Load() {}
};

#endif

