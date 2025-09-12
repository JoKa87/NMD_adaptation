#include "correct_cnvs.h"


vector<string> chrs {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                     "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                     "chr21", "chr22", "chrX", "chrY"};


vector<vector<vector<double> > > Correct_cnvs::convert_cnv(const Data& cnv)
{
vector<vector<vector<double> > > converted_cnv (24);

    for(unsigned int i = 0; i < cnv.descriptors.size(); i++)
    {
    int match_index = get_match_index(chrs, {cnv.descriptors[i][0]});
    if(match_index != -1) {converted_cnv[match_index].push_back(cnv.values[i]);}
    }

return converted_cnv;
}


CNV Correct_cnvs::find_cnv(const vector<vector<vector<double> > >& converted_cnv, const vector<int>& position, const string& chr)
{
CNV cnv;
int match_index = get_match_index(chrs, {chr});
int i = 0;

    if(match_index != -1 && converted_cnv[match_index].size() > 0)
    {
        do
        {
            if(position[0] >= converted_cnv[match_index][i][0] && position[1] <= converted_cnv[match_index][i][1])
            {
            cnv.total = converted_cnv[match_index][i][2];
            cnv.minor = converted_cnv[match_index][i][3];
            }

            if(cnv.minor > cnv.total)
            {
            cout << "< error @correct_cnvs @find_cnv. inconsistent sizes of total and minor CNV." << endl;
            }

        i++;
        } while(i < converted_cnv[match_index].size() && cnv.total == -1);
    }

    else if(match_index == -1 && params->verbosity == 1)
    {
    cout << "< error @correct_cnvs @find_cnv. no match for chromosome " << chr << endl;
    }

return cnv;
}


float Correct_cnvs::get_cnv_average(const Data& cnv)
{
float aggregated_cnv = 0, aggregated_range = 0;

    for(unsigned int i = 0; i < cnv.values.size(); i++)
    {
    aggregated_cnv   += (cnv.values[i][1]-cnv.values[i][0])*cnv.values[i][2];
    aggregated_range += (cnv.values[i][1]-cnv.values[i][0]);
    //cout << setprecision(13) << "i " << i << " aggregated_cnv " << aggregated_cnv << " aggregated_range " << aggregated_range << " ";
    //cout << cnv.values[i][0] << " " << cnv.values[i][1] << " " << cnv.values[i][2] << endl;
    }

if(aggregated_range == 0) {return -1;}
else                      {return aggregated_cnv/aggregated_range;}
}


void Correct_cnvs::print_rna(const Data& rna, const string& project)
{
ofstream outfile;
string outdir = params->data_dir + "\\" + project + "\\" + params->folder_tag;

bool dir_exists = false;
if(CreateDirectory(outdir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()) {dir_exists = true;}

outfile.open((outdir+"\\"+rna.fname).c_str());

    for(unsigned int i = 0; i < rna.header.size(); i++)
    {
    if(params->print_cnv_column == true && i == 1)                                        {outfile << rna.header[i] << "\t" << "cnv_total" << "\t" << "cnv_minor" << "\t" << "cnv_avg" << endl;}
    if(params->print_cnv_column == false || (params->print_cnv_column == true && i != 1)) {outfile << rna.header[i] << endl;}
    }

    for(unsigned int i = 0; i < rna.descriptors.size(); i++)
    {
        for(unsigned int j = 0; j < rna.descriptors[i].size(); j++)
        {
        outfile << rna.descriptors[i][j] << "\t";
        }

        for(unsigned int j = 0; j < rna.values[i].size(); j++)
        {
        if(params->print_cnv_column == false && j < rna.values[i].size()-1) {outfile << rna.values[i][j] << "\t";}
        else if(params->print_cnv_column == false)                          {outfile << rna.values[i][j] << endl;}
        else if(params->print_cnv_column == true)                           {outfile << rna.values[i][j] << "\t";}
        }

    if(params->print_cnv_column == true) {outfile << cnv_columns[0][i] << "\t" << cnv_columns[1][i] << "\t" << cnv_average << endl;}
    }

outfile.close();
return;
}


void Correct_cnvs::run(Data& rna, const vector<Data>& cnvs, const vector<vector<int> >& positions, const vector<string>& chrs,
                       const vector<string>& gene_ids, const string& project)
{
vector<vector<vector<vector<double> > > > converted_cnvs;
vector<float> cnv_averages;

    for(auto& cnv : cnvs)
    {
    converted_cnvs.push_back(convert_cnv(cnv));

    //*determine average CNVs per covered range
    cnv_averages.push_back(get_cnv_average(cnv));
    }

cnv_average = get_avg(cnv_averages);

if(params->print_cnv_column == true) {cnv_columns = vector<vector<double> > (2, vector<double> (rna.values.size(), -1));}

    for(unsigned int i = 0; i < rna.descriptors.size(); i++)
    {
        //*proceed only if mapping was successful for given gene
        if(chrs[i].size() > 0)
        {
            if(read_line(rna.descriptors[i][0], '.')[0].compare(gene_ids[i]) == 0)
            {
            vector<vector<double> > temp_cnvs (2);
            CNV last_cnv;
            bool cnv_mismatch = false, cnv_missing = false;

                //for(unsigned int j = 0; j < converted_cnvs.size() && cnv_mismatch == false; j++)
                for(unsigned int j = 0; j < converted_cnvs.size(); j++)
                {
                CNV cnv = find_cnv(converted_cnvs[j], positions[i], chrs[i]);
                temp_cnvs[0].push_back(cnv.total);
                temp_cnvs[1].push_back(cnv.minor);

                if(cnv.total != last_cnv.total && last_cnv.total != -1) {cnv_mismatch = true;}
                if(cnv.total == -1)                                     {cnv_missing = true;}
                else                                                    {last_cnv = cnv;}
                //cout << "j " << j << " " << chrs[i] << " " << positions[i][0] << " " << positions[i][1] << " " << last_cnv << " " << cnv << endl;
                }

                //*case that no mismatch is found and all values are covered
                if(cnv_mismatch == false && cnv_missing == false)
                {
                    for(unsigned int j = 0; j < rna.values[i].size(); j++)
                    {
                    //* crucial error corrected on 250311: scaling should be in opposite direction; now omitted completely as correction can be conducted later
                    //*all rna values can also be altered at a later stage (if correct_rnas=false)

                        if(params->correct_rnas == true)
                        {
                        if(get_avg(temp_cnvs[0]) > 0)  {rna.values[i][j] *= 2/float(get_avg(temp_cnvs[0]));}
                        if(get_avg(temp_cnvs[0]) == 0) {rna.values[i][j] = -1;}
                        }

                        if(params->print_cnv_column == true)
                        {
                        cnv_columns[0][i] = float(last_cnv.total)/2;
                        cnv_columns[1][i] = float(last_cnv.minor)/2;
                        }
                    //cout << "i " << i << " " << rna.descriptors[i][0] << " last cnv " << last_cnv << " converted " << float(last_cnv/2) << " " << float(last_cnv) / 2 << endl;
                    }
                }

                //*case that mismatch is found and all values are covered
                else if(cnv_mismatch == true && cnv_missing == false && params->average_mismatch == true)
                {
                    for(unsigned int j = 0; j < rna.values[i].size(); j++)
                    {
                    //*crucial error corrected on 250311: scaling should be in opposite direction
                    //*all rna values can also be altered at a later stage (if correct_rnas=false)

                        if(params->correct_rnas == true)
                        {
                        if(get_avg(temp_cnvs[0]) > 0)  {rna.values[i][j] *= 2/float(get_avg(temp_cnvs[0]));}
                        if(get_avg(temp_cnvs[0]) == 0) {rna.values[i][j] = -1;}
                        }

                        if(params->print_cnv_column == true)
                        {
                        cnv_columns[0][i] = float(get_avg(temp_cnvs[0]))/2;
                        cnv_columns[1][i] = float(get_avg(temp_cnvs[1]))/2;
                        }
                    }
                }

                //*case that mismatch is found and all values are covered but values should not be averaged
                else if(cnv_mismatch == true && cnv_missing == false && params->average_mismatch == false && params->ignore_misses == false)
                {
                    for(unsigned int j = 0; j < rna.values[i].size(); j++)
                    {
                        //*changed on 250311, all rna values can also be altered at a later stage (if correct_rnas=false)
                        if(params->correct_rnas == true)
                        {
                        rna.values[i][j] = -1;
                        }

                        if(params->print_cnv_column == true)
                        {
                        cnv_columns[0][i] = -1;
                        cnv_columns[1][i] = -1;
                        }
                    }
                }

                //*case that not all values are covered
                else if(cnv_missing == true && params->ignore_misses == false)
                {
                    for(unsigned int j = 0; j < rna.values[i].size(); j++)
                    {
                        //*changed on 250311, all rna values can also be altered at a later stage (if correct_rnas=false)
                        if(params->correct_rnas == true)
                        {
                        rna.values[i][j] = -1;
                        }

                        if(params->print_cnv_column == true)
                        {
                        cnv_columns[0][i] = -1;
                        cnv_columns[1][i] = -1;
                        }
                    }
                }
            }

            else
            {
            cout << "< error @correct_cnvs @run. mismatching descriptor with " << read_line(rna.descriptors[i][0], '.')[0] << " and " << gene_ids[i] << endl;
            }
        }

        //*changed on 250311, all rna values can also be altered at a later stage (if correct_rnas=false)
        else if(params->correct_rnas == true && params->ignore_misses == false)
        {
            for(unsigned int j = 0; j < rna.values[i].size(); j++)
            {
            rna.values[i][j] = -1;
            }
        }
    }

print_rna(rna, project);
return;
}
