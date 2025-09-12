#include <iostream>
#include <mutex>

#include "correct_cnvs.h"
#include "load.h"
#include "params.h"
#include "utils.h"

using namespace std;


void __run__(mutex& mtx, int& file_count, Data known_genes, Dict known_gene_dict, vector<vector<vector<string> > > type_map, vector<unsigned int> thread_index, Params params)
{
Correct_cnvs correct_cnvs(&params);
Load load(&params);
Utils utils;
bool known_gene_dict_mapped = false;

vector<string> chrs, gene_ids;
vector<vector<int> > positions;

    for(unsigned int i = thread_index[0]; i <= thread_index[1]; i++)
    {
    mtx.lock();
    cout << "converting file " << file_count << " of " << type_map.size() << "            \r";
    mtx.unlock();

    string project   = type_map[i][0][0];
    bool files_exist = false;

        //*if overwrite mode is false, check whether files were already created
        if(params.overwrite == false)
        {
        files_exist = true;

            for(auto& rna_fname : type_map[i][2])
            {
            ifstream file((params.data_dir+"\\"+project+"\\"+params.folder_tag+"\\"+rna_fname).c_str());
            if(file.good() == false) {files_exist = false;}
            }
        }

        if(files_exist == false && (utils.get_match(params.selected_projects, {project}) == true || params.selected_projects.size() == 0))
        {
        vector<Data> cnvs;

            //*proceed only if both CNV and RNA files are available for case id
            if(type_map[i][1].size() > 0 && type_map[i][2].size() > 0)
            {
                for(auto& cnv_fname : type_map[i][1])
                {
                cnvs.push_back(load.data(params.cnv_selectors, params.data_dir+"\\"+project+"\\CNV_ranges\\", cnv_fname, 1, 0, '\t'));
                //cout << cnv_fname << " cnv " << cnvs[0].descriptors.size() << " " << cnvs[0].values.size() << endl;
                //if(cnv_fname.compare("c7ba090c-98e4-4f0d-bf9a-e7ec40bd7487") == 0) cout << cnv_fname << " cnv " << cnvs[0].descriptors.size() << " " << cnvs[0].values.size() << endl;
                }

            vector<Data> rnas; bool match = false;

                for(auto& rna_fname : type_map[i][2])
                {
                rnas.push_back(load.data(params.rna_selectors, params.data_dir+"\\"+project+"\\RNA\\", rna_fname, 6, 1, '\t'));
                //cout << rna_fname << " rna " << rnas[0].descriptors.size() << " " << rnas[0].values.size() << endl;
                //if(rna_fname.compare("c7ba090c-98e4-4f0d-bf9a-e7ec40bd7487") == 0) cout << rna_fname << " rna " << rnas[0].descriptors.size() << " " << rnas[0].values.size() << endl;
                }

                if(rnas[0].values.size() > 0 && known_gene_dict_mapped == false)
                {
                    for(unsigned int j = 0; j < rnas[0].descriptors.size(); j++)
                    {
                    int match_index = utils.search_dictionary(known_gene_dict, utils.read_line(rnas[0].descriptors[j][0], '.')[0], 3, true);

                        if(match_index == -1)
                        {
                        if(params.verbosity == 1) {cout << "< error. gene id: " << utils.read_line(rnas[0].descriptors[j][0], '.')[0] << " could not be mapped to dictionary." << endl;}
                        chrs.push_back("");
                        gene_ids.push_back("");
                        positions.push_back({-1, -1});
                        }

                        else
                        {
                        chrs.push_back(known_genes.descriptors[match_index][1]);
                        gene_ids.push_back(known_genes.descriptors[match_index][0]);
                        positions.push_back({known_genes.values[match_index][0], known_genes.values[match_index][1]});
                        }
                    }

                //cout << "sizes " << chrs.size() << " " << gene_ids.size() << " " << positions.size() << endl;
                known_gene_dict_mapped = true;
                }

                for(auto& rna : rnas)
                {
                correct_cnvs.run(rna, cnvs, positions, chrs, gene_ids, project);
                }
            }
        }

    file_count++;
    }

return;
}


int main()
{
Params params;
Load load(&params);
Utils utils;

vector<vector<vector<string> > > type_map = load.type_map(params.type_map_selectors, params.type_map_fname);
Data known_genes                          = load.data(params.known_gene_selectors, params.data_dir, params.known_gene_fname, 1, 0);

vector<string> gene_ids;

    for(unsigned int i = 0; i < known_genes.descriptors.size(); i++)
    {
    gene_ids.push_back(known_genes.descriptors[i][0]);
    }

Dict known_gene_dict = utils.create_dict(gene_ids, 3, true, true);
cout << "< known genes dictionary created." << endl;

vector<vector<unsigned int> > thread_index = utils.split_index(params.threads, type_map.size()); //thread=1 instead of params->threads (empirical finding)

vector<thread> threads;
mutex mtx;
int file_count;

    for(unsigned int i = 0; i < thread_index.size(); i++)
    {
    threads.push_back(thread(__run__, ref(mtx), ref(file_count), known_genes, known_gene_dict, type_map, thread_index[i], params));
    }

    for(auto& th: threads)
    {
    th.join();
    }

cout << "< calculation finished." << endl;
utils.stop();
return 0;
}
