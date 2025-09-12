#include "load.h"


vector<unsigned int> Load::parse_columns(const vector<string>& entries, const vector<string>& selectors, const string& path)
{
vector<unsigned int> selected_index;

    for(unsigned int i = 0; i < selectors.size(); i++)
    {
    int match_index = get_match_index(entries, {selectors[i]});
    if(match_index == -1) {cout << "< error occurred @load @parse_columns. selector " << selectors[i] << " was not found in ";
                           cout << path << endl; exit(0);}
    else                  {selected_index.push_back(match_index);}
    }

return selected_index;
}


Data Load::data(const vector<vector<string> >& selectors, const string& dir, const string& fname, const int& skiprows, const int& classifier_row, const char& separator)
{
Data data;

ifstream file;
file.open((dir+"\\"+fname).c_str());
data.fname = fname;

    if(!(file.is_open()))
    {
    cout << "< error. " << dir+"\\"+fname << " could not be opened." << endl;
    }

	if(file.is_open())
    {
    vector<unsigned int> selected_index;

        for(int i = 0; i < skiprows; i++)
        {
        string firstline;
        getline(file, firstline);
        data.header.push_back(firstline);
        if(i == classifier_row) {selected_index = parse_columns(read_line(firstline, separator), selectors[0], dir+"\\"+fname);}
        }

    string line;

		while(getline(file, line))
        {
        vector<string> entries = read_line(line, separator);
        vector<string> temp_descriptors;
        vector<double> temp_values;

            for(unsigned int i = 0; i < selected_index.size(); i++)
            {
                if(selected_index[i] < entries.size())
                {
                if(selectors[1][i].compare("DESC") == 0) {temp_descriptors.push_back(entries[selected_index[i]]);}
                if(selectors[1][i].compare("VAL") == 0)  {temp_values.push_back(atof(entries[selected_index[i]].c_str()));}
                }

                else
                {
                cout << "< error @load. selected index exceeds entries." << endl;
                cout << "  " << line << endl;
                }
            }

        data.descriptors.push_back(temp_descriptors);
        data.values.push_back(temp_values);
        }

    file.close();
    }

return data;
}


vector<vector<vector<string> > > Load::type_map(const vector<string>& selectors, const string& path)
{
vector<vector<vector<string> > > type_map;

ifstream file;
file.open((params->data_dir+"\\"+path).c_str());

    if(!(file.is_open()))
    {
    cout << "< error. " << path << " could not be opened." << endl;
    }

	if(file.is_open())
    {
    string firstline;
    getline(file, firstline);
    vector<unsigned int> selected_index = parse_columns(read_line(firstline, ','), selectors, params->data_dir+"\\"+path);

    string line;

		while(getline(file, line))
        {
        vector<string> entries = read_line(line, ',');
        //cout << line << " " << entries.size() << endl;

        vector<vector<string> > temp;
        if(entries.size() != 6) {cout << line << endl;}

            for(auto& i : selected_index)
            {
                if(i < entries.size())
                {
                vector<string> sub_entries = read_line(entries[i], ' ');
                temp.push_back(sub_entries);
                //for(auto& entry: sub_entries) {cout << entry << " ";}
                //cout << endl;
                }

                else
                {
                cout << "< error @load. selected index exceeds entries." << endl;
                cout << line << endl;
                }
            }

        type_map.push_back(temp);
        }

    file.close();
    }

return type_map;
}

