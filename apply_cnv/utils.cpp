#include "utils.h"


vector<vector<float> >& Utils::add(vector<vector<float> >& matrix1, const vector<vector<float> >& matrix2)
{
    for(unsigned int i = 0; i < matrix1.size(); i++)
    {
        for(unsigned int j = 0; j < matrix1[0].size(); j++)
        {
        matrix1[i][j] += matrix2[i][j];
        }
    }

return matrix1;
}


Dict Utils::create_dict(const vector<string>& input, const int& dict_size, const bool& redundant, const bool& read_reverse, const int& start_index)
{
Dict dict;
vector<string> alphanumericals_lower {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                      "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
                                      "k", "l", "m", "n", "o", "p", "q", "r", "s", "t",
                                      "u", "v", "w", "x", "y", "z", "-"};
vector<string> alphanumericals_upper {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                      "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
                                      "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                                      "U", "V", "W", "X", "Y", "Z", "-"};

dict.index = vector<vector<int> > (pow(37, dict_size));
dict.keys  = vector<vector<string> > (pow(37, dict_size));

    for(unsigned int i = 0; i < input.size(); i++)
    {
    vector<int> sizes {int(input[i].size())-start_index, dict_size};
    int min_it = minimum(sizes);
    int index = 0;

        if(min_it <= 0)
        {
        cout << "< error occurred @create_dict. selected start index (" << start_index << ") exceeds input string size for input " << input[i] << "." << endl;
        cout << "  type '0' to proceed." << endl;
        stop();
        }

        if(read_reverse == false)
        {
            for(unsigned int j = 0; j < min_it; j++)
            {
            int match_index = get_match_index(alphanumericals_lower, {input[i].substr(start_index+j, 1)});
            if(match_index == -1) {match_index = get_match_index(alphanumericals_upper, {input[i].substr(start_index+j, 1)});}
            index += pow(37, min_it-1-j)*match_index;
            //if(input[i].compare("MCD") == 0) cout << "input " << input[i] << " i " << i << " j " << j << " index " << index << " match_index " << match_index << endl;
            }
        }

        else
        {
            for(unsigned int j = 0; j < min_it; j++)
            {
            int match_index = get_match_index(alphanumericals_lower, {input[i].substr(input[i].size()-1-j-start_index, 1)});
            if(match_index == -1) {match_index = get_match_index(alphanumericals_upper, {input[i].substr(input[i].size()-1-j-start_index, 1)});}
            index += pow(37, min_it-1-j)*match_index;
            //if(input[i].compare("MCD") == 0) cout << "input " << input[i] << " i " << i << " j " << j << " index " << index << " match_index " << match_index << endl;
            }
        }

        if(redundant == true)
        {
        dict.index[index].push_back(i);
        dict.keys[index].push_back(input[i]);
        }

        else if(redundant == false && get_match_index(dict.keys[index], {input[i]}) == -1)
        {
        dict.index[index].push_back(i);
        dict.keys[index].push_back(input[i]);
        }

    //if(input[i].compare("MCD") == 0) {cout << "index size " << dict.index[index].size() << endl;
    //cout << endl; int x = 0; do{cin>>x;} while(x== 0);}
    }

return dict;
}


vector<vector<float> >& Utils::divide(vector<vector<float> >& matrix1, const float& divisor)
{
    if(divisor != 0)
    {
        for(size_t i = 0; i < matrix1.size(); i++)
        {
            for(size_t j = 0; j < matrix1[0].size(); j++)
            {matrix1[i][j] /= divisor;}
        }
    }

    else
    {cout << "< warning. zero division not possible!" << endl;}

return matrix1;
}


int Utils::faculty(const int& value)
{
int out = 0;

    for(unsigned int i = 0; i < value; i++)
    {
    out += i+1;
    }

return out;
}


string Utils::get_name_from_path(const string& path)
{
string filename;
string teststr, testletter;
unsigned int i = 0;
    do
    {
    testletter = path.substr(i, 1);

        if(testletter.compare("\\") != 0 && testletter.compare(".") != 0)
        {
        teststr += testletter;
        }

        if(testletter.compare("\\") == 0)
        {
        teststr = "";
        }

    i++;
    } while(i < path.size() && testletter.compare(".") != 0);

return teststr;
}


long long int Utils::get_time()
{
const auto now = chrono::system_clock::now();
const auto nowAsTimeT = chrono::system_clock::to_time_t(now);
const auto nowMs = chrono::duration_cast<chrono::milliseconds>(now.time_since_epoch());
stringstream nowSs;
nowSs << nowMs.count();
string nowMs_str = nowSs.str();
return atoll(nowMs_str.c_str());
}


void Utils::get_performance(const long long int& time1, const long long int& time2,
                            const int& verbosity, const string& message)
{
    if(verbosity == 1)
    {cout << "< time elapsed for " << message << " " << time2-time1 << " ms." << endl;}
}


vector<string> Utils::read_line(string& line, const char& separator, const unsigned int& skipcount, const unsigned int& index)
{
vector<string> row;
string teststr = "";
unsigned int entries = 0;
bool show = false;

    for(unsigned int i = 0; i < line.size(); i++)
    {
        if(*(&(line[i])) != separator)
        {
        teststr += *(&(line[i]));
        }

        if((*(&(line[i])) == *(&(separator)) || i == line.size()-1) && entries >= skipcount)
        {
        row.push_back(teststr);
        //cout << "line " << line << " ro " << row.size() << " teststr " << teststr << "#" << endl;
        }

        if(*(&(line[i])) == *(&(separator)))
        {
        entries++;
        teststr = "";
        }
    }

if(*(&(line[line.size()-1])) == separator) {row.push_back("");}
return row;
}


int Utils::search_dictionary(const Dict& dict, const string& target, const int& dict_size, const bool& read_reverse, const int& start_index)
{
int search_index = -1;

vector<string> alphanumericals_lower {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                      "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
                                      "k", "l", "m", "n", "o", "p", "q", "r", "s", "t",
                                      "u", "v", "w", "x", "y", "z", "-"};
vector<string> alphanumericals_upper {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                      "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
                                      "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                                      "U", "V", "W", "X", "Y", "Z", "-"};

vector<int> sizes {target.size()-start_index, dict_size};
int min_it = minimum(sizes);
int index = 0;

    if(min_it <= 0)
    {
    cout << "< error occurred @search_dictionary. selected start index (" << start_index << ") exceeds input string size for " << target << "." << endl;
    cout << "  type '0' to proceed." << endl;
    stop();
    }

    if(read_reverse == false)
    {
        for(unsigned int i = 0; i < min_it; i++)
        {
        int match_index = get_match_index(alphanumericals_lower, {target.substr(start_index+i, 1)});
        if(match_index == -1) {match_index = get_match_index(alphanumericals_upper, {target.substr(start_index+i, 1)});}
        index += pow(37, min_it-1-i)*match_index;
        //if(input[i].compare("MCD") == 0) cout << "input " << input[i] << " i " << i << " j " << j << " index " << index << " match_index " << match_index << endl;
        }
    }

    else
    {
        for(unsigned int i = 0; i < min_it; i++)
        {
        int match_index = get_match_index(alphanumericals_lower, {target.substr(target.size()-1-i-start_index, 1)});
        if(match_index == -1) {match_index = get_match_index(alphanumericals_upper, {target.substr(target.size()-1-i-start_index, 1)});}
        index += pow(37, min_it-1-i)*match_index;
        //cout << "target " << target << " i " << i << " " << target.substr(target.size()-1-i-start_index, 1) << " index " << index << " match_index " << match_index << endl;
        }
    }

//cout << "index " << index << " " << dict.keys[index].size() << endl;
int match_index = -1;
if(index < dict.keys.size()) {match_index = get_match_index(dict.keys[index], {target});}
else                         {cout << "< dictionary error occurred." << endl;}
if(match_index != -1)        {search_index = dict.index[index][match_index];}
//cout << dict.index[index][0] << " " << search_index << endl;
//cout << endl; int x = 0; do{cin>>x;} while(x== 0);
return search_index;
}


vector<vector<unsigned int> > Utils::split_index(const int& splits, const unsigned int& total_entries)
{
vector<vector<unsigned int> > splitted_index;
int usable_splits = min(int(splits), int(total_entries));

unsigned int split_size = int(total_entries/usable_splits), thread_index = 0;
splitted_index = vector<vector<unsigned > > (usable_splits, vector<unsigned > (2, 0));
int last_thread_index = -1;

    for(int i = 0; i < usable_splits-1; i++)
    {
    if(last_thread_index == i*split_size) {splitted_index[i][0] = i*split_size+1;}
    else                                  {splitted_index[i][0] = i*split_size;}
    splitted_index[i][1] = (i+1)*split_size-1;
    last_thread_index = splitted_index[i][1];
    cout << "i " << i << " " << splitted_index[i][0] << " " << splitted_index[i][1] << endl;
    }

if(last_thread_index != -1 && last_thread_index == usable_splits-1*split_size) {splitted_index[usable_splits-1][0] = (usable_splits-1)*split_size+1;}
else                                                                           {splitted_index[usable_splits-1][0] = (usable_splits-1)*split_size;}
splitted_index[usable_splits-1][1] = total_entries-1;
cout << splitted_index[usable_splits-1][0] << " " << splitted_index[usable_splits-1][1] << endl;
return splitted_index;
}


void Utils::stop()
{
string i = "0";

    do
    {
    cin >> i;
    } while(i.compare("0") == 0);

return;
}
