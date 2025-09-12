#ifndef UTILS_H
#define UTILS_H
#include <bits/stdc++.h>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "params.h"

using namespace std;


struct Dict
{
vector<vector<string> > keys;
vector<vector<int> > index;
};


class Utils
{
public:
vector<vector<float> >& add(vector<vector<float> >& matrix1, const vector<vector<float> >& matrix2);
Dict create_dict(const vector<string>& input, const int& dict_size=3, const bool& redundant=true, const bool& read_reverse=false, const int& start_index=0);
vector<vector<float> >& divide(vector<vector<float> >& matrix1, const float& divisor);
int faculty(const int& value);
template<typename T> vector<unsigned int> fastsort(const vector<T>& values, bool* indexwise=nullptr, unsigned int* max_memorysize=nullptr);
template<typename T> float get_avg(const vector<T>& values, vector<T>* exceptions=nullptr);
template<typename T> bool get_match(const vector<T>& templates, const vector<T>& targets);
template<typename T> int get_match_index(const vector<T>& templates, const vector<T>& targets);
string get_name_from_path(const string& path);
long long int get_time();
void get_performance(const long long int& time1, const long long int& time2,
                     const int& verbosity = 0, const string& message = "");
template<typename T> float get_sem(const vector<T>& values, const T& avg, vector<T>* exceptions=nullptr);
template<typename T> vector<unsigned int> get_sorted_indices(vector<T> values);
template<typename T> T get_sum(const vector<T>& values, vector<T>* exceptions=nullptr);
template<typename T> float get_var(const vector<T>& values, const T& avg, vector<T>* exceptions=nullptr);
template<typename T> T maximum(const vector<T>& values);
template<typename T> T minimum(const vector<T>& values);
template<typename T> unsigned int get_max_index(const vector<T>& values);
template<typename T> unsigned int get_min_index(const vector<T>& values);
template<typename T> void normalize(vector<T>& values, T* value=nullptr);
vector<string> read_line(string& line, const char& separator=',', const unsigned int& skipcount=0, const unsigned int& index=0);
template<typename T> bool redundancy(const vector<T>& input, const T& value);
template<typename T> vector<T> remove_redundancies(const vector<T>& input);
int search_dictionary(const Dict& dict, const string& target, const int& dict_size=3, const bool& read_reverse=false, const int& start_index=0);
vector<vector<unsigned int> > split_index(const int& splits, const unsigned int& total_entries);
void stop();
template<class T> string to_string(const T& input, bool* is_double=nullptr);
template<class T> vector<vector<T> > transpose(const vector<vector<T> >& matrix);
template<class T> vector<T> unify(const vector<vector<T> >& values);
};


template<typename T> bool Utils::get_match(const vector<T>& templates, const vector<T>& targets)
{
bool match = false;

    for(unsigned int i = 0; i < targets.size() && match == false; i++)
    {
        for(unsigned int j = 0; j < templates.size() && match == false; j++)
        {
            if(targets[i] == templates[j])
            {
            match = true;
            }
        }
    }

return match;
}


template<typename T> int Utils::get_match_index(const vector<T>& templates, const vector<T>& targets)
{
int match_index = -1;

    for(unsigned int i = 0; i < targets.size() && match_index == -1; i++)
    {
        for(unsigned int j = 0; j < templates.size() && match_index == -1; j++)
        {
            if(targets[i] == templates[j])
            {
            match_index = j;
            }
        }
    }

return match_index;
}


template<typename T> unsigned int Utils::get_max_index(const vector<T>& values)
{
unsigned int result;
T max_value;

    for(unsigned int i = 0; i < values.size(); i++)
    {
        if(i == 0 || max_value < values[i])
        {max_value = values[i];
        result = i;}
    }

return result;
}


template<typename T> unsigned int Utils::get_min_index(const vector<T>& values)
{
unsigned int result;
T min_value;

    for(T i = 0; i < values.size(); i++)
    {
        if(i == 0 || result > values[i])
        {min_value = values[i];
        result = values[i];}
    }

return result;
}


template<typename T> T Utils::maximum(const vector<T>& values)
{
T result;
    for(T i = 0; i < values.size(); i++)
    {
        if(i == 0 || result < values[i])
        {result = values[i];}
    }

return result;
}


template<typename T> T Utils::minimum(const vector<T>& values)
{
T result;
    for(T i = 0; i < values.size(); i++)
    {
        if(i == 0 || result > values[i])
        {result = values[i];}
    }

return result;
}


template<typename T> vector<unsigned int> Utils::fastsort(const vector<T>& values, bool* indexwise, unsigned int* max_memorysize)
{
vector<unsigned int> sorted_indices;
T minval, maxval;

    for(unsigned int i = 0; i < values.size(); i++)
    {
    if(i == 0 || minval > values[i]) {minval = values[i];}
    if(i == 0 || maxval < values[i]) {maxval = values[i];}
    }

unsigned int memorysize;
//cout << "minval " << minval << " maxval " << maxval << endl;
    if(indexwise == nullptr)
    {
    memorysize = int(0.1*values.size());
    if(memorysize == 0) {memorysize = 1;}

        if(max_memorysize != nullptr)
        {if(memorysize > *max_memorysize) {memorysize = *max_memorysize;}}
    }

    else if(indexwise != nullptr)
    {
    if(*indexwise == true) {memorysize = maxval-minval;}
    }

vector<vector<T> > sort_matrix (vector<vector<T> > (memorysize+1, (vector<T> (0, 0))));

vector<vector<unsigned int> > sort_matrix_indices (vector<vector<unsigned int> > (memorysize+1, (vector<unsigned int> (0, 0))));

float granularity = float(memorysize)/float(maxval-minval);
if(granularity == 0) {granularity = 1;}

    for(unsigned int i = 0; i < values.size(); i++)
    {
    int index = int((values[i]-minval)*granularity);
    //cout << i << " value " << values[i] << " index " << index <<" " << "granularity " << granularity << " " << maxval << " " << minval << " " << memorysize << endl;
    sort_matrix[index].push_back(values[i]);
    sort_matrix_indices[index].push_back(i);
    }

vector<unsigned int> temp;
    for(unsigned int i = 0; i < sort_matrix.size(); i++)
    {
        if(sort_matrix[i].size() > 1)
        {
        temp = Utils::get_sorted_indices(sort_matrix[i]);

            for(unsigned int j = 0; j < temp.size(); j++)
            {sorted_indices.push_back(sort_matrix_indices[i][temp[j]]);}
        }

        if(sort_matrix[i].size() == 1)
        {sorted_indices.push_back(sort_matrix_indices[i][0]);}
    }

    //for(unsigned int i = 0; i < sorted_indices.size(); i++)
    //{cout << values[sorted_indices[i]] << " ";}
//cout << endl;

//string i = "0"; do{cin>>i;}while(i.compare("0") == 0);

return sorted_indices;
}


template<typename T> float Utils::get_avg(const vector<T>& values, vector<T>* exceptions)
{
T avg = 0;
T sum = 0;

    if(exceptions == nullptr)
    {
        for(unsigned int i = 0; i < values.size(); i++)
        {
        avg += values[i];
        sum++;
        }
    }

    if(exceptions != nullptr)
    {
        for(unsigned int i = 0; i < values.size(); i++)
        {
        bool match = false;

            for(unsigned int j = 0; j < (*exceptions).size() && match == false; j++)
            {
            if((*exceptions)[j] == values[i]) {match = true;}
            }

            if(match == false)
            {avg += values[i];
            sum++;}
        }
    }

if(sum > 0) {avg /= sum;}

return avg;
}


template<typename T> float Utils::get_sem(const vector<T>& values, const T& avg, vector<T>* exceptions)
{
T sem = 0;
T sum = 0;
T var = 0;

    if(exceptions == nullptr)
    {
        for(size_t i = 0; i < values.size(); i++)
        {var += pow((avg-values[i]),2);
        sum++;}
    }

    if(exceptions != nullptr)
    {
        for(size_t i = 0; i < values.size(); i++)
        {
        bool match = false; //cout << "i1 " << i << " " << (*exceptions).size() << endl;

            for(size_t j = 0; j < (*exceptions).size() && match == false; j++)
            {
            //cout << "i " << i << " j " << j << " " << (*exceptions)[j] << " " << values[i] << endl;
                if((*exceptions)[j] == values[i])
                {match = true;}
            }

            if(match == false)
            {var += pow((avg-values[i]),2);
            sum++;}
        }
    }

if(sum > 1) {sem = sqrt(var/(sum-1))/sqrt(sum);}
if(sum <= 1) {sem = 0;}
return sem;
}


template<typename T> vector<unsigned int> Utils::get_sorted_indices(vector<T> values)
{
vector<unsigned int> sorted_indices;

    for(unsigned int i = 0; i < values.size(); i++)
    {sorted_indices.push_back(i);}

bool swapped = true;

    for(unsigned int i = 0; i < values.size() && swapped == true; i++)
    {
    swapped = false;

        for(unsigned int j = 0; j < values.size()-1; j++)
        {
        T value1 = values[values.size()-1-j];
        T value2 = values[values.size()-2-j];
        unsigned int index1 = sorted_indices[sorted_indices.size()-1-j];
        unsigned int index2 = sorted_indices[sorted_indices.size()-2-j];

            if(value1 < value2)
            {values[values.size()-1-j] = value2;
            values[values.size()-2-j] = value1;
            sorted_indices[sorted_indices.size()-1-j] = index2;
            sorted_indices[sorted_indices.size()-2-j] = index1;
            swapped = true;}
        }
    }

return sorted_indices;
}


template<typename T> T Utils::get_sum(const vector<T>& values, vector<T>* exceptions)
{
T sum = 0;

    if(exceptions == nullptr)
    {
        for(unsigned int i = 0; i < values.size(); i++)
        {
        sum += values[i];
        }
    }

    if(exceptions != nullptr)
    {
        for(unsigned int i = 0; i < values.size(); i++)
        {
        bool match = false;

            for(unsigned int j = 0; j < (*exceptions).size() && match == false; j++)
            {
                if((*exceptions)[j] == values[i])
                {
                match = true;
                }
            }

            if(match == false)
            {
            sum += values[i];
            }
        }
    }

return sum;
}


template<typename T> void Utils::normalize(vector<T>& values, T* value)
{
T divisor = 0;

    if(value == nullptr)
    {
        for(unsigned int i = 0; i < values.size(); i++)
        {
        divisor += values[i];
        }
    }

    if(value != nullptr)
    {
        if(*value == 0)
        {
        cout << "< zero division error." << endl;
        }

    divisor = *value;
    }

    for(unsigned int i = 0; i < values.size(); i++)
    {
    if(divisor != 0) {values[i] / divisor;}
    }

return;
}


template<typename T> float Utils::get_var(const vector<T>& values, const T& avg, vector<T>* exceptions)
{
T var = 0;
T sum = 0;

    if(exceptions == nullptr)
    {
        for(size_t i = 0; i < values.size(); i++)
        {
        var += pow((avg-values[i]),2);
        sum++; //cout << "var " << var << " " << avg << " " << values[i] << " " << sum << endl;
        }
    }

    if(exceptions != nullptr)
    {
        for(size_t i = 0; i < values.size(); i++)
        {
        bool match = false;

            for(size_t j = 0; j < (*exceptions).size() && match == false; j++)
            {
                if((*exceptions)[j] == values[i])
                {match = true;}
            }

            if(match == false)
            {var += pow((avg-values[i]),2);
            sum++;}
        }
    }

if(sum > 1) {var /= (float) (sum-1);}
if(sum <= 1) {var = 0;}
return var;
}


template<typename T> bool Utils::redundancy(const vector<T>& input, const T& value)
{
bool redundance = false;

    for(size_t i = input.size()-1; i >= 0 && redundance == false; i--)
    {
        if(input[i] == value)
        {redundance = true;}
    }

return redundance;
}


template<typename T> vector<T> Utils::remove_redundancies(const vector<T>& input)
{
vector<T> reduced_input;
vector<string> converted_input;

    for(auto& i : input)
    {
    converted_input.push_back(to_string(i));
    }

Dict dict = create_dict(converted_input, 3, false);
vector<int> matches;

    for(unsigned int i = 0; i < dict.index.size(); i++)
    {
        for(unsigned int j = 0; j < dict.index[i].size(); j++)
        {
        reduced_input.push_back(input[dict.index[i][j]]);
        }
    }

return reduced_input;
}


template<class T> string Utils::to_string(const T& input, bool* is_double)
{
stringstream inputstream;

    if(is_double == nullptr)
    {inputstream << input;}

    if(is_double != nullptr)
    {inputstream << setprecision(13) << input;}

return inputstream.str();
}


template<class T> vector<vector<T> > Utils::transpose(const vector<vector<T> >& matrix)
{
vector<vector<T> > out;
T value;

    if(matrix.size() > 0 && matrix[0].size() > 0)
    {
    out = vector<vector<T> > (matrix[0].size(), (vector<T> (matrix.size(), value)));

        for(unsigned int i = 0; i < matrix.size(); i++)
        {
            for(unsigned int j = 0; j < matrix[i].size(); j++)
            {
            out[j][i] = matrix[i][j];
            }
        }
    }

    else
    {
    cout << "< transposing matrix not possible. matrix is empty." << endl;
    }

return out;
}


template<class T> vector<T> Utils::unify(const vector<vector<T> >& values)
{
vector<T> unified_vector;
vector<T> temp;

    for(unsigned int i = 0; i < values.size(); i++)
    {
        for(unsigned int j = 0; j < values[i].size(); j++)
        {temp.push_back(values[i][j]);}
    }

vector<unsigned int> sorted_indices = Utils::get_sorted_indices(temp);

    if(temp.size() > 1 && temp[sorted_indices[0]] != temp[sorted_indices[1]])
    {unified_vector.push_back(temp[sorted_indices[0]]);}

    for(unsigned int i = 1; i < sorted_indices.size(); i++)
    {
        if(temp[sorted_indices[i]] != temp[sorted_indices[i-1]])
        {unified_vector.push_back(temp[sorted_indices[i]]);}
    }

return unified_vector;
}

#endif
