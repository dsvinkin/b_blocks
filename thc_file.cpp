#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "thc_file.h"

using namespace std;

string getline(FILE* f)
{
    string s;
    const size_t BUFSIZE = 10;
    char buf[BUFSIZE];

    for (bool end = false; !end;)
    {
        if (!fgets(buf, BUFSIZE, f)) break;
        if (buf[strlen(buf)-1] == '\n')
        {
        	buf[strlen(buf)-1] = '\0';
        	end = true;
        }
        s += buf;
    }
    return s;
}

void split_double(
    const string& str, 
    vector<double>& arr)
{
    const char* p = str.c_str();
    char *q;
    arr.clear();
    for(;;)
    {
        double x = strtod(p, &q);
        if (p==q) break;
        p = q;
        arr.push_back(x);
        if (*p == ':') ++p; // support for colon-delimited fields
   }
}

void split_string(
    const string& str, 
    vector<string>& arr)
{
//разбиваем строку, в которой ключи разделены пробелами
	istringstream iss(str);
	copy(istream_iterator<string>(iss),
		 istream_iterator<string>(),
	     back_inserter<vector<string> >(arr));
}

bool ReadCharRow(
    FILE* f, 
    vector<string>& arr, 
    bool do_rewind=false)
{
	if (feof(f)) return false;
    string s = getline(f);
	split_string(s, arr);
	if (do_rewind) rewind(f);
    return true;
}

bool ReadNumericRow(
    FILE* f, 
    vector<double>& arr, 
    bool do_rewind=false)
{
    if (feof(f)) return false;
    string s = getline(f);
    split_double(s, arr);
    if (do_rewind) rewind(f);
    return true;
}

bool resolve_file_structure(
    FILE* f,
    unsigned &n_channels)
{
// определяем структуру файла
// если есть заголовок, ставим указатель на след строку
// если заголовка нет, просто вычисляем количество каналов

	vector<double> arr_from_line;
	vector<string> arr_header_line;

	bool file_structure_is_resolved = false;
	if (ReadNumericRow(f, arr_from_line, 1) && arr_from_line.size() != 0)
	{
		file_structure_is_resolved = true;
	}
	else if(ReadCharRow(f, arr_header_line, 0) && arr_header_line.size() != 0)
	{
		file_structure_is_resolved = true;
	}
    n_channels = arr_from_line.size() - 2;
	return file_structure_is_resolved;
}

void read_thc_ascii(
    light_curve& lc,
	string file_name,
    bool thiFlag)
{
	/*
     Файл должен содержать три колонки: ti tf и Counts
    */

    FILE* fp=fopen(file_name.c_str(), "r");
	if (fp == NULL)	throw thc_file_exception("Can not open " + file_name);

    unsigned n_channels_nomin = 1;
    if(thiFlag) n_channels_nomin = 2;
    unsigned n_channels = 0;
	if(!resolve_file_structure(fp, n_channels) || n_channels!=n_channels_nomin) 
    {
        string str = "file: " + file_name;
        printf("File structure was not resolved or file contains more then one channel %s\n", str.c_str());
        exit(1);
     //  throw thc_file_exception(str);
    }

	while(1)
	{
		vector<double> arr_from_line;
	    if (!ReadNumericRow(fp, arr_from_line) || arr_from_line.size()!=n_channels_nomin+2)  break;

        lc.arr_t_i.push_back(arr_from_line[0]);
        lc.arr_t_f.push_back(arr_from_line[1]);

        lc.arr_counts.push_back(arr_from_line[2]);
        if(thiFlag){
            lc.arr_D_counts.push_back(arr_from_line[3]); }
        else{ lc.arr_D_counts.push_back(arr_from_line[2]); }
	}
	fclose(fp);
}