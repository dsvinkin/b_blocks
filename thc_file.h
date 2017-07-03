#ifndef THC_FILE_H_
#define THC_FILE_H_

#include <sstream>
#include <string>
#include <exception>
#include "light_curve.h"

using namespace std;

class thc_file_exception: public exception
{
    string s;
    const char* what() const throw() { return s.c_str(); }

public:
    thc_file_exception(string str) : s(str) {}
    ~thc_file_exception() throw () {} 
};

inline string int2str(int x)
{
    ostringstream buff;
    buff << x;
    return buff.str();   
}

void read_thc_ascii(
    light_curve& lc,
	string file_name,
    bool thcFlag);

#endif /* THC_FILE_H_ */
