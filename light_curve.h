#ifndef _LIGHT_CURVE_H_
#define _LIGHT_CURVE_H_

#include <vector>
#include <string>

using namespace std;

struct b_block
{
	double Ti,
		   Tf,
		   dT,
		   rate,
		   rate_error;
};

struct light_curve {
	vector<double> arr_t_i,
				   arr_t_f,
				   arr_counts,
                   arr_D_counts;
    unsigned size() const;
    double get_min_res() const;
	void print(const string& fname) const;
    void print_thi(const string& fname) const;
	void print_rate(const string& fname) const;
};

struct cells{
	vector<unsigned> arr_cell_size,
				     arr_cell_counts;
	void print(const string& fname) const;
};

struct gaussian_cells{
	vector<double> arr_cell_a,
				   arr_cell_b,
				   arr_cell_c;
	void print(const string& fname) const;
};

#endif
