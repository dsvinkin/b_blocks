#ifndef _BAYESIAN_BLOCKS_G_H_
#define _BAYESIAN_BLOCKS_G_H_

#include "light_curve.h"

struct bb_output_params{
	vector<double> prior_ratio;
	vector<unsigned> number_of_intervals;
   
	void print(const char*);
};

vector<unsigned> do_bayesian_blocks_gaussian(
    light_curve& lc, 
	double ncp_prior, 
	unsigned& n_int);

void rebin_lc_gauss(
    const light_curve& lc,
    const vector<unsigned>& arr_of_change_points,
	light_curve& lc_rebin);

#endif
