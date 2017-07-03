#include <math.h>
#include <stdio.h>

#include "bayesian_blocks_gaussian.h"

using namespace std;

/*
Файл содержит функции для разбиения временной истории со значениями, распределёнными по Гауссу
в формате gaussian_data (определён в light_curve.h и light_curve.cpp)
на Байесовы блоки.
*/

void lcgauss2cells(
    const light_curve& lc, 
	gaussian_cells& cells, 
	double& data_range)
{
	unsigned nbins=lc.size();
	
	double min_val = lc.arr_counts[0],
		   max_val = lc.arr_counts[0];
		   
	for(unsigned i=0; i<nbins; i++)
	{
		double bin_val = lc.arr_counts[i], 
			   bin_sig = lc.arr_D_counts[i];
		
		bin_sig = bin_sig*bin_sig;
		
		cells.arr_cell_c.push_back( 0.5*(bin_val*bin_val/bin_sig) );
		cells.arr_cell_b.push_back( -(bin_val/bin_sig) );
		cells.arr_cell_a.push_back( 0.5*(1/bin_sig) );
		
		if (bin_val > max_val) max_val = bin_val;
		if (bin_val < min_val) min_val = bin_val;
	}
	
	data_range = max_val - min_val;
}

void logprob_gaussian_lc(
    const vector<double>& arr_sum_a,
	const vector<double>& arr_sum_b,
	const vector<double>& arr_sum_c,
    unsigned n_cells,
	vector<double>& arr_merged_log_prob, 
	double ncp_prior,
	double data_range)
{
	const double logpi = log(M_PI);
	double log_cc = 0;
	if (data_range > 0.0) log_cc = -log(data_range);
	
	//Based on revised Scargle paper, as of 29 Nov 2003 
	for (unsigned i=0; i < n_cells; i++) 
	{
		double log_prob_ratio = log_cc + 0.5*(logpi - log(arr_sum_a[i])) +
			(arr_sum_b[i]*arr_sum_b[i]) / (4.0*arr_sum_a[i]) -
			arr_sum_c[i] - ncp_prior;
			
		arr_merged_log_prob.push_back(log_prob_ratio);
	}
}

void do_find_cp_by_bb_algorithm_for_gaussian_data(
    const gaussian_cells& cells, 
	double ncp_prior,
	double data_range,
	vector<unsigned>& arr_of_change_points,
	vector<double>&   arr_best_log_prob, 
	vector<unsigned>& arr_last_cell_start)
{
	const unsigned n_cells = cells.arr_cell_a.size();
		
	vector<double> arr_sum_a(n_cells, 0.0); 
	vector<double> arr_sum_b(n_cells, 0.0);
	vector<double> arr_sum_c(n_cells, 0.0);
	vector<double> arr_merged_log_prob;  
	
	for(unsigned i=0; i < n_cells; i++)
	{
		/* Accumulate the parameters */
		for(unsigned j=0; j < i; j++) 
		{
			arr_sum_a[j] += cells.arr_cell_a[i];
			arr_sum_b[j] += cells.arr_cell_b[i];
			arr_sum_c[j] += cells.arr_cell_c[i];
		}
		arr_sum_a[i] = cells.arr_cell_a[i];
		arr_sum_b[i] = cells.arr_cell_b[i];
		arr_sum_c[i] = cells.arr_cell_c[i];
		
		arr_merged_log_prob.clear();
		logprob_gaussian_lc(arr_sum_a, arr_sum_b, arr_sum_c, i+1, 
							arr_merged_log_prob, ncp_prior, data_range);
							
		unsigned idx_of_max_probability = 0;
		arr_best_log_prob.push_back(arr_merged_log_prob[0]);
		
		if (i > 0) 
		{
			for(unsigned j=1; j<=i; j++) 
			{
				double temp = arr_best_log_prob[j-1] + arr_merged_log_prob[j];
				if (temp > arr_best_log_prob[i]) 
				{
					arr_best_log_prob[i] = temp;
					idx_of_max_probability = j;
				}
			}
		}
//printf("i=%d cp=%d\n",i, imaxer); 
        // Record the new best position
		arr_last_cell_start.push_back(idx_of_max_probability);
	}

#if 0
  /* Debugging output to a file */
  { 
    FILE *out = fopen("bb_data.dat", "w");
	fprintf(out, "CellIdx LastStart BestLP Mergered\n");
    for (unsigned i=0; i<n_cells; i++) {
      fprintf(out, "%d %d %8.3f %8.3f\n", i, arr_last_cell_start[i], arr_best_log_prob[i], arr_merged_log_prob[i]);
    }
    fclose(out);
  }
#endif

	// Calculate number of change points 
	unsigned ncp = 2;
	unsigned index = arr_last_cell_start[n_cells-1];
	
//printf("index=%u\n", index);
	while (index > 1) 
	{
		ncp ++;
		index = arr_last_cell_start[index-1];
	}

	// Create output array of change points  
	unsigned i_cp = ncp-1;
	arr_of_change_points.assign(ncp, 0);
	
	arr_of_change_points[i_cp--] = n_cells;
	index = arr_last_cell_start[n_cells-1];
	while (index > 1) 
	{
		arr_of_change_points[i_cp--] = index;
		index = arr_last_cell_start[index-1];
	}
}

vector<unsigned> do_bayesian_blocks_gaussian(
    light_curve & lc, 
	double ncp_prior, 
	unsigned& number_of_change_points)
{
	gaussian_cells lc_cells;
	double data_range = 0.0;
	lcgauss2cells(lc, lc_cells, data_range);
	
	vector<unsigned> arr_of_change_points,
					 arr_last_cell_start;
				  
	vector<double>   arr_best_log_prob;
	
	do_find_cp_by_bb_algorithm_for_gaussian_data(
	    lc_cells, ncp_prior, data_range, // inputs
		arr_of_change_points, // outputs
		arr_best_log_prob, 
		arr_last_cell_start);

	number_of_change_points=arr_of_change_points.size();
	
    return arr_of_change_points;
}

void rebin_lc_gauss(const light_curve& lc,
                    const vector<unsigned>& arr_of_change_points,
					light_curve& lc_rebin)
{
	unsigned number_of_cp = arr_of_change_points.size();
	
	for (unsigned j=0; j < number_of_cp - 1; j++) 
	{
		/* Compute cumulants needed to get unbiased estimators for the
		   mean and sigma of the points within the block. */
		lc_rebin.arr_t_i.push_back(lc.arr_t_i[arr_of_change_points[j]]);
		lc_rebin.arr_t_f.push_back(lc.arr_t_f[arr_of_change_points[j+1]-1]);
		
		double sum_of_dispersions=0.0,
		       sum_of_values=0.0;
		
		for (unsigned i = arr_of_change_points[j]; i < arr_of_change_points[j+1]; i++) 
		{
			sum_of_dispersions+=lc.arr_D_counts[i]*lc.arr_D_counts[i];
			sum_of_values+=lc.arr_counts[i];
		}

		unsigned n_bins = arr_of_change_points[j+1] - arr_of_change_points[j];
		lc_rebin.arr_counts.push_back( sum_of_values/double(n_bins) );
		lc_rebin.arr_D_counts.push_back( sqrt(sum_of_dispersions)/double(n_bins) );
	}
}
