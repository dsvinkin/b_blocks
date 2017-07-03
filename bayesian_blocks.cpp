#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>
#include <assert.h>

#include "bayesian_blocks.h"

/*
Файл содержит функции для разбиения временной истории с Пуассоновской статистикой отсчётов
в формате light_curve (определён в light_curve.h и light_curve.cpp)
на Байесовы блоки.  
*/

void lc2cells(
    const vector<double>& arr_t_i,
	const vector<double>& arr_t_f,
	const vector<double>& arr_counts,
	vector<unsigned>& arr_cellsizes,
	vector<unsigned>& arr_cellpops,  
	double time_del)
{
	/*
    тут не нужен массив ошибок, тк отсчеты считаются пуассоновскими
	новый массив создается прямым переписыванием
	размеры ячейки измеряются в величинах наименьшего размера бина
	*/
	for (unsigned i=0; i<arr_t_i.size(); i++) 
	{
		arr_cellsizes.push_back(lround((arr_t_f[i] - arr_t_i[i]) / time_del));
		arr_cellpops.push_back(lround(arr_counts[i]));
	}
}

void logprob_lc(
    const vector<unsigned>& arr_cellsizes,
	const vector<unsigned>& arr_cellpops,
	unsigned n_cells,
	vector<double>&   arr_logprob, 
	const double& ncp_prior)
{
	for (unsigned i=0; i<n_cells; i++) 
	{
		double log_prob_ratio = gsl_sf_lngamma(arr_cellpops[i]+1.0) - (arr_cellpops[i]+1.0)*log(double(arr_cellsizes[i]));
		log_prob_ratio -= ncp_prior;
		
//printf("arr_cellpops[i]=%u arr_cellsizes[i]=%u gsl_sf_lngamma(arr_cellpops[i]+1.0)=%8.3lf log_prob_ratio=%8.3lf\n",
//			    arr_cellpops[i], arr_cellsizes[i], gsl_sf_lngamma(arr_cellpops[i]+1.0),  log_prob_ratio);
				
		arr_logprob.push_back(log_prob_ratio);
	}
}

void do_find_cp_by_bb_algorithm(
    const vector<unsigned>& arr_cellsizes, 
	const vector<unsigned>& arr_cellpops, 
	double ncp_prior, 
	vector<unsigned>& arr_of_change_points,
	vector<double>&   arr_best_log_prob, 
	vector<unsigned>& arr_last_cell_start)
{
	
    const unsigned n_cells=arr_cellsizes.size();
  
	vector <unsigned> arr_cumsizes(n_cells,0); 
	vector <unsigned> arr_cumpops (n_cells,0);  
	vector<double> arr_merged_log_prob;  
   
	for (unsigned i=0; i<n_cells; i++) 
	{
		for(unsigned j=0; j<i; j++) 
		{
			arr_cumsizes[j] += arr_cellsizes[i];
			arr_cumpops[j]  += arr_cellpops[i];
//printf("arr_cumsizes[j]=%u arr_cumpops[j]=%u\n",arr_cumsizes[j], arr_cumpops[j]);
		}
		arr_cumsizes[i] = arr_cellsizes[i];
		arr_cumpops[i] = arr_cellpops[i];
//printf("arr_cumsizes[j]=%u arr_cumpops[j]=%u\n\n",arr_cumsizes[i], arr_cumpops[i]);      

		// Compute the cost function for the cumulants 
		arr_merged_log_prob.clear();
		logprob_lc(arr_cumsizes, arr_cumpops, i+1, arr_merged_log_prob, ncp_prior);
//printf("\n\n");
		// Where is the maximum probability in the joint best|merged arrays? 
		unsigned idx_of_max_probability = 0;
		arr_best_log_prob.push_back(arr_merged_log_prob[0]);
		
		if (i > 0) 
		{
			for(unsigned j=1; j<=i; j++) 
			{
				double temp = arr_best_log_prob[j-1] + arr_merged_log_prob[j];
				
				assert(! isnan(temp));
				
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
    for (unsigned i = 0; i < n_cells; i++) {
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
	arr_of_change_points.assign(ncp, 0);
	
	unsigned i_cp = ncp-1;
	arr_of_change_points[i_cp] = n_cells;
//printf("making cp arr: i_cp=%u arr_of_change_points[i_cp]=%u\n", i_cp, arr_of_change_points[i_cp]);
	i_cp--;
	index = arr_last_cell_start[n_cells-1];
	while (index > 1) 
	{
		arr_of_change_points[i_cp] = index;
//printf("making cp arr: i_cp=%u arr_of_change_points[i_cp]=%u\n", i_cp, arr_of_change_points[i_cp]);
		i_cp--;
		index = arr_last_cell_start[index-1];
	}
}

void do_bayesian_blocks(
    light_curve lc, 
	double ncp_prior, 
	light_curve& lc_blocks)
{
	cells lc_cells;
	lc2cells(lc.arr_t_i, lc.arr_t_f, lc.arr_counts, // input
			  lc_cells.arr_cell_size, lc_cells.arr_cell_counts, //output
			  lc.get_min_res()); // input
#if 0  
  //test cells
  printf("lc.min_res= %8.3lf\n", lc.get_min_res());
  lc_cells.print("cells_poisson.dat");
#endif  
	vector<unsigned> arr_of_change_points,
					 arr_last_cell_start;
				  
	vector<double> arr_best_log_prob;
	
	
	do_find_cp_by_bb_algorithm(lc_cells.arr_cell_size, lc_cells.arr_cell_counts, ncp_prior, // inputs
						arr_of_change_points, // outputs
						arr_best_log_prob, 
						arr_last_cell_start);

	unsigned number_of_change_points=arr_of_change_points.size();
//printf("before for 1: number_of_change_points=%u\n", number_of_change_points);	
	for (unsigned i = 0; i < number_of_change_points - 1; i++)
	{
		double cell_size=0.0;
		double cell_counts = 0.0;
//printf("in for 1: i=%u arr_of_change_points[i]=%u arr_of_change_points[i+1]=%u\n", 
//	   i, arr_of_change_points[i],  arr_of_change_points[i+1]);

		for (unsigned j = arr_of_change_points[i]; j < arr_of_change_points[i+1]; j++) 
		{
			  cell_counts += lc.arr_counts[j];
			  cell_size += lc.arr_t_f[j] - lc.arr_t_i[j];
//printf("in for 2: cell_counts=%8.3lf cell_size=%8.3lf\n",cell_counts,  cell_size);
		}
		
		assert(i+1 <= number_of_change_points-1);
		assert(! isnan(cell_counts));
		assert(! isnan(cell_size));
//printf("cell_counts=%8.3lf cell_size=%8.3lf\n",cell_counts,  cell_size);
		assert(! isnan(cell_counts/cell_size));
		
		lc_blocks.arr_t_i.push_back( lc.arr_t_i[arr_of_change_points[i]] );
		lc_blocks.arr_t_f.push_back( lc.arr_t_i[arr_of_change_points[i+1]] );
		lc_blocks.arr_counts.push_back(cell_counts/cell_size);
        lc_blocks.arr_D_counts.push_back(cell_counts/cell_size/cell_size);
//printf("block idx_ti=%u idx_tf=%u\n", arr_of_change_points[i], arr_of_change_points[i+1]);
//printf("block ti=%8.3lf tf=%8.3lf\n", lc.arr_t_i[arr_of_change_points[i]], lc.arr_t_f[arr_of_change_points[i+1]]);
	}
	lc_blocks.arr_t_f[ number_of_change_points-2 ] = lc.arr_t_f[ arr_of_change_points[ number_of_change_points-1 ] - 1 ];
//printf("last block ti=%8.3lf tf=%8.3lf\n", lc.arr_t_i[arr_of_change_points[number_of_change_points-1]], lc.arr_t_f[arr_of_change_points[number_of_change_points-1]]);
}