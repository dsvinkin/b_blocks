using namespace std;

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "light_curve.h"

void light_curve::print(const string& fname) const
{
	FILE *fp=fopen(fname.c_str(),"w");
//	fprintf(fp, "TBinStart  TBinEnd Counts\n");

	for(unsigned i=0; i<size(); i++)
	{
		fprintf(fp, "%8.3lf %8.3lf %8.3lf %8.3lf\n", 
            arr_t_i[i], arr_t_f[i], arr_counts[i], sqrt(arr_D_counts[i]));
	}
	fclose(fp);
}

unsigned light_curve::size() const
{
    return arr_t_i.size();
}

double light_curve::get_min_res() const
{
    double min_res = arr_t_f[0] - arr_t_i[0];
    for(unsigned i=0; i<size(); i++)
    {
        if(min_res > arr_t_f[i] - arr_t_i[i]) {min_res = arr_t_f[i] - arr_t_i[i];}  
    }
    return min_res;
}

void light_curve::print_rate(const string& fname) const
{
	FILE *fp=fopen(fname.c_str(),"w");
//	fprintf(fp, "Ti  Tf  Rate  RateErr\n");

	for(unsigned i=0; i<size(); i++)
	{
		fprintf(fp, "%8.3lf %8.3lf %8.3lf %8.3lf\n", 
            arr_t_i[i], arr_t_f[i], 
            arr_counts[i]/(arr_t_f[i]-arr_t_i[i]), 
            sqrt(arr_counts[i])/(arr_t_f[i]-arr_t_i[i]) );
	}
	fclose(fp);
}

void light_curve::print_thi(const string& fname) const
{
    FILE *fp=fopen(fname.c_str(),"w");
//    fprintf(fp, "Ti  Tf  Rate  RateErr\n");

    for(unsigned i=0; i<size(); i++)
    {
        fprintf(fp, "%8.3lf %8.3lf %8.3lf %8.3lf\n",
            arr_t_i[i], arr_t_f[i],
            arr_counts[i],
            arr_D_counts[i]);
    }
    fclose(fp);
}

void cells::print(const string& fname) const
{
	FILE *fp=fopen(fname.c_str(),"w");
	fprintf(fp, "CellIdx CellSize CellCounts\n");
	
	unsigned n_bins=arr_cell_size.size();
	for(unsigned i=0; i<n_bins; i++)
	{
		fprintf(fp, "%u %u %u\n", i+1, arr_cell_size[i], arr_cell_counts[i]);
	}
	fclose(fp);
}

void gaussian_cells::print(const string& fname) const
{
	FILE *fp=fopen(fname.c_str(),"w");
	fprintf(fp, "CellIdx  Ca Cb Cc\n");
	
	unsigned n_bins=arr_cell_a.size();
	for(unsigned i=0; i<n_bins; i++)
	{
		fprintf(fp, "%u %8.3lf %8.3lf %8.3lf\n", 
            i+1, arr_cell_a[i], arr_cell_b[i], arr_cell_c[i]);
	}
	fclose(fp);
}