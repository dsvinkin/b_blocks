#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <getopt.h>

#include "thc_file.h"
#include "bayesian_blocks_gaussian.h"
#include "bayesian_blocks.h"

using namespace std;

struct input_parameters
{
    double ncp_prior;
    string 
        inp_file_name, 
        out_file_name,
        lc_type;
        
        
    void print()
    {
        printf(
            "Input parameters:\nprior ratio:%6.1lf\n"
            "input file: %s\n"
            "output file: %s\n"
            "lc type: %s\n", 
            ncp_prior, inp_file_name.c_str(), out_file_name.c_str(), lc_type.c_str());
    }
};

void read_options(
    int argc, 
    char * const argv[], 
    input_parameters& par_struct);

void usage();

int main(int argc, char *argv[])
{
	input_parameters pars;
	read_options(argc, argv, pars);
	pars.print();

	if(pars.lc_type=="Gauss")
    {
        light_curve lc;
	    printf("reading input file...\n");
	    read_thc_ascii(lc, pars.inp_file_name, 1);
	
	    unsigned number_of_cp = 0;
	    printf("making Bayesian blocks...\n");
	    vector<unsigned> arr_of_change_points = 
            do_bayesian_blocks_gaussian(lc, pars.ncp_prior, number_of_cp);
	    printf("number of change points found: %u\n", number_of_cp - 2);
	
	    light_curve lc_rebin;
	    rebin_lc_gauss(lc, arr_of_change_points, lc_rebin);
	    lc_rebin.print_thi(pars.out_file_name.c_str());
    }
    else if(pars.lc_type=="Poisson")
    {
        light_curve lc;
        read_thc_ascii(lc, pars.inp_file_name, 0);
        printf("making Bayesian blocks...\n");
        light_curve bb_blocks;
        do_bayesian_blocks(lc, pars.ncp_prior, bb_blocks);
        
        bb_blocks.print(pars.out_file_name.c_str());
    }
	
	return EXIT_SUCCESS;
}

void read_options(
    int argc, 
    char * const argv[], 
    input_parameters& par_struct) 
{    
    for (;;)
    {              
        static struct option long_options[] = 
		{            
			{"help",        no_argument, 0, 'h'},  
			{"PriorRatio",  required_argument, 0, 'p'},
			{"FileName",    required_argument, 0, 'f'},
			{"bbFileName",  required_argument, 0, 'o'},
            {"LcType",      required_argument, 0, 't'},                             
			{0, 0, 0, 0}        
		}; 
		
		int c = getopt_long(argc,argv, "hp:f:o:t:", long_options, NULL);        
		if (c==-1) break;    
		
        switch (c) 
		{            
			case 'h': usage();
                      exit(EXIT_SUCCESS);
                      break;
            case 'p': par_struct.ncp_prior = atof(optarg); 
					  break;
			case 'f': par_struct.inp_file_name = optarg;
			          break;
			case 'o': par_struct.out_file_name = optarg;
			          break;
            case 't': par_struct.lc_type = optarg;
                      break;
		}    
	}
	
    if (optind < argc || argc < 5) 
    {
        usage();
        exit(EXIT_FAILURE);
    }
}

void usage()
{
    fprintf(stdout,
        "Usage: bb_gauss [OPTS]\n"  
        "-h, --help         print this help and exit\n"
		"-p, --PriorRatio   prior ratio\n"
		"-f, --FileName     input THC (LcType=Poisson) or THI (LcType=Gauss) file name\n"
		"-o, --bbFileName   output THI file name\n"
        "-t, --LcType       Gauss|Poisson\n"
        );
}


