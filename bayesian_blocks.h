#ifndef _BAYESIAN_BLOCKS_H_
#define _BAYESIAN_BLOCKS_H_

#include "light_curve.h"

void do_bayesian_blocks(struct light_curve lc, 
						double ncp_prior, 
						light_curve& lc_blocks);

#endif
