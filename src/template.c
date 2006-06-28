#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/-test.h"


int _test(seq_t *seq, double *pvalue)
{
	unsigned int i;
	int S;

	if ( seq->n < _TEST_LENGTH ) {
		fprintf(stderr, "Error[Test]: Sequence length too short\n");
		return -1;
	}


	//*pvalue = 2*gsl_cdf_ugaussian_Q(fabs());
	//*pvalue = gsl_cdf_chisq_Q(X, k);

	return 0;
}
