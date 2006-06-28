#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/cumsum-test.h"


int cumsum_test(seq_t *seq, double *pvalue, double *param)
{
	int i, z, k, S, s, f;

	double sum1, sum2, sqrtn, t;

	if ( seq->n < CUMSUM_TEST_LENGTH ) {
		fprintf(stderr, "Error[CumSum Test]: Sequence length too short\n");
		return -1;
	}

	for ( i = 0, S = 0, z = 0; i < seq->n; i++ ) {

		S += 2*SEQ(seq, i) - 1;

		if ( abs(S) > z )
			z = abs(S);
	}

	sqrtn = sqrt(seq->n);

	t = ((-1.0*seq->n)/((double)z)+1.0)/4.0;
	s = (int)t;
	t = (seq->n/((double)z)-1.0)/4.0;
	f = (int)t;
	for ( k = s, sum1 = 0.0; k <= f; k++ )
		sum1 += gsl_cdf_ugaussian_P(((double)((4.0*k+1.0)*z))/sqrtn) - gsl_cdf_ugaussian_P(((double)((4.0*k-1)*z))/sqrtn);

	t = ((-1.0*seq->n)/((double)z)-3.0)/4.0;
	s = (int)t;
	t = (((double)seq->n)/((double)z)-1.0)/4.0;
	f = (int)t;
	for ( k = s, sum2 = 0.0; k <= f; k++ )
		sum2 += gsl_cdf_ugaussian_P(((double)((4.0*k+3.0))*z)/sqrtn) - gsl_cdf_ugaussian_P(((double)((4.0*k+1.0)*z))/sqrtn);


	*pvalue = 1.0 - sum1 + sum2;

	return 0;
}
