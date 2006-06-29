#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/notm-test.h"


int notm_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int i, j, W, N, M;
	double mu, sigma2, X;

	seq_t *B = ((notm_param_t *)param)->seq;

	M = ((notm_param_t *)param)->M;

	if ( seq->n < NOTM_TEST_LENGTH ) {
		fprintf(stderr, "Error[NOTM Test]: Sequence length too short\n");
		return -1;
	}

	if ( M <= 0.01*seq->n ) {	
		fprintf(stderr, "Error[NOTM Test]: Invalid parameter\n");
		return -1;
	}

	mu = (double)(M - B->n + 1)/pow(2.0, B->n);
	sigma2 = ((double)M)*(1.0/pow(2.0, B->n) - (double)(2*B->n-1)/pow(2.0, 2.0*B->n));
	N = (int)(seq->n/M);

	for ( i = 0, X = 0.0; i < N; i++ ) {

		for ( j = 0, W = 0; j <= M - B->n; j++ ) {

			if ( seq_cmp(B, seq, i*M + j, B->n) == 0 ) {
				W++;
				j += B->n - 1;			
			}
		}
		
		X += pow(((double)W) - mu, 2.0)/sigma2;
	}

	*pvalue = gsl_cdf_chisq_Q(X, N);

	return 0;
}
