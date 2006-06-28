#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/frequencywb-test.h"


int frequencywb_test(seq_t *seq, double *pvalue, double *param)
{
	//unsigned int M = param[0];
	unsigned int M = 100;
	unsigned int i, j, N;
	int S;
	double X;

	if ( seq->n < FREQUENCYWB_TEST_LENGTH ) {
		fprintf(stderr, "Error[FrequencyWB Test]: Sequence length too short\n");
		return -1;
	}

	N = (int)(seq->n/M);

	for ( i = 0, X = 0; i < N; i++ ) {
		S = 0;
		for ( j = 0; j < M; j++)
			S += SEQ(seq, i*M + j);

		X += pow(((double)S)/((double)M) - 0.5, 2.0);
	}

	pvalue[0] = gsl_cdf_chisq_Q(4.0*((double)M)*X, N);

	return 0;
}
