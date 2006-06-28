#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/frequency-test.h"


int frequency_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i;
	int S;

	if ( seq->n < FREQUENCY_TEST_LENGTH ) {
		fprintf(stderr, "Error[Frequency Test]: Sequence length too short\n");
		return -1;
	}

	for ( i = 0, S = 0; i < seq->n; i++ )
		S += 2*SEQ(seq, i) - 1;

	pvalue[0] = 2.0*gsl_cdf_ugaussian_Q(fabs(S)/sqrt(seq->n));

	return 0;
}
