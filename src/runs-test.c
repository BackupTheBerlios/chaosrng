#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/runs-test.h"


int runs_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int i;
	int ones, V;
	double prop;

	if ( seq->n < RUNS_TEST_LENGTH ) {
		fprintf(stderr, "Error[Runs Test]: Sequence length too short");
		return -1;
	}

	for ( i = 0, ones = 0, V = 1; i < (seq->n - 1); i++ ) {
		ones += SEQ(seq, i);
		V += SEQ(seq, i)^SEQ(seq, (i + 1));
	}

	ones += SEQ(seq, i);
	prop = (double)ones/seq->n;

	*pvalue = 2*gsl_cdf_ugaussian_Q(fabs(V-2*seq->n*prop*(1-prop))/(2*sqrt(seq->n)*prop*(1-prop)));

	return 0;
}
