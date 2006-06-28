#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/ranex-test.h"


int ranex_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i, J;
	unsigned int V[6] = { 0, 0, 0, 0, 0, 0 };
	int S, N, etat = 4;

	double pi[6] = {.8750000000,.01562500000,.01367187500,.01196289063,.01046752930,.0732727051};
	double X;

	if ( seq->n < RANEX_TEST_LENGTH ) {
	   fprintf(stderr, "Error[RanEx Test]: Sequence length too short\n");
	   *pvalue = 0.0;
	   return -1;
   	}

	for ( i = 0, J = 0, S = 0, N = 0; i < seq->n; i++ ) {

		S += 2*SEQ(seq, i) - 1;

		if ( S == 0 ) {
			if ( N >= 5 )
				V[5]++;
			else
				V[N]++;

			N = 0;
			J++;
		}
		else if ( S == etat )
			N++;
	}

	if ( J < 500 )
		fprintf(stderr, "Warning[RanEx Test]: J=%d<500 Not enough excursions to compute test statistic\n", J);
	fprintf(stdout, " %d", J);
	if ( S != 0 ) {
		J++;
		if ( N >= 5 )
			V[5]++;
		else
			V[N]++;
	}

	for ( i = 0, X = 0.0; i < 6; i++ )
		X += pow(((double)V[i]) - ((double)J)*pi[i], 2.0)/(((double)J)*pi[i]);

	*pvalue = gsl_cdf_chisq_Q(X, 5);

	return 0;
}
