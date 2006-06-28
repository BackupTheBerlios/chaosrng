#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/vranex-test.h"


int vranex_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i, J;
	int S, N, etat = 4;

	double X;

	if ( seq->n < VRANEX_TEST_LENGTH ) {
	   fprintf(stderr, "Error[VRanEx Test]: Sequence length too short\n");
	   *pvalue = 0.0;
	   return -1;
   	}

	for ( i = 0, J = 0, S = 0, N = 0; i < seq->n; i++ ) {

		S += 2*SEQ(seq, i) - 1;

		if ( S == 0 )
			J++;
		else if ( S == etat )
			N++;
	}

	if ( S != 0 )
		J++;

	X = fabs(((double)X) - ((double)J))/sqrt(((double)J)*2.0*(4.0*fabs(etat)-2.0));


	*pvalue = 2.0*gsl_cdf_ugaussian_Q(X);

	return 0;
}
