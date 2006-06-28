#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/notm-test.h"

/*
 * TODO:
 * - Passer B en argument
 * - M > 0.01*n
 */

int notm_test(seq_t *seq, double *pvalue)
{
	unsigned int i, j, W, N, M = 10, m = 3;
	double mu, sigma2, X;

	unsigned char buff;
	unsigned char tab[3] = {0,0,1};
	seq_t B = {&buff, 3, 1};
	seq_init_with_uchar(&B, tab, 3);
	seq_print(&B, 3);

	/*
	if ( seq->n < NOTM_TEST_LENGTH ) {
		fprintf(stderr, "Error[NOTM Test]: Sequence length too short\n");
		return -1;
	}
	*/

	mu = (double)(M - m + 1)/pow(2.0, m);
	sigma2 = ((double)M)*(1.0/pow(2.0, m) - (double)(2*m-1)/pow(2.0, 2.0*m));
	N = (int)(seq->n/M);

	for ( i = 0, X = 0.0; i < N; i++ ) {

		for ( j = 0, W = 0; j <= M - m; j++ ) {

			if ( seq_cmp(&B, seq, i*M + j, m) == 0 ) {
				W++;
				j += m - 1;			
			}
		}
		
		X += pow(((double)W) - mu, 2.0)/sigma2;
	}

	*pvalue = gsl_cdf_chisq_Q(X, N);

	return 0;
}
