#include <stdio.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/otm-test.h"

/*
 * TODO:
 * - Passer B en argument
 * - verification des parametres
 */

int otm_test(seq_t *seq, double *pvalue)
{
	unsigned int i, j, W, N, K = 5, M = 10, m = 2;
	double eta, X;
	double pi[K+1], Spi;
	unsigned int V[K+1];

	unsigned char buff;
	unsigned char tab[2] = {1,1};
	seq_t B = {&buff, 2, 1};
	seq_init_with_uchar(&B, tab, 2);

	/*
	if ( seq->n < OTM_TEST_LENGTH ) {
		fprintf(stderr, "Error[NOTM Test]: Sequence length too short\n");
		return -1;
	}
	*/

	N = (int)(seq->n/M);
	memset(V, 0, (K + 1)*sizeof(unsigned int));

	for ( i = 0; i < N; i++ ) {

		for ( j = 0, W = 0; j <= M - m; j++ ) {

			if ( seq_cmp(&B, seq, i*M + j, m) == 0 )
				W++;				
		}

		if ( W >= K )
			V[K]++;
		else
			V[W]++;
	}
	
	eta = (double)(M - m + 1)/pow(2.0, m + 1);
	Spi = pi[0] = exp(-1.0*eta);
	X = pow(((double)V[0]) - ((double)N)*pi[0], 2.0)/(((double)N)*pi[0]);

	for ( i = 1; i < K; i++ ) {


		for ( j = 1, pi[i] = 0.0; j <= i; j++ )
			pi[i] += pow(eta, j)*tgamma(i)/(tgamma(i-j+1)*pow(tgamma(j), 2.0)*j);

		pi[i] *= exp(-1.0*eta)/pow(2.0, i); 
		Spi += pi[i];

		X += pow(((double)V[i]) - ((double)N)*pi[i], 2.0)/(((double)N)*pi[i]);
	}

	pi[i] = 1.0 - Spi;
	X += pow(((double)V[i]) - ((double)N)*pi[i], 2.0)/(((double)N)*pi[i]);

	*pvalue = gsl_cdf_chisq_Q(X, K);

	return 0;
}
