#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/approxentropy-test.h"


static unsigned int next_motif(seq_t *motif, unsigned int N)
{
	if ( motif->buff[N] == UCHAR_MAX )
		motif->buff[++N] = 1;
	else
		motif->buff[N]++;

	return N;
} 


int approxentropy_test(seq_t *seq, double *pvalue, void *param)
{
	//unsigned int M = param[0];
	unsigned int M = 3;
	unsigned int i, j, k,  m, N, vm;
	seq_t *motif;

	double phi[2], stat;


	
	if ( M >= (int)(log((double)seq->n)/log(2.0)) - 2 ) {
		fprintf(stderr, "Error[ApproxEntropy Test]: Invalid parameter\n");
		return -1;
	}

	motif = seq_new(M);
	if ( !motif ) { 
		fprintf(stderr, "Error[ApproxEntropy Test]: Memory allocation failed\n");
		return -2;
	}

	for ( k = 0; k < 2; k++ ) {
		
		m = M + 1 - k;
		motif->n = m;
		i = 0;
		N = 0;
		phi[k] = 0.0;

		memset(motif->buff, 0, motif->size*sizeof(unsigned char)); 

		do {
			for ( j = 0, vm = 0; j < seq->n; j++) {
				if ( seq_cmp(motif, seq, j, m) == 0 )
					vm++;
			}

			phi[k] += ((double)vm)/((double)seq->n)*log(((double)vm)/((double)seq->n));

			if ( i + 1 < pow(2.0, (double)m) ) {
				i++;
				N = next_motif(motif, N);
			}
			else
				break;
		} while ( 1 );
	}

	stat = 2.0*((double)seq->n)*(log(2.0) - phi[1] + phi[0]);
	*pvalue = gsl_cdf_chisq_Q(stat, (int)pow(2.0, M));

	seq_free(motif);

	return 0;
}
