#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/serial-test.h"


static unsigned int next_motif(seq_t *motif, unsigned int N)
{
	if ( motif->buff[N] == UCHAR_MAX )
		motif->buff[++N] = 1;
	else
		motif->buff[N]++;

	return N;
} 


int serial_test(seq_t *seq, double *pvalue, double *param)
{
	//unsigned int M = param[0];
	unsigned int M;
	unsigned int i, j, k,  m, N, vm;
	unsigned int V[3] = { 0, 0, 0 };
	seq_t *motif;

	double psi[3];
	double stat1, stat2;

	//M = ((int)log2(seq->n)) - 3; 12.15 pour n = 2000*24
	M = 3;
	
/*
	if ( M >= log((double)seq->n)/log(2.0) - 2.0 ) {
		fprintf(stderr, "Error[Serial Test]: Sequence length too short for block size\n");
		return -1;
	}
*/
	//M = log((double)seq->n)/log(2.0) - 3.0;
	motif = seq_new(M);
	if ( !motif ) { 
		fprintf(stderr, "Error[Serial Test]: Memory allocation failed\n");
		return -2;
	}

	for ( k = 0; k < 3; k++ ) {
		
		m = M - k;
		motif->n = m;
		i = 0;
		N = 0;

		memset(motif->buff, 0, motif->size*sizeof(unsigned char)); 

		do {
			for ( j = 0, vm = 0; j < seq->n; j++) {
				if ( seq_cmp(motif, seq, j, m) == 0 )
					vm++;
			}

			V[k] += vm*vm;

			if ( i + 1 < pow(2.0, (double)m) ) {
				i++;
				N = next_motif(motif, N);
			}
			else
				break;
		} while ( 1 );

		
		psi[k] = pow(2.0, m)/((double)seq->n)*((double)V[k]) - ((double)seq->n);
	}

	stat1 = psi[0] - psi[1]; 
	pvalue[0] = gsl_cdf_chisq_Q(stat1,(int) pow(2.0, M - 1));
	stat2 = psi[0] - 2.0*psi[1] + psi[2];
	pvalue[0] = gsl_cdf_chisq_Q(stat2, (int)pow(2.0, M - 2));

	seq_free(motif);

	return 0;
}
