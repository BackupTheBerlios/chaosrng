#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/maurer-test.h"


static unsigned int dec(seq_t *seq, unsigned int offset, unsigned int L)
{
	unsigned int i, ret;

	for ( i = 0, ret = 0; i < L; i++ )
		ret |= SEQ(seq, offset + i ) << i;
	
	return ret;
}

int maurer_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i, ps, L = 6, Q = 10*((unsigned int)pow(2,L)), K;
	unsigned int *p;
	double sum;

	if ( seq->n < 1010*((unsigned int)pow(2,L)) ) { //64640
		fprintf(stderr, "Error[Maurer Test]: Sequence length too short\n");
		return -1;
	}

	ps =(int)pow(2,L);
	p = (unsigned int *)malloc(ps*sizeof(unsigned int));
	if ( !p ) {
		fprintf(stderr, "Error[Maurer Test]: Memory allocation failed\n");
		return -2;
	}

	memset(p, 0, ps*sizeof(unsigned int));

	for ( i = 0; i < Q; i++ )
		*(p + dec(seq, i*L, L)) = i + 1;

	K = (int)((seq->n - Q*L)/L);

	for ( sum = 0; i < Q + K; i ++ ) {
		sum += log( (double)(i + 1) - (double)(*(p + dec(seq, i*L, L))))/log(2);
		*(p + dec(seq, i*L, L)) = i + 1;
	} 

	sum /= (double)K;

	*pvalue = 2*gsl_cdf_ugaussian_Q(fabs((sum - 5.2177052)/sqrt(2.954)));

	free(p);

	return 0;
}
