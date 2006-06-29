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
		ret |= SEQ(seq, offset + L - i - 1 ) << i;
	
	return ret;
}

int maurer_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int i, ps, L, Q, K, *p; 
	double sum, sigma;

	double Lparams[11][2] = {
			{5.2177052, 2.954},
			{6.19662507, 3.125},
			{7.1836656, 3.258},
			{8.1764248, 3.311},
			{9.1723243, 3.356},
			{10.170032, 3.384},
			{11.168765, 3.401},
			{12.168070, 3.410},
			{13.167693, 3.416},
			{14.167488, 3.419},
			{15.167379, 3.421}
	};

	if ( seq->n < 1010*6*((unsigned int)pow(2,6)) ) {
		fprintf(stderr, "Error[Maurer Test]: Sequence length too short\n");
		return -1;
	}

	for ( L = 6; L <= 16; L++ ) {
		if ( seq->n <= 1010*(L + 1)*((unsigned int)pow(2,L + 1)) )
			break;
	}

	ps =(int)pow(2,L);
	p = (unsigned int *)malloc(ps*sizeof(unsigned int));
	if ( !p ) {
		fprintf(stderr, "Error[Maurer Test]: Memory allocation failed\n");
		return -2;
	}

	memset(p, 0, ps*sizeof(unsigned int));

	Q = 10*((unsigned int)pow(2,L));

	for ( i = 0; i < Q; i++ )
		*(p + dec(seq, i*L, L)) = i + 1;

	K = (int)(seq->n/L) - Q;

	for ( sum = 0; i < Q + K; i ++ ) {
		sum += log( (double)(i + 1) - (double)(*(p + dec(seq, i*L, L))))/log(2);
		*(p + dec(seq, i*L, L)) = i + 1;
	} 

	sum /= (double)K;
	sigma = (0.7 - 0.8/((double)L) + (4.0 + 32.0/((double)L))*pow(K, -3.0/((double)L))/15.0)*sqrt(Lparams[6 - L][1]/((double)K));

	*pvalue = 2*gsl_cdf_ugaussian_Q(fabs((sum - Lparams[6 - L][0])/sigma));

	free(p);

	return 0;
}
