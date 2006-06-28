#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl_complex.h>
#include <gsl_complex_math.h>
#include <gsl_fft_real.h>
#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/dft-test.h"

void  __ogg_fdrffti(int n, double* wsave, double* ifac);
void  __ogg_fdrfftf(int n, double* X, double* wsave, double* ifac);

int dft_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i, N;
	double *X, *pX, *wsave, *ifac;
	double d, T;

	/*
	if ( seq->n < DFT_TEST_LENGTH ) {
		fprintf(stderr, "Error[Frequency Test]: Sequence length too short\n");
		return -1;
	}
	*/

	X = (double *)malloc(seq->n*sizeof(double));
	if ( !X ) {
		fprintf(stderr, "Error[DFT Test]: Memory allocation failed\n");
		return -2;
	}

	memset(X, 0, seq->n*sizeof(double));

	wsave = (double *)malloc((2*seq->n+15)*sizeof(double));
	if ( !wsave ) {
		fprintf(stderr, "Error[DFT Test]: Memory allocation failed\n");
		free(X);
		return -2;
	}
	memset(wsave, 0, (2*seq->n+15)*sizeof(double));

        ifac = (double *)malloc(15*sizeof(double));
	if ( !ifac ) {
		fprintf(stderr, "Error[DFT Test]: Memory allocation failed\n");
		free(X);
		free(wsave);	
		return -2;
	}
	memset(ifac, 0, 15*sizeof(double));

	for ( i = 0, pX = X; i < seq->n; i++, pX++ ){	
		 X[i] = (double)(2*SEQ(seq, i) - 1);
	}

        __ogg_fdrffti(seq->n, wsave, ifac);

        __ogg_fdrfftf(seq->n,X,wsave,ifac);

	T = sqrt((double)(3*seq->n));

	N = (fabs(X[0]) < T) ? 1 : 0 ;
	fprintf(stdout, "%f\n", X[0]);
	for ( i = 0; i < seq->n/2-1; i++ ) {
	fprintf(stdout, "%f %f\n", X[2*i+1], X[2*(i+1)]);
		if ( sqrt(pow(X[2*i+1], 2.0) + pow(X[2*i+2], 2.0)) < T)
			N++;
	}

	fprintf(stdout, "%d\n", N);
	
	
	d = (((double)N)-(0.95*seq->n)/2.0)/sqrt(seq->n*0.95*0.05/2.0);


	*pvalue = 2.0*gsl_cdf_ugaussian_Q(fabs(d));

	free(X);
	free(wsave);
	free(ifac);

	return 0;
}
