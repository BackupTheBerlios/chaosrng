#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>
#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/dft-test.h"


int dft_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int i, N;
	double *X, *pX, d, T;

	fftw_complex *out;
	fftw_plan p;

	if ( seq->n < DFT_TEST_LENGTH ) {
		fprintf(stderr, "Error[DFT Test]: Sequence length too short\n");
		return -1;
	}

	X = (double *)malloc(seq->n*sizeof(double));
	if ( !X ) {
		fprintf(stderr, "Error[DFT Test]: Memory allocation failed\n");
		return -2;
	}

	memset(X, 0, seq->n*sizeof(double));

	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(seq->n/2 + 1));
	if ( !out ) {
		free(X);
		fprintf(stderr, "Error[DFT Test]: Memory allocation failed\n");
		return -2;
	}

	for ( i = 0, pX = X; i < seq->n; i++, pX++ )
		 X[i] = (double)(2*SEQ(seq, i) - 1);

	p = fftw_plan_dft_r2c_1d(seq->n, X, out, FFTW_ESTIMATE);

	fftw_execute(p);

	T = 0.95*sqrt((double)(3*seq->n));

	for ( i = 0, N = 0; i < seq->n/2; i++ ) {
		fprintf(stderr, "%f %f\n", out[i][0], out[i][1]);
		if ( sqrt(pow(out[i][0], 2.0) + pow(out[i][1], 2.0)) < T)
			N++;
	}

	d = (((double)N)-(0.95*seq->n)/2.0)/sqrt(seq->n*0.95*0.05/2.0);

	*pvalue = 2.0*gsl_cdf_ugaussian_Q(fabs(d));

	free(X);
	fftw_destroy_plan(p);
	fftw_free(out);

	return 0;
}
