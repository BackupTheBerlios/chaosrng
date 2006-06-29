#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/lincomplex-test.h"


static int BerlekampMassey(seq_t *seq, unsigned int offset, unsigned int len)
{
	int L = 0, N = 0, m = -1, d, i, j;
	unsigned char cbuff[(int)(len/8)], bbuff[(int)(len/8)], tbuff[(int)(len/8)];
	seq_t cseq, bseq, tseq;

	cseq.buff = cbuff;
	cseq.n = len;
	cseq.size = (int)(len/8);

	bseq.buff = bbuff;
	bseq.n = len;
	bseq.size = (int)(len/8);

	tseq.buff = tbuff;
	tseq.n = len;
	tseq.size = (int)(len/8);

	setSEQ(&cseq, 0, 1);
	setSEQ(&bseq, 0, 1);

	while ( N < len ) {
		d = SEQ(seq, offset + N);
		for ( i = 1; i <= L; i++ )
			d ^= SEQ(&cseq,i)&SEQ(seq, offset + N-i);
		if ( d == 1 ) {
			seq_init_with_seq(&tseq, &cseq, len, 0);
			for ( j = 0; (j+N-m) < len; j++ )
				setSEQ(&cseq, j+N-m,  SEQ(&cseq, j+N-m)^SEQ(&bseq, j));
			if ( L <= (N>>1) ) {
				L = N + 1 - L;
				m = N;
				seq_init_with_seq(&bseq, &tseq, len, 0);
			}
		}
		N++;
	}

	return L;
}

int lincomplex_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int M = (unsigned int)param;
	unsigned int i, j, L, N = (int)(seq->n/M);
	int  V[7] = { 0, 0, 0, 0, 0, 0, 0 };
	double pi[7] = { 0.01047, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833 };
	double mu, T, X;
			
	if ( seq->n < LINCOMPLEX_TEST_LENGTH ) {
		fprintf(stderr, "Error[LinComplex Test]: Sequence length too short\n");
		return -1;
	}

	if ( M < 500 || M > 5000 ) {
		fprintf(stderr, "Error[LinComplex Test]: Invalid parameter\n");
		return -1;
	}

	mu = ((double)M)/2.0 + (9.0 + pow(-1.0, M + 1))/36.0 - (((double)M)/3.0 + 2.0/9.0)/pow(2.0, M);

	for ( i = 0; i < N; i++ ) {
		L = BerlekampMassey(seq, i*M, M);
		T = pow(-1, M)*((double)(L) - mu) + 2.0/9.0;
		if ( T <= -2.5 )
			V[0]++;
		else if ( T > 2.5 )
			V[6]++;
		else {
			for ( j = 1; j < 6; j++ ) {
				if ( (-2.5*((double)j) + (double)(j-1) < T) && (-2.5*((double)j) + (double)(j) >= T) ) {
					V[j]++;
					break;
				}
			}
		}
	}

	for ( i = 0, X = 0.0; i < 7; i++ )
		X += pow(((double)V[i]) - ((double)N)*pi[i], 2)/(((double)N)*pi[i]);

	*pvalue = gsl_cdf_chisq_Q(X, 6);

	return 0;
}
