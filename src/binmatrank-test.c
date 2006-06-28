#include <stdio.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/binmatrank-test.h"


static void row_xor(seq_t *mat, int n1, int n2, unsigned int m) 
{
	int i;

	for ( i = 0; i < m; i++ )
		setSEQ(mat, n2*m + i, SEQ(mat, n2*m + i)^SEQ(mat, n1*m + i));
}

static void row_swap(seq_t *mat, int n1, int n2, unsigned int m) 
{
	int i, tmp;

	for ( i = 0; i < m; i++ ) {
		tmp = SEQ(mat, n2*m + i);
		setSEQ(mat, n2*m + i, SEQ(mat, n1*m + i));
		setSEQ(mat, n1*m + i, tmp);
	}
}

static void forward_row_operations(seq_t *mat, unsigned int m) 
{
	int i, j;

	for( i = 0; i < m - 1	; i++) {
		if ( !SEQ(mat, i*(m + 1)) ) {
			for ( j = i + 1; j < m; j++ ) {
				if ( SEQ(mat, j*m +i) ) {
					row_swap(mat, i, j, m);
					break;
				}
			}
			for ( j++; j < m; j++) {
				if ( SEQ(mat, j*m + i) )
					row_xor(mat, i, j, m);
			}
		}
		else {
			for ( j = i + 1; j < m; j++ ) {
				if ( SEQ(mat, j*m + i) )
					row_xor(mat, i, j, m);
			}
		}
	}
}

static void backward_row_operations(seq_t *mat, unsigned int m) 
{
	int i, j, r;

	for( r = 0, i = m - 1; i > 0; i--) {
		if ( !SEQ(mat, i*(m+1)) ) {
			for ( j = i - 1; j > -1; j-- ) {
				if ( SEQ(mat, j*m + i) ) {
					row_swap(mat, i, j, m);
					r++;
					break;
				}
			}
			for ( j--; j > -1; j--) {
				if ( SEQ(mat, j*m + i) )
					row_xor(mat, i, j, m);
			}
		}
		else {
			r++;
			for ( j = i - 1; j > -1; j-- ) {
				if ( SEQ(mat, j*m + i) )
					row_xor(mat, i, j, m);
			}
		}
	}
}

static unsigned int rank(seq_t *mat, unsigned int m)
{
	int i, r;

	forward_row_operations(mat, m);
	backward_row_operations(mat, m);

	for ( i = 0, r = 0; i < m; r += SEQ(mat, i*(m + 1)), i++);

	return r;
}

int binmatrank_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int offset, R, R1, R2, N;
	double stat;

	unsigned char buff[BINMATRANK_SIZE];
	seq_t mat;

	mat.buff = buff;
	mat.size = BINMATRANK_SIZE;

	if ( seq->n < BINMATRANK_TEST_LENGTH ) {
		fprintf(stderr, "Error[BinMatRank Test]: Sequence length too short\n");
		return -1;
	}

	for ( offset = 0, R = 0, R1 = 0, R2 = 0, N = 0; (offset + BINMATRANK_N) < seq->n; offset += BINMATRANK_N ) {
		
		N++;

		seq_init_with_seq(&mat, seq, BINMATRANK_N, offset);
		
		R = rank(&mat, BINMATRANK_M);

		if ( R == BINMATRANK_M )
			R1++;
		else if ( R == BINMATRANK_M - 1 )
			R2++;
	}

	stat = pow(R1 - 0.2888*N, 2.0)/(0.2888*N) + pow(R2 - 0.5776*N, 2.0)/(0.5776*N) +pow(N - (R2 + R1) - 0.1336*N, 2.0)/(0.1336*N); 

	*pvalue = gsl_cdf_chisq_Q(stat, 2);

	return 0;
}
