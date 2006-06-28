#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "include/seq.h"


seq_t *seq_new(unsigned int n)
{
	seq_t *seq;

	seq = (seq_t *)malloc(sizeof(seq_t));
	if ( !seq )
		return NULL;

	seq->n = n;
	seq->size = (int)(n/8) + 1;

	seq->buff = (unsigned char *)malloc(seq->size*sizeof(unsigned char));
	if ( !seq->buff ) {
		free(seq);
		return NULL;
	}

	memset(seq->buff, 0, seq->size*sizeof(unsigned char));
	
	return seq;
}

void seq_free(seq_t *seq) 
{
	if ( seq ) {

		if ( seq->buff )
			free(seq->buff);
		
		free(seq);
	}	
}

int seq_cmp(seq_t *seq1, seq_t *seq2, unsigned int offset, unsigned int len)
{
	int i;

	if ( seq1->n != len )
		return -1;

	for ( i = 0; i < len; i++ ) {
		if ( SEQ(seq1, i) != SEQ(seq2, (offset + i)%seq2->n) )
			return -1;
	}

	return 0;
}

int seq_cmp2(seq_t *seq1, seq_t *seq2, unsigned int offset1, unsigned int offset2, unsigned int len)
{
	int i;

	for ( i = 0; i < len; i++ ) {
		if ( SEQ(seq1, offset1 + i) != SEQ(seq2, offset2 + i) )
			return -1;
	}

	return 0;
}

void seq_init_with_uchar(seq_t *seq, unsigned char *tab, int n)
{
	int i,j,k;
	unsigned char tmp = 0;

	memset(seq->buff, 0, seq->size*sizeof(unsigned char));

	for ( i = 0, j = 0, k = 0; i < n; i++, k++ ) {

			if ( k == 8 ) {
					k = 0;
					seq->buff[j] = tmp;
					tmp = 0;
					j++;
			}

			tmp |= tab[i] << k; 
	}

	if ( k > 0 )
		seq->buff[j] = tmp;

	seq->n = n;
}

void seq_init_with_seq(seq_t *seq1, seq_t *seq2, unsigned int n, unsigned int offset)
{
	int i;

	seq1->n = n;

	for ( i = 0; i < n; i++ )
		setSEQ(seq1, i, SEQ(seq2, offset + i));
}

void seq_print(seq_t *seq, unsigned int npl)
{
	int i, j;

	for ( i = 0, j = 0; i < seq->n; i++) {

		fprintf(stdout, "%d ", SEQ(seq, i));
		j++;

		if ( npl != 0 && j == npl ) {
			fprintf(stdout, "\n");
			j = 0;
		}
	}

	if ( j  > 0 )
			fprintf(stdout, "\n");
}


void seq_add_double(seq_t *seq, unsigned int offset, double f, int k, int l)
{
	int i;
	unsigned long long d;
	double m;
	unsigned char c;

	m = frexp(f, &i);
	memcpy(&d, &m, sizeof(double));

	d = d >> 8;
	for ( i = 0; i < 3; i++ ) {
		c = d&0xff;
		d = d >> 8;
		memcpy(&(seq->buff[offset*3 + i]), &c, sizeof(unsigned char));
	}
}

void seq_add_iter(seq_t *seq, unsigned int niter, double iter, int k, int l) 
{
	unsigned int i, j, offset_seq, offset_tmp;

	seq_t tmp = { (unsigned char *)&iter, 64, 8};

	offset_seq = niter*l;
	offset_tmp = (int)(k/8);

	for ( i = 0, j = 7 - (k%8); i < l; i++, j-- ) {
		setSEQ(seq, offset_seq + i, SEQ(&tmp, (offset_tmp*8 + j)));		
		if ( j == 0 ) {
			j = 8;
			offset_tmp++;
		}
	}
}
