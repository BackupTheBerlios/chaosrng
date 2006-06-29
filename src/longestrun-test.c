#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/longestrun-test.h"


static unsigned int longestrun(seq_t *seq, unsigned int i, unsigned int M)
{
	unsigned int r, j = 0,lr = 0;

	while ( j < M ) {

		r = 0;
		
		while ( j < M && SEQ(seq, i*M + j) == 1 ) {	
			r++;
			j++;
		}

		if ( lr < r )
			lr = r;

		j++;
	}
	
	return lr;
}

int longestrun_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int N;
	typedef struct {
		unsigned int M;
		unsigned int nbClass;
		unsigned int firstClass;
		unsigned int nMin;
		double *p;
	} longestrun_data_t;

      	static double pM8[4] = {0.21484375,0.3671875,0.23046875,0.1875};
	static double pM128[6] = {0.1174035788,0.242955959,0.249363483,0.17517706, 0.102701071, 0.112398847};
	static double pM10000[7] = {0.0882,0.2092,0.2483,0.1933,0.1208,0.0675,0.0727};

	static longestrun_data_t data[3] = {
		{8, 4, 1, 128, pM8},
		{128, 6, 4, 6272, pM128},
		{10000, 7, 10, 750000, pM10000}	
	};
	
	int indM;
	unsigned int i, lr;
	unsigned int *class;
	
	double X;


	if ( seq->n < data[0].nMin ) {
		fprintf(stderr, "Error[LongestRun Test]: Sequence length too short\n");
		return -1;
	}
	else if ( seq->n < data[1].nMin )
		indM = 0;
	else if ( seq->n < data[2].nMin )
		indM = 1;
	else
		indM = 2;

	class = (unsigned int *)malloc(data[indM].nbClass*sizeof(unsigned int));
	if ( !class ) {
		fprintf(stderr, "Error[LongestRun Test]: Memory allocation failed\n");
		return -2;
	}

	memset(class, 0, data[indM].nbClass*sizeof(unsigned int));


	N = (int)(seq->n/data[indM].M); 

	for ( i = 0; i < N; i++ ) {

		lr = longestrun(seq, i, data[indM].M);

		if ( lr <= data[indM].firstClass )
			class[0]++;
		else if ( lr >= (data[indM].firstClass + data[indM].nbClass - 1) )
			class[data[indM].nbClass - 1]++;
		else
			class[lr - data[indM].firstClass]++;
		
	}

	for ( i = 0, X = 0; i < data[indM].nbClass; i++ )
		X += pow(((double)class[i]) - ((double)N)*data[indM].p[i], 2.0)/(((double)N)*data[indM].p[i]);

	*pvalue = gsl_cdf_chisq_Q(X, (data[indM].nbClass - 1));

	free(class);

	return 0;
}
