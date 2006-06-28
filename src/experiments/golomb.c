#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "include/seq.h"


int golomb_postulate_1(seq_t *seq)
{
	unsigned int i, ones;
	int res;

	for ( i = 0, ones = 0; i < seq->n; i++ )
		ones += SEQ(seq, i);

	fprintf(stdout, "Golomb's ramdomness postulate 1\n");
	fprintf(stdout, "Balance property:\n");
	fprintf(stdout, "--------------------------------\n");
	fprintf(stdout, "Number of 0's: %d\n", seq->n - ones);
	fprintf(stdout, "Number of 1's: %d\n", ones);

	switch ( seq->n - 2*ones ) {
		case -1:
		case 0:
		case 1:
			fprintf(stdout, "Result: Accepted\n\n");
			res = 0;
			break;
		default:
			fprintf(stdout, "Result: Rejected\n\n");
			res = 1;
	}
	
	return res;
}

static int add_run(unsigned int ***runs, unsigned int *size, 
								unsigned int type, unsigned int len)
{
	unsigned int j;
	unsigned int **pruns, *run;

	for ( j = 0, pruns = *runs; j < *size; j++, pruns++ )
		if ( **pruns == len) {
				(*(*pruns + 1 + type))++;
				break;
	}

	if ( j == *size ) {
						
		pruns = (unsigned int **)realloc(*runs, (*size+1)*sizeof(unsigned int *));
		if ( !pruns ) {
			fprintf(stderr, "Error: Memory allocation failed.\n");
			return -1;
		}

		*runs = pruns;

		run = *(pruns + (*size)++) = malloc(3*sizeof(unsigned int));
		if ( !run ) {
				fprintf(stderr, "Error: Memory allocation failed.\n");
				return -1;
		}

		memset(run, 0, 3*sizeof(unsigned int));
		*run = len;
		*(run + 1 + type) = 1;
	}

	return 0;
}

static void sort_runs(unsigned int **runs, unsigned int size)
{
	unsigned int i, nperm;
	unsigned int tmp[3];

	do {
		for ( i = 0, nperm = 0; i < (size-1); i++ ) {
				if ( **(runs + i) > **(runs + i + 1) ) {
					memcpy(tmp, *(runs + i + 1), 3*sizeof(unsigned int));
					memcpy(*(runs + i + 1), *(runs + i), 3*sizeof(unsigned int));
					memcpy(*(runs + i), tmp, 3*sizeof(unsigned int));
					nperm++;
				}
		}
	} while ( nperm > 0 );
}

static int check_runs(unsigned int **runs, unsigned int size, unsigned int nruns)
{
	unsigned int i;
	unsigned int **pruns;

	for (i = 1, pruns = runs; i <= size; i++, pruns++) {

		if ( **pruns != i )
			return 1;

		if ( (*(*pruns + 1) + *(*pruns + 2)) < 1/pow(2,i)*nruns )
			return 1;
	}

	return 0;
}

int golomb_postulate_2(seq_t *seq)
{
	int res;
	unsigned int **runs = NULL, **pruns;
	unsigned int size = 0;
	unsigned int type = 2;
	unsigned int i, b, len, nruns;

	fprintf(stdout, "Golomb's ramdomness postulate 2\n");
	fprintf(stdout, "Run property:\n");
	fprintf(stdout, "--------------------------------\n");

	type = SEQ(seq, 0);

	for ( i = 1, len = 1, nruns = 1, res = 0; i < seq->n; i++ ) {

		b = SEQ(seq, i);
		if ( b == type )
			len++;
		else {
			
			if ( add_run(&runs, &size, type, len) == -1 ) {
				res = -1;
				break;
			}

			type = b;
			len = 1;
			nruns++;
		} 
	} 

	if ( res == 0 ) {
		if ( add_run(&runs, &size, type, len) == -1 )
			res = -1;
		else {

			sort_runs(runs, size);

			fprintf(stdout, "Number of runs: %d\n", nruns);
			fprintf(stdout, "Length\t  Gap\tBlock\n");
			for ( i = 0, pruns = runs; i < size; i++, pruns++ )
				fprintf(stdout, "%5d\t%5d\t%5d\n", **pruns, *(*pruns+1), *(*pruns+2));

			if ( check_runs(runs, size, nruns) == 0 )
				fprintf(stdout, "Result: Accepted\n\n");
			else
				fprintf(stdout, "Result: Rejected\n\n");
		}
	}
	
	if ( runs ) {
		for ( i = 0, pruns = runs; i < size; i++, pruns++ )
			if ( *pruns )
				free(*pruns);
		free(runs);
	}

	return res;
}

static int autocorrelation(seq_t *seq, int t)
{
	int i, C;

	for ( i = 0, C = 0; i < seq->n; i++ ) 
		C +=  (2*PNSEQ(seq, i) - 1)*(2*PNSEQ(seq, i + t) - 1);

	return C;
}

int golomb_postulate_3(seq_t *seq)
{
	int C, K;
	unsigned int t;

	fprintf(stdout, "Golomb's ramdomness postulate 3\n");
	fprintf(stdout, "Autocorrelation property:\n");
	fprintf(stdout, "--------------------------------\n");

	K = autocorrelation(seq, 1);
	
	for ( t = 2; t < seq->n; t++ ) {
		C = autocorrelation(seq, t);
		if ( C != K ) {
				fprintf(stdout, "Result: Rejected\n\n");
				return 1;
		}
	}
			
	fprintf(stdout, "K = %d\n", K);
	fprintf(stdout, "Result: Accepted\n\n");

	return 0;
}

int main(int argc, char **argv)
{
	unsigned char buff[2] = {0x26, 0x5e};
	seq_t seq = {buff, 15, 2};
	int res, res2;

	fprintf(stdout, "S = ");
	print_seq(&seq, 0);
	
	res = golomb_postulate_1(&seq);

	res2 = golomb_postulate_2(&seq);
	if ( res2 == -1 )
		exit(EXIT_FAILURE);

	res += golomb_postulate_3(&seq) + res2;

	if ( res == 0 )
		fprintf(stdout, "Every golomb's postulate is satisfied, S is a pn-sequence.\n");
	else
		fprintf(stdout, "Golom's postulates are not satisfied, S is not a pn-sequence.\n");

	return 0;
}
