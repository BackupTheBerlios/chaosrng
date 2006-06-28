#include <stdio.h>
#include <string.h>
#include <math.h>

#include "include/seq.h"
#include "include/gen.h"
#include "include/frequency-test.h"
#include "include/frequencywb-test.h"
#include "include/runs-test.h"
#include "include/longestrun-test.h"
#include "include/binmatrank-test.h"
#include "include/dft-test.h"
#include "include/notm-test.h"
#include "include/otm-test.h"
#include "include/maurer-test.h"
#include "include/lempelziv-test.h"
#include "include/lincomplex-test.h"
#include "include/serial-test.h"
#include "include/approxentropy-test.h"
#include "include/cumsum-test.h"
#include "include/ranex-test.h"
#include "include/vranex-test.h"


typedef struct {
	const char *name;
	int (*test)(seq_t*, double*, double*);
	unsigned int npvalue;
	unsigned int nparam;
} battery_t;

static const battery_t battery[] = {
	{"Frequency", &frequency_test, 1, 0},
	//{"FrequencyWB", &frequencywb_test, 1, 1},
	//{"Runs", &runs_test, 1, 0},
	//{"LongestRun", &longestrun_test, 1, 1},
	//{"BinMatRank", &binmatrank_test, 1, 0},
	//{"DFT", &dft_test, 1, 0}, 
	//{"Maurer", &maurer_test, 1, 1},
	{"LempelZiv", &lempelziv_test, 1, 0}, // 10^6
	//{"LinComplex", &lincomplex_test, 1, 1}, // 10^6
	//{"Serial", &serial_test, 2, 1},
	//{"ApproxEntropy", &approxentropy_test, 1, 1},
	//{&cumsum_test, 1, 0},
	//{&ranex_test, 1, 1}, // 10^6
	//{&vranex_test, 1, 1}, // 10^6
	{NULL, 0, 0}
};


static int battery_test_unidim(seq_t *seq, int ntest, double *param)
{
	int i, ret;
	double pvalue[2];

	ret = battery[ntest].test(seq, pvalue, param);
	if ( ret ==  0 ) {
		for ( i = 0; i < battery[ntest].npvalue; i++ )
			fprintf(stdout, " %f", pvalue[i]);
	}
	else
		fprintf(stdout, " %d", ret);

	fflush(stdout);

	return ret;	
}


static void battery_test(seq_t **gen, unsigned int dim) 
{
	int i, j, ret;
	double param = 1;

	for ( i = 0; battery[i].test; i++ ) {
		
		fprintf(stdout, "%-10s\t", battery[i].name);
		fflush(stdout);

		for ( j = 0; j < dim; j++ ) {
			ret = battery_test_unidim(*(gen + j), i, &param);
		}

		fprintf(stdout, "\n");
	}
}

int main(int argc, char **argv)
{
	seq_t **gen;
	double param[2] = {2.2, -0.61};
	double seed[2] = {0.7, 0.5};

	gen = gen3_new(100, 100000, param, seed, 12, 10);

	battery_test(gen, 2);

	gen3_free(gen);

	return 0;
}
