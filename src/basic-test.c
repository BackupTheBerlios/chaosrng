#include <stdio.h>
#include <string.h>
#include <math.h>

#include "include/seq.h"
#include "include/gen.h"
#include "include/gen-urandom.h"
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


static struct {
	unsigned int id;
	const char *name;
	int (*test)(seq_t*, double*, void*);
	unsigned int npvalue;
	unsigned int param;
} def_battery[] = {
	{0, "Frequency", &frequency_test, 1, 0},
	{1, "FrequencyWB", &frequencywb_test, 1, 1},
	{2, "Runs", &runs_test, 1, 0},
	{3, "LongestRun", &longestrun_test, 1, 0},
	{4, "BinMatRank", &binmatrank_test, 1, 0},
	{5, "DFT", &dft_test, 1, 0}, 
	{6, "NOTM", &notm_test, 1, 1}, 
	{7, "OTM", &otm_test, 1, 1}, 
	{8, "Maurer", &maurer_test, 1, 0},
	{9, "LempelZiv", &lempelziv_test, 1, 0},
	{10, "LinComplex", &lincomplex_test, 1, 1},
	{11, "Serial", &serial_test, 2, 1},
	{12, "ApproxEntropy", &approxentropy_test, 1, 1},
	{13, "CumSum", &cumsum_test, 1, 0},
	{14, "RanEx", &ranex_test, 1, 0},
	{15, "VRanEx", &vranex_test, 1, 0},
};

typedef struct {
	unsigned int id;
	void *param;
} battery_t;


static int battery_test_unidim(seq_t *seq, int ntest, void *param)
{
	int i, ret;
	double pvalue[2];

	ret = def_battery[ntest].test(seq, pvalue, param);
	if ( ret ==  0 ) {
		for ( i = 0; i < def_battery[ntest].npvalue; i++ )
			fprintf(stdout, "%-15f ", pvalue[i]);
	}
	else
		fprintf(stdout, "%-15d ", ret);

	fflush(stdout);

	return ret;	
}

static void battery_test(seq_t **gen, unsigned int dim, unsigned char mask, battery_t *bat, unsigned int n) 
{
	int i, j, ret;

	for ( j = 0; j < dim; j++ ) {

		if ( mask & (1<<j) ) {
			for ( i = 0; i < n; i++ ) {
				ret = battery_test_unidim(*(gen + j), bat[i].id, bat[i].param);
			}

			fprintf(stdout, "\n");
			fflush(stdout);
		}
	}
}

static void battery_print_header(battery_t *bat, unsigned int n)
{
	unsigned int i;

	for ( i = 0; i < n; i++ )
		fprintf(stdout, "%-15s ", def_battery[bat[i].id].name);
	fprintf(stdout, "\n");
	fflush(stdout);
}

int main(int argc, char **argv)
{
	/*
	seq_t **gen;
	//double param[2] = {2.2, -0.61};
	double param;
	double seed[2] = {0.5, 0.6};
	unsigned int M1 = 11000, M10 = 2250, M11 = 16, M12 = 16;
	battery_t bat[] = {
		{0, NULL},
		{1, (void *)M1},
		{2, NULL},
		{3, NULL},
		{4, NULL},
		{5, NULL},
		{8, NULL},
		{9, NULL},
		//{10, (void *)M10},
		//{11, (void *)M11},
		//{12, (void *)M12},
		{13, NULL}
		//{14, NULL},
		//{15, NULL}
	};


	fprintf(stdout, "%-11s", "lambda");
	battery_print_header(bat, 8);

	//param = 0.99989;
	//param = 0.883131;
	param = 0.7;

	gen = gen1_new(10000, 250000, &param, seed, 12, 4);

	do {
		fprintf(stdout, "%.8f ", param);
		battery_test(gen, 2, 1, bat, 8);

		param += 0.00000001;


		if ( param >= 1.0 ) 		
			break;

		gen1_init(gen, 10000, 250000, &param, seed, 12, 4);

	} while ( 1 );

	gen1_free(gen);
	*/

	/*
	seq_t *seq;
	unsigned char tab[20] = {0,1,0,1,1,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1};
	double pvalue;

	seq = seq_new(20);
	seq_init_with_uchar(seq, tab, 20);

	maurer_test(seq, &pvalue, NULL);
	fprintf(stdout, "%f\n", pvalue);
	*/

	/*
	seq_t *seq;
	double pvalue;

	seq = seq_urandom_new(125000);

	dft_test(seq, &pvalue, NULL);
	fprintf(stdout, "%f\n", pvalue);
	*/

	return 0;
}
