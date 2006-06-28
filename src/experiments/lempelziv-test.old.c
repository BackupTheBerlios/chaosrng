#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/lempelziv-test.h"


int lempelziv_test(seq_t *seq, double *pvalue, double *param)
{
	unsigned int i, j, k, W;
	seq_t **dico, **tmp;

	if ( seq->n == LEMPELZIV_TEST_LENGTH ) {
		fprintf(stderr, "Error[LempelZiv Test]: Sequence length too short\n");
		return -1;
	}

	dico = (seq_t **)malloc(sizeof(seq_t *));
	if ( !dico ) {
		fprintf(stderr, "Error[LempelZiv Test]: Memory allocation failed\n");
		return -2;
	}

	*dico = seq_new(1);
	if ( !*dico ) {
		fprintf(stderr, "Error[LempelZiv Test]: Memory allocation failed\n");
		free(dico);
		return -2;
	}

	setSEQ(*dico, 0, SEQ(seq, 0));
	(*dico)->n = 1;

	for ( i = 1, j = 1, W = 1; i < seq->n; i++ ) {
		
		for ( k = 0; k < W; k++ ) {
			if ( seq_cmp(*(dico+k), seq, j, i - j + 1) == 0 )
				break;
		}

		if ( k == W ) {

			tmp = (seq_t **)realloc(dico, (W+1)*sizeof(seq_t *));
			if ( !tmp ) {
				fprintf(stderr, "Error[LempelZiv Test]: Memory allocation failed\n");
				for ( i = 0; i < W; i++ )
					seq_free(*(dico + i));
				free(dico);
				return -2;
			}

			
			dico = tmp;

			*(dico + W) = seq_new(i - j + 1);
			if ( !*(dico + W) ) {
				fprintf(stderr, "Error[LempelZiv Test]: Memory allocation failed\n");
				for ( i = 0; i < W; i++ )
					seq_free(*(dico + i));
				free(dico);
				return -2;
			}
			
			seq_init_with_seq(*(dico + W), seq, i - j + 1, j);
			seq_print(*(dico + W), 100000);
			W++;
			j = i + 1;
		}
	}

	for ( i = 0; i < W; i ++ )
		seq_free(*(dico + i));
	free(dico);
	
	*pvalue = gsl_cdf_ugaussian_Q((69586.25 - ((double)W))/70.448718);

	return 0;
}
