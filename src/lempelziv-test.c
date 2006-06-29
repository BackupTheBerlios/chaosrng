#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

#include "include/seq.h"
#include "include/lempelziv-test.h"


typedef struct bintree bintree_t;

struct bintree {
	bintree_t *nodes[2];
};


static void bintree_free(bintree_t *b)
{
	if ( b->nodes[0] ) {
		bintree_free(b->nodes[0]);
		free(b->nodes[0]);
	}
	if ( b->nodes[1] ) {
		bintree_free(b->nodes[1]);
		free(b->nodes[1]);
	}
}

int lempelziv_test(seq_t *seq, double *pvalue, void *param)
{
	unsigned int i, W;
	unsigned char b;
	bintree_t root, *curnode, *tmp;

	if ( seq->n < LEMPELZIV_TEST_LENGTH ) {
		fprintf(stderr, "Error[LempelZiv Test]: Sequence length too short\n");
		return -1;
	}
	else if ( seq->n > LEMPELZIV_TEST_LENGTH )
		fprintf(stderr, "Warning[LempelZiv Test]: %d bits discarded\n", seq->n - LEMPELZIV_TEST_LENGTH);
	
	root.nodes[0] = NULL;
	root.nodes[1] = NULL;

	for ( i = 0, W = 0, curnode = &root; i < LEMPELZIV_TEST_LENGTH; i++ ) {

		b = SEQ(seq, i);
		if ( !curnode->nodes[b] ) {

			tmp = (bintree_t *)malloc(sizeof(bintree_t));
			if ( !tmp ) {
				fprintf(stderr, "Error[LempelZiv Test]: Memory allocation failed\n");
				bintree_free(&root);
				return -2;
			}

			tmp->nodes[0] = NULL;
			tmp->nodes[1] = NULL;
			curnode->nodes[b] = tmp;
			curnode = &root;
			W++;
		}
		else
			curnode = curnode->nodes[b];
	}
	
	*pvalue = gsl_cdf_ugaussian_Q((69586.25 - ((double)W))/70.448718);

	return 0;
}
