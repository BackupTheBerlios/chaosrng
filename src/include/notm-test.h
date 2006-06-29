#ifndef HAVE_NOTM_TEST_H
#define HAVE_NOTM_TEST_H

#define NOTM_TEST_LENGTH 1000000

typedef struct {
	unsigned int M;
	seq_t *seq;	
} notm_param_t;

int notm_test(seq_t *seq, double *pvalue, void *param);

#endif /* !HAVE_NOTM_TEST_H */
