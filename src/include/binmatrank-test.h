#ifndef HAVE_BINMATRANK_TEST_H
#define HAVE_BINMATRANK_TEST_H

#define BINMATRANK_M 32

#define BINMATRANK_N BINMATRANK_M*BINMATRANK_M

#define BINMATRANK_TEST_LENGTH 38*BINMATRANK_M*BINMATRANK_M

#define BINMATRANK_SIZE (int)((BINMATRANK_M*BINMATRANK_M)/8) + 1

int binmatrank_test(seq_t *seq, double *pvalue, void *param);

#endif /* !HAVE_BINMATRANK_TEST_H */
