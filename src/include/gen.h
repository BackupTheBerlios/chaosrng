#ifndef HAVE_GEN_H
#define HAVE_GEN_H


#define GEN_INIT(vgen) \
void gen##vgen##_init(seq_t **seq, unsigned int offsetNiter, unsigned int Niter, double *param, double *seed, unsigned int k, unsigned int l)

#define GEN_NEW(vgen) \
seq_t **gen##vgen##_new(unsigned int offsetNiter, unsigned int Niter, double *param, double *seed, unsigned int k, unsigned int l)

#define GEN_FREE(vgen) \
void gen##vgen##_free(seq_t **seq)


/*
 * Double Logistic map
 */
GEN_INIT(1);
GEN_NEW(1);
GEN_FREE(1);

/*
 * Two-dimensional non invertible map
 */
GEN_INIT(2);
GEN_NEW(2);
GEN_FREE(2);

/*
 * Cubic type map
 */
GEN_INIT(3);
GEN_NEW(3);
GEN_FREE(3);


/*
 * 
 */
GEN_INIT(4);
GEN_NEW(4);
GEN_FREE(4);

/*
 * Order two DPCM system
 */
GEN_INIT(5);
GEN_NEW(5);
GEN_FREE(5);

/*
 * Opto electronic system
 */
GEN_INIT(6);
GEN_NEW(6);
GEN_FREE(6);

/*
 * Order three DPCM system
 */
GEN_INIT(7);
GEN_NEW(7);
GEN_FREE(7);

/*
 * Triple Logistic map
 */
GEN_INIT(8);
GEN_NEW(8);
GEN_FREE(8);


#endif /* !HAVE_GEN_H */
