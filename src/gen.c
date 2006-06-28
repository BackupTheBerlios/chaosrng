#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "include/seq.h"
#include "include/gen.h"


static seq_t **gen_new(unsigned int dim, unsigned int size)
{
	unsigned int i, j;
	seq_t **ret, *tmp;

	ret = malloc(dim*sizeof(seq_t *));
	if ( !ret )
		return NULL;

	for ( i = 0; i < dim; i++ ) {

		tmp = seq_new(size);
		if ( !tmp ) {

			for ( j = 0; j < i; j++ )
				seq_free(*(ret+j));

			free(ret);

			return NULL;
		}

		*(ret + i) = tmp;
	}
	
	return ret;
}


#define GEN_NEW_F(vgen, vdim) \
GEN_NEW(vgen)															\
{																\
	seq_t **seq;														\
	seq = gen_new(vdim, Niter*l);												\
	if ( !seq ) {														\
		fprintf(stderr, "Error[Gen]: Memory allocation failed\n");							\
		return NULL;													\
	}															\
	gen##vgen##_init(seq, offsetNiter, Niter, param, seed, k, l);								\
	return seq;														\
}

#define GEN_NEW_F2(vgen) GEN_NEW_F(vgen, 2)

#define GEN_NEW_F3(vgen) GEN_NEW_F(vgen, 3)


static void gen_free(seq_t **seq, unsigned int dim)
{
	unsigned int i;

	for ( i = 0; i < dim; i++ )
		seq_free(*(seq+i));

	free(seq);
}

#define GEN_FREE_F(vgen, vdim) \
GEN_FREE(vgen) { gen_free(seq, vdim); }

#define GEN_FREE_F2(vgen) GEN_FREE_F(vgen, 2)

#define GEN_FREE_F3(vgen) GEN_FREE_F(vgen, 3)

/*
 * Double Logistic map
 */

GEN_INIT(1)
{
	unsigned int i;
	double x0, y0, x, y, a;

	x0 = *seed;
	y0 = *(seed + 1);

	a = *param;

	for ( i = 0; i < offsetNiter; i++ ) {
		x = (1.0 - a)*x0 + 4.0*a*y0*(1.0 - y0);
		y = (1.0 - a)*y0 + 4.0*a*x0*(1.0 - x0);
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = (1.0 - a)*x0 + 4.0*a*y0*(1.0 - y0);
		y = (1.0 - a)*y0 + 4.0*a*x0*(1.0 - x0);

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(1)
GEN_FREE_F2(1)

/*
 * Two-dimensional non invertible map
 */

GEN_INIT(2)
{
	unsigned int i;
	double x0, y0, x, y, a, b;

	x0 = *seed;
	y0 = *(seed + 1);

	a = *param;
	b = *(param + 1);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = a*x0 + y0;
		y = pow(x0, 2.0) + b;
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = a*x0 + y0;
		y = pow(x0, 2.0) + b;

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(2)
GEN_FREE_F2(2)


/*
 * Cubic type map
 */

GEN_INIT(3)
{
	unsigned int i;
	double x0, y0, x, y, a, b;

	x0 = *seed;
	y0 = *(seed + 1);

	a = *param;
	b = *(param + 1);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = y0;
		y = a*(-1.0*pow(x0, 3.0) + x0) + b*(-1.0*pow(y0, 3.0) + y0);
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = y0;
		y = a*(-1.0*pow(x0, 3.0) + x0) + b*(-1.0*pow(y0, 3.0) + y0);

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(3)
GEN_FREE_F2(3)


/*
 * 
 */

GEN_INIT(4)
{
	unsigned int i;
	double x0, y0, x, y, a, b;

	x0 = *seed;
	y0 = *(seed + 1);

	a = *param;
	b = *(param + 1);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = y0;
		y = a*pow(sin(x0 + b), 2.0);
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = y0;
		y = a*pow(sin(x0 + b), 2.0);

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(4)
GEN_FREE_F2(4)


/*
 * Order two DPCM system
 */

#define Q(t) tanh(p*t)

GEN_INIT(5)
{
	unsigned int i;
	double x0, y0, x, y, a1, a2, p, s;

	x0 = *seed;
	y0 = *(seed + 1);

	a1 = *param;
	a2 = *(param + 1);
	p = *(param + 2);
	s = *(param + 3);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = y0;
		y = a1*(y0 + Q(s - y0)) + a2*(x0 + Q(s - x0));
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = y0;
		y = a1*(y0 + Q(s - y0)) + a2*(x0 + Q(s - x0));

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(5)
GEN_FREE_F2(5)


/*
 * Opto electronic system
 */

GEN_INIT(6)
{
	unsigned int i;
	double x0, y0, x, y, a, b;

	x0 = *seed;
	y0 = *(seed + 1);

	a = *param;
	b = *(param + 1);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = y0;
		y = a*pow(sin(b + x0 - y0), 2.0);
		x0 = x;
		y0 = y;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = y0;
		y = a*pow(sin(b + x0 - y0), 2.0);

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		

		x0 = x;
		y0 = y;
	}
}

GEN_NEW_F2(6)
GEN_FREE_F2(6)


/*
 * Order three DPCM system
 */

GEN_INIT(7)
{
	unsigned int i;
	double x0, y0, z0, x, y, z, a, b, c, p, s;

	x0 = *seed;
	y0 = *(seed + 1);
	z0 = *(seed + 2);

	a = *param;
	b = *(param + 1);
	c = *(param + 2);
	p = *(param + 3);
	s = *(param + 4);

	for ( i = 0; i < offsetNiter; i++ ) {
		x = y0;
		y = z0;
		y = b*(z0 + Q(s - z0)) + a*(y0 + Q(s - y0)) + c*(x0 + Q(s - x0));
		x0 = x;
		y0 = y;
		z0 = z;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = y0;
		y = z0;
		y = b*(z0 + Q(s - z0)) + a*(y0 + Q(s - y0)) + c*(x0 + Q(s - x0));

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		
		seq_add_iter(*(seq+2), i, z, k, l);		

		x0 = x;
		y0 = y;
		z0 = z;
	}
}

GEN_NEW_F3(7)
GEN_FREE_F3(7)


/*
 * Triple Logistic map
 */

GEN_INIT(8)
{
	unsigned int i;
	double x0, y0, z0, x, y, z, landa;

	x0 = *seed;
	y0 = *(seed + 1);
	z0 = *(seed + 2);

	landa = *param;

	for ( i = 0; i < offsetNiter; i++ ) {
		x = landa*(3.0*y0 + 1.0)*x0*(1.0 - x0);
		y = landa*(3.0*z0 + 1.0)*y0*(1.0 - y0);
		z = landa*(3.0*x0 + 1.0)*z0*(1.0 - z0);
		x0 = x;
		y0 = y;
		z0 = z;
	}

	for ( i = 0; i < Niter; i++ ) {

		x = landa*(3.0*y0 + 1.0)*x0*(1.0 - x0);
		y = landa*(3.0*z0 + 1.0)*y0*(1.0 - y0);
		z = landa*(3.0*x0 + 1.0)*z0*(1.0 - z0);

		seq_add_iter(*seq, i, x, k, l);		
		seq_add_iter(*(seq+1), i, y, k, l);		
		seq_add_iter(*(seq+2), i, z, k, l);		

		x0 = x;
		y0 = y;
		z0 = z;
	}
}

GEN_NEW_F3(8)
GEN_FREE_F3(8)
