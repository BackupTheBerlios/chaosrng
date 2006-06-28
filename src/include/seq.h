#ifndef HAVE_SEQ_H
#define HAVE_SEQ_H


typedef struct {
	unsigned char *buff;
	unsigned int n;
	unsigned int size;
} seq_t;


#define SEQ(vs, vn) ((((vs)->buff[(int)((vn)/8)]) >> ((vn)%8))&0x01)

#define setSEQ(vs, vn, vv) do {\
		if ( (((vs)->buff[(int)((vn)/8)]&(0x01<<((vn)%8)))>>((vn)%8))!=(vv) ) {\
				if ( SEQ(vs, vn) == 0 )\
					(vs)->buff[(int)((vn)/8)] |= (0x01<<((vn)%8));\
				else\
					(vs)->buff[(int)((vn)/8)] ^= (0x01<<((vn)%8));\
		}\
	} while ( 0 )



seq_t *seq_new(unsigned int n);

void seq_free(seq_t *seq);

int seq_cmp(seq_t *seq1, seq_t *seq2, unsigned int offset, unsigned int len);

int seq_cmp2(seq_t *seq1, seq_t *seq2, unsigned int offset1, unsigned int offset2, unsigned int len);

void seq_init_with_uchar(seq_t *seq, unsigned char *tab, int n);

void seq_init_with_seq(seq_t *seq1, seq_t *seq2, unsigned int n, unsigned int offset);

void seq_print(seq_t *seq, unsigned int npl);

void seq_add_double(seq_t *seq, unsigned int offset, double f, int k, int l);

void seq_add_iter(seq_t *seq, unsigned int niter, double iter, int k, int l);


#endif /* !HAVE_SEQ_H */
