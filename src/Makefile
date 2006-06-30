all : basic-test

CC = gcc

CFLAGS = -Wall -Wno-long-double -I/sw/include/gsl

LIBS =  -L/sw/lib -lfftw3 -lgsl -lgslcblas -lm

SRC = seq.c gen.c gen-urandom.c frequency-test.c\
	frequencywb-test.c\
	runs-test.c\
	longestrun-test.c\
	binmatrank-test.c\
	dft-test.c\
	notm-test.c\
	otm-test.c\
	maurer-test.c\
	lincomplex-test.c\
	lempelziv-test.c\
	serial-test.c\
	approxentropy-test.c\
	ranex-test.c\
	vranex-test.c\
	cumsum-test.c\
	basic-test.c\

DEP = Makefile

OBJ = $(SRC:.c=.o)

.SUFFIXES :
.SUFFIXES : .c .o
.c.o:
	$(CC) $(CFLAGS) -c $<

basic-test: $(OBJ) $(DEP)
	$(CC) $(CFLAGS) $(LIBS) $(OBJ) -o $@

clean :
	rm -f *.o
