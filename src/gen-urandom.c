#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <fcntl.h>

#include "include/seq.h"
#include "include/gen-urandom.h"


seq_t *seq_urandom_new(unsigned int n)
{
	int fd, ret;
	seq_t *seq;

	seq = seq_new(n*8);
	if ( !seq )
		return NULL;

	fd = open("/dev/urandom", O_RDONLY, 0);
	if ( fd < 0 ) {
		seq_free(seq);
		return NULL;
	}

	ret = read(fd, seq->buff, n);
	if ( ret < 0 ) {
		seq_free(seq);
		seq = NULL;
		fprintf(stderr, "Error[URANDOM]\n");
	}

	close(fd);

	return seq;
}
