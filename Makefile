CC=gcc
ACC_CC=pgcc

COMMON_CFLAGS=-O2 -g
BASE_CFLAGS=-Wall -pedantic -std=gnu99
LDFLAGS=-lcrypto

ACC_CFLAGS=-acc -Minfo=accel -ta=tesla:managed
#-noswitcherror -Mlarge_arrays -ta=tesla:cc60
C_DEPS=random.c

.PHONY: all clean verify

all: ca_seq

ca: ca.cu $(C_DEPS)
	nvcc $(COMMON_CFLAGS) $(LDFLAGS) $^ -o $@
ca_seq: ca_seq.c $(C_DEPS)
	$(CC) $(COMMON_CFLAGS) $(BASE_CFLAGS) $(LDFLAGS) $^ -o $@

ca_openacc: ca_openacc.c $(C_DEPS)
	$(ACC_CC) $(COMMON_CFLAGS) $(ACC_CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f ca_seq *.o