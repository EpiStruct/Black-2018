CFLAGS=-I/usr/local/include -Wall -L/usr/local/lib -std=c99 -O3 -DHAVE_INLINE -D_GNU_SOURCE=1
LDFLAGS= -lgsl -lgslcblas -lm 
CC=gcc

mcmc_SEIAR: mcmc_SEIAR.c SEIAR_model.c particle_mcmc.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@
