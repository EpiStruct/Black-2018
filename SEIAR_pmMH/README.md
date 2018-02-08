# pmMH algorithm

This code implements the particle marginal Metropolis Hastings algorithm using an alive particle filter and importance sampling, as described in Section 8 of the paper. This is coded in C and requires the [GNU Scientific Library](https://www.gnu.org/software/gsl/).

## compiling

just run make.

## running the mcmc

Once compiled, the program is run with the following arguments

./mcmc_SEIAR particles thin runtime name algo

- particles - the number of particles to use.
- thin - thinning of the chain.
- runtime - length of time to run for in minutes.
- name - name of file to save samples and stats to.
- algo - 1 or 2. 1 uses importance sampling, 2 uses the alive sampling. 

For example: ./mcmc_SEIAR 50 10 15 some_output 1

Note that the other mcmc parameters and the data are hard coded in the file mcmc_SEIAR.c
There is very little checking done by this program so it will not detect incorrect parameters range / types or invalid file names, so be careful.

## implementation

There are three things to watch out for in this particular implementation. 
1. The events are indexed from 0 instead of 1, in contrast to the matlab code.
2. The order of the events has been rearranged as described in the discussion of the paper. This means that the observed event rate, which is always 0, is the last element in the vector of events. This allows for a little bit of optimisation. Saying that, this code is not particularly optimised. 
3. The variable names a and b for the original and modified event rates have been swapped. 