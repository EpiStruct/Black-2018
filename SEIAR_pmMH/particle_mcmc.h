//
//  particle_mcmc.h
//  part_test
//
//  Created by Andrew Black on 26/10/17.
//  Copyright © 2017 Andrew Black. All rights reserved.
//

#ifndef particle_mcmc_h
#define particle_mcmc_h

#include <stdio.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"


// Struct to hold particle stuff.
typedef struct
{
    int N;
    int part;
    int events;
    gsl_vector *par;
    int *Z;
    int *Z_tmp;
    double *w;
} PartStruct;


PartStruct* init_part(int,int,int);
void resample_part(gsl_rng*, PartStruct* );
void free_PartStruct(PartStruct *p);

// time-series data struct
typedef struct 
{
    gsl_vector_int* vec;
    int NF; // final size.
} data_s;



// configuration data
typedef struct
{
    int part;
    int dim;
    int events;

    double runTime;
    int max_samples;
    int thin;
    
    char* out_file;
    char* stats_file;
    char* cov_file;
    
    // data on time series and size of population.
    data_s* ts;
    int N;
 
    gsl_vector* theta_init;
    gsl_vector *low, *high;
    gsl_matrix* sigma;
    double scale;
    
    // pointer to the likelihood function
    double (*likePtr)(gsl_rng*, PartStruct *, gsl_vector* , data_s*);
    
    
} Config;

void mcmc_time(Config*);

void print_state(PartStruct*);

double test_likelihood( Config* m, gsl_vector* theta);


#endif /* particle_mcmc_h */
