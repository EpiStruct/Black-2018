//
//  SEIAR_mcmc.c
//
//  Created by Andrew Black on 26/10/17.
//  Copyright © 2017 Andrew Black. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "particle_mcmc.h"
#include "SEIAR_model.h"

int main(int argc, const char * argv[]) {
    
    if (argc != 6){
        
        printf("arguments: <part> <thin> <runtime (mins)> <output name> <algo 1/2>\n");
        exit(EXIT_FAILURE);
        
    }
    
    
    // config setup for the SEIAR model.
    Config *m;
    m = malloc(sizeof(Config));
    
    // dimension of the parameter space and number of events.
    m->dim = 5;
    m->events = 5;
    
    m->low = gsl_vector_alloc(m->dim);
    m->high = gsl_vector_alloc(m->dim);
    m->theta_init = gsl_vector_alloc(m->dim);
    
    // boundaries and initial point.
    m->low->data = (double[]){0.1, 0.1, 0.5, 0.5, 0.0};
    m->high->data = (double[]){15.0, 15.0, 15.0, 1.0, 1.0};
    m->theta_init->data = (double[]){8.0, 1.0, 1.0, 0.9, 0.5};
    
    // set the covariance matrix
    m->sigma = gsl_matrix_alloc(m->dim,m->dim);
    
    FILE * f = fopen ("SEIAR_cov_mat.txt", "r");
    if(f){
        gsl_matrix_fscanf(f,m->sigma);
        fclose(f);
     
    }
    else{
        printf("Cant read matrix\n");
        exit(EXIT_FAILURE);
    }
    
//    gsl_matrix_set_identity(m->sigma);
    m->scale = 0.5;
    
    m->part = atoi(argv[1]);
    m->thin = atoi(argv[2]);
    m->runTime = atof(argv[3]);
    
    switch (atoi(argv[5])) {
        case 1:
            m->likePtr = SEIAR_likelihood_is;
            printf("Exact importance sampler\n");
            break;
        case 2:
                m->likePtr = SEIAR_likelihood_alive;
            printf("Alive filter\n");
            break;
        default:
            printf("type must be 1 or 2");
            exit(EXIT_FAILURE);
            break;
    }
    
    
    asprintf(&m->out_file, "%s%s", argv[4], ".txt");
    asprintf(&m->stats_file, "%s_stats%s", argv[4], ".txt");
    
    printf("part = %d,\noutput file: %s\n",m->part,m->out_file);
    printf("thinning = %d,\n", m->thin);    
    
    // time-series of cases.
    int length = 15;

    // allocate struct to hold time-series data.
    m->ts = malloc(sizeof(data_s));
    m->ts->vec = gsl_vector_int_alloc(length);
    m->ts->vec->data = (int[]){1,0,3,11,27,40,24,15,15,4,4,0,1,0,1};
    m->ts->NF = 146;
    m->N = 164;
    
    
    m->max_samples = pow(10,7);
    
    // run the mcmc for given length of time.
    mcmc_time(m);

    return 0;
}



