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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "particle_mcmc.h"
#include "SEIAR_model.h"

double prior_eval(gsl_vector* theta){
    // evaluate the prior at theta.
    double R0_p = 1.0; //log(gsl_ran_gamma_pdf(theta->data[0],10.0,2.0/10.0));
    double sig_p = log(gsl_ran_gamma_pdf(theta->data[1],10.0,1.0/10.0));
    double gam_p = log(gsl_ran_gamma_pdf(theta->data[2],10.0,1.0/10.0));

    return R0_p + sig_p +  gam_p;
}

int main(int argc, const char * argv[]) {
    
    if (argc != 7){
        
        printf("arguments: <part> <thin> <runtime (mins)> <output name> <algo 1/2> <ts 1,2,3 or 4>\n");
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
    m->high->data = (double[]){8.0, 15.0, 15.0, 1.0, 1.0};
    m->theta_init->data = (double[]){2.0, 1.0, 1.0, 0.9, 0.5};
    
    // set the covariance matrix
    m->sigma = gsl_matrix_alloc(m->dim,m->dim);
    
    FILE * f = fopen ("cov_mat/cov_p.txt", "r");
    if(f){
        gsl_matrix_fscanf(f,m->sigma);
        fclose(f);
     
    }
    else{
        printf("Cant read matrix\n");
        exit(EXIT_FAILURE);
    }
    
//    gsl_matrix_set_identity(m->sigma);
    m->scale = 0.3;
    
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
    

    m->prior = prior_eval;

    asprintf(&m->out_file, "%s%s", argv[4], ".txt");
    asprintf(&m->stats_file, "%s_stats%s", argv[4], ".txt");
    
    printf("part = %d,\noutput file: %s\n",m->part,m->out_file);
    printf("thinning = %d,\n", m->thin);    
    
    // time-series
    m->ts = malloc(sizeof(data_s));


    switch (atoi(argv[6])) {
        case 1:
            printf("Small boat\n");

            m->ts->vec = gsl_vector_int_alloc(18);
            m->ts->vec->data = (int[]) {0, 0, 1, 8, 6, 14, 12, 15, 15, 14, 10, 5, 3, 8, 4, 2, 2, 2};
            m->ts->NF = 121;
            m->N = 150;
            break;
        case 2:
            printf("Medium boat\n");
            // med boat.
            m->ts->vec = gsl_vector_int_alloc(26);
            m->ts->vec->data = (int[]) {0, 0, 1, 8, 11, 10, 21, 17, 24, 32, 31, 29, 21, 19, 21, 7, 9, 8, 6, 3, 2, 4, 1,
                                        1, 1, 1};
            m->ts->NF = 288;
            m->N = 350;
            break;
        case 3:
            printf("Large boat\n");
            m->ts->vec = gsl_vector_int_alloc(26);
            m->ts->vec->data = (int[]) {0, 1, 1, 3, 4, 6, 3, 5, 11, 23, 22, 27, 18, 34, 39, 37, 28, 30, 28, 21, 12, 16,
                                        7, 4, 1, 2};
            m->ts->NF = 383;
            m->N = 500;
            break;
        case 4:
            printf("A bigger boat\n");
            m->ts->vec = gsl_vector_int_alloc(40);
            m->ts->vec->data = (int[]) {0, 0, 0, 0, 7, 2, 8, 7, 9, 17, 13, 38, 54, 48, 66, 82, 78, 81, 68, 54, 49, 37, 17, 14, 13, 5, 9, 3, 1, 1, 2, 1, 1, 1, 1, 2, 0, 0, 0, 1};
            m->ts->NF = 790;
            m->N = 1000;
            break;
        default:
            printf("Which time series to use? 1,2,3 or 4");
            exit(EXIT_FAILURE);
            break;
    }


    m->max_samples = pow(10,7);
    
    // run the mcmc for given length of time.
    mcmc_time(m);

    return 0;
}



