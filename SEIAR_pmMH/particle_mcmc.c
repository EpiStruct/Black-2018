//
//  particle_mcmc.c
//  part_test
//
//  Created by Andrew Black on 26/10/17.
//  Copyright © 2017 Andrew Black. All rights reserved.
//

#include "particle_mcmc.h"

#include <math.h>
#include <gsl/gsl_linalg.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>



PartStruct *init_part(int part,int events, int parameters){
    // initialise the particle array and return a pointer to it.
    
    PartStruct *p;
    p = malloc(sizeof(PartStruct));
    p->part = part;
    p->events = events;
    p->Z = calloc( part * events, sizeof(p->Z));
    p->Z_tmp = calloc( part * events, sizeof(p->Z_tmp));
    
    p->par = gsl_vector_alloc(parameters);
    p->w = calloc(part,sizeof(p->w));
    
    return p;
}

void free_PartStruct(PartStruct *p){
    // free the particle array memory.
    
    free(p->Z);
    free(p->Z_tmp);
    free(p->w);
    gsl_vector_free(p->par);
    free(p);
    
}

void print_state(PartStruct *p){
    
    for (int i=0; i<p->part; i++) {
        
        for (int j=0; j< p->events; j++) {
            
            printf("%d ",p->Z[i*p->events+j]);
            
        }
        
        printf("w=%f \n",p->w[i]);
    }
    
}


double test_likelihood( Config* m, gsl_vector* theta){
    // given some parameters calculate the likelihood for our time series
    
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,mySeed);
    
    // setup the space for the particles.
    PartStruct *X = init_part(m->part,m->events,m->dim);
    X->N = m->N;
    
    double like = m->likePtr(r,X,theta,m->ts);

    free_PartStruct(X);
    gsl_rng_free(r);
    
    return like;
   
}






void resample_part(gsl_rng *r, PartStruct *X){
    // Systematic re-sampling.

    const int part = X->part;
    const int events = X->events;
    double cumsum[part];
    int tmp_part[part][events];
    
    // calculate the cumulative sum of the weights
    cumsum[0] = X->w[0];
    for (int i = 1; i < part; ++i)
    {
        cumsum[i] = X->w[i] + cumsum[i-1];
    }
    
    for (int i = 0; i < part; ++i)
    {
        cumsum[i] = cumsum[i]/cumsum[part-1];
    }
    
    
    unsigned int new_samp[part];
    
    double inc = 1.0/part;
    double u = gsl_rng_uniform(r)*inc ;
    int i = 0;
    for (int j = 0; j < part; ++j)
    {
        while(u > cumsum[i]){
            ++i;
        }
        new_samp[j] = i;
        
        for (int k=0; k<events; k++) {
            tmp_part[j][k] = X->Z[i*events+k];
        }
        
        u = u + inc;
    }
    
    for (int j = 0; j < part; j++)
    {
        // copy over the old particles.
        for (int i=0; i<events; i++) {
            X->Z[j*events+i] = tmp_part[new_samp[j]][i];
        }
        
    }
    
    
}

void multinomial_resample(gsl_rng *r, PartStruct *X){

    const int part = X->part;
    const int events = X->events;
    unsigned int n[part]; // holds the samples.

    gsl_ran_multinomial(r,part,part,X->w,n);

    // sample particles into the tmp array
    int count = 0;

    for (int k = 0; k < part; ++k)
    {
        for(int l = 0; l < n[k]; l++){

            for (int i=0; i<events; i++) {
                X->Z_tmp[count*events+i] = X->Z[k*events+i];
            }
            count ++;
        }
    }

    // copy back into the original array.
    for (int k = 0; k < part; ++k)
    {
        for (int i=0; i<events; i++) {
            X->Z[k*events+i] = X->Z_tmp[k*events+i];
        }

    }

}


int check_support(Config *m, gsl_vector *can){
    // check that candidate is in the support
    
    for (int i = 0; i < m->dim; ++i)
    {
        if (gsl_vector_get(can,i) <= gsl_vector_get(m->low,i))
        {
            return 0;
        }
        
        if (gsl_vector_get(can,i) >= gsl_vector_get(m->high,i))
        {
            return 0;
        }
        
    }
    
    return 1;
}



void fprintf_matrix(gsl_matrix *X, const char *f_name, int rows){
    // print a matrix to a file.
    // bum is the name of the file to output
    
    FILE * f = fopen (f_name, "w");
    
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < X->size2; ++j)
        {
            fprintf(f, "%f ",gsl_matrix_get(X,i,j));
        }
        fprintf(f, "\n");
    }
    
    fclose (f);
}




void mcmc_time(Config *m){
    

    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,mySeed);
    
    // setup the space for the particles.
    PartStruct *X = init_part(m->part,m->events,m->dim);
    X->N = m->N;
    
    const int dim = m->dim;
    
    gsl_vector *theta = gsl_vector_calloc(dim);
    gsl_vector *theta_can = gsl_vector_calloc(dim);
    
    
    // convert to the cholesky form for the mvn algorithm.
    gsl_matrix *sigma = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(sigma,m->sigma);
    gsl_matrix_scale(sigma, m->scale);
    gsl_linalg_cholesky_decomp1(sigma);
        
    // this matrix will hold the samples from the mcmc.
    gsl_matrix *output = gsl_matrix_calloc(m->max_samples,dim);
        
    double can_like = 0.0;
    
    // set initial parameters.
      gsl_vector_memcpy(theta,m->theta_init);

    for(int i=0; i<dim; i++){

        gsl_matrix_set(output,0,i,theta->data[i]);  
    }

    
    double like = m->likePtr(r,X,theta,m->ts);
    double curr_prior = m->prior(theta);

    int ar = 0, os = 0;
    
    time_t time_now = time(NULL);
    time_t start_time = time_now;
    time_t endTime = time_now + (m->runTime)*60;
    
    printf("started at : %s",ctime(&time_now));
    printf("should end at : %s",ctime(&endTime));
    fflush(stdout);
    
    int j = 1;
    int k = 0;
    
    while((time_now < endTime) && (j < m->max_samples))
    {
        
        // generate the candidate point.
        gsl_ran_multivariate_gaussian(r,theta,sigma,theta_can);
        
        // if theta_can is within the support then estimate likelihood.     
        if (check_support(m,theta_can))
        {
            os++;
            can_like = m->likePtr(r,X,theta_can,m->ts);
            double prior_can = m->prior(theta_can);

             // accept / reject step.
            if( log(gsl_rng_uniform(r)) < can_like - like + prior_can - curr_prior){
                
                gsl_vector_memcpy(theta,theta_can);
                like = can_like;
                curr_prior = prior_can;
                ar++;
            }
            
            // update the time, but only in the slow loop where we evaluate the likelihood.
            
            time(&time_now);
        }
        
        k ++;
        
        if(k > (m->thin)-1){

            gsl_matrix_set_row(output,j,theta);
            j ++;
            k = 0;
            
        }
  
    }
    
    
    printf("\nended at : %s",ctime(&time_now));
    printf("took %ld seconds",time_now-start_time);
    
    printf("\nSamples = %d, AR = %f, OS = %f \n\n",j,(double)ar/os,(double)os/j/m->thin);
    
    // write the output file
    fprintf_matrix(output,m->out_file,j);
    
    // write the stats file
    FILE * f = fopen (m->stats_file, "w");
    
    fprintf(f, "%ld %f",time_now-start_time,(double)ar/os);
    
    fclose (f);
 
    gsl_vector_free(theta);
    gsl_vector_free(theta_can);
    gsl_matrix_free(output);
    free_PartStruct(X);
}





