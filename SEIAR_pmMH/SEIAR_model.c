//
//  SEIAR_model.c
//  part_test
//
//  Created by Andrew Black on 25/10/17.
//  Copyright © 2017 Andrew Black. All rights reserved.
//

#include "SEIAR_model.h"

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>



void SEIAR_is(gsl_rng *r, PartStruct *X, int NR, int NF){
    // simulate the state forward 1 day using importance sampling.
    
    double b[5], a[5];
    double b0, a0;
    double t;
   	int Z[5];
    double orig, is;
    double lambda;
    double r2;
    
    double t_next[NR+2]; // holds the times of the detection events.
    int type_next[NR+2];
    double L_imp;
    double t_dash;
    
   	// derefernce pointers
   	const int N = X->N;
    const double gamma = X->par->data[2];
    const double sigma = X->par->data[1];
    const double p = X->par->data[3];
    const double betaP = (X->par->data[0])/(N-1);
    const double betaS = (X->par->data[4])/(N-1);
    
    double L_init = -gsl_sf_lnfact(NR);

    for (int i = 0; i < X->part; ++i)
    {
        
        Z[0] = X->Z[i*5];
        Z[1] = X->Z[i*5+1];
        Z[2] = X->Z[i*5+2];
        Z[3] = X->Z[i*5+3];
        Z[4] = X->Z[i*5+4];
        
        L_imp = L_init;
        
        // generate the order statistics over a day.
        double acm = 0.0;
        
        for (int j = 0; j<NR; j++){
            acm = acm - log(gsl_rng_uniform_pos(r));
            t_next[NR-j-1] = acm;
        }
        acm = acm - log(gsl_rng_uniform_pos(r));
        
        // normalise the times;
        for (int j = 0; j<NR; j++){
            t_next[j] = t_next[j]/acm;
            type_next[j] = 4;
        }
        
        t_next[NR] = 0.0;
        t_next[NR+1] = 0.0;
        type_next[NR] = 8; // this doesn't exist so will thrown an error is done.
        type_next[NR+1] = 8;

        int next = NR-1;
        
        t = 0;
        
        while(next >= 0){

            // calculate the rates of the original process.
            int E = (Z[0]-Z[1]-Z[2]);
            
            b[0] = (N-Z[0])*(betaP*(Z[1]-Z[4])+betaS*(Z[4]-Z[3]))*(Z[1]-Z[3]);
            b[1] = sigma*p*E;
            b[2] = sigma*(1-p)*E;
            b[3] = gamma*(Z[4]-Z[3]);
            b[4] = gamma*(Z[1]-Z[4]);
            
            b0 = b[0]+b[1]+b[2]+b[3]+b[4];
            
            // decide if we need to force an extra event or not.
            if(type_next[next]==4 && Z[1]-Z[4]==0){
                // next event is a 4, but Ip=0.
                
                if(E==0){
                    // E=0, so need to put in a type 0 first.
                    
                    lambda = b[0];
                    
                    double tE = -log(1-gsl_rng_uniform_pos(r)*(1-exp(-lambda*(t_next[next]-t))))/lambda;
                    double is_contr = log(lambda) - lambda*tE - log(1-exp(-lambda*(t_next[next]-t)));
                    L_imp = L_imp - is_contr;
                    next = next +1;
                    type_next[next] = 0;
                    t_next[next] = t + tE;
                }
                else{
                    // E > 0, so put in a type 1
                    lambda = b[1];
                    
                    double tE = -log(1-gsl_rng_uniform_pos(r)*(1-exp(-lambda*(t_next[next]-t))))/lambda;
                    double is_contr = log(lambda) - lambda*tE - log(1-exp(-lambda*(t_next[next]-t)));
                    L_imp = L_imp - is_contr;
                    next = next +1;
                    type_next[next] = 1;
                    t_next[next] = t + tE;
                }
            }
            else if(type_next[next]==1 && E==0){
                // next event is type 1, but E=0, so force type 0.
                
                lambda = b[0];
                
                double tE = -log(1-gsl_rng_uniform_pos(r)*(1-exp(-lambda*(t_next[next]-t))))/lambda;
                double is_contr = log(lambda) - lambda*tE - log(1-exp(-lambda*(t_next[next]-t)));
                
                L_imp = L_imp - is_contr;
                next = next +1;
                type_next[next] = 0;
                t_next[next] = t + tE;
            }

            
            // set the modified rates.
            
            a[0] = b[0];
            a[1] = b[1];
            a[2] = b[2];
            a[3] = b[3];
            a[4] = 0.0;

            
            if(type_next[next] == 0){
                a[0] = 0.0;
            }
            
            
            if(type_next[next] == 1){
                a[1] = 0.0;
            }
            
            if(Z[2] == N - NF){ a[2]=0.0;}
            
            
            if(Z[0]-Z[2]-Z[3] == 1){
                a[3] = 0.0;
                a[2] = 0.0;
            }
            
            if(Z[1] == NF){a[1]=0.0;}
            
            if(p==1 && Z[0]== NF){
                a[0]=0.0;
            }
            
            a0 = a[0]+a[1]+a[2]+a[3];
           
            
            t_dash = gsl_ran_exponential(r,1/a0);
            
            if(t_dash < t_next[next]-t){
                // put the event in
                
                r2 = gsl_rng_uniform(r);
                int index = 0;
                double sum = a[0];
                while(sum <= r2*a0)
                {
                    index ++;
                    sum += a[index];
                }
                
                t = t+t_dash;
                Z[index] = Z[index] +1;
                
                orig = log(b[index]) - b0*t_dash;
                is = log(a[index]) - a0*t_dash;
                
                L_imp = L_imp + orig - is;
                
            }
            else{
                // next forced event is implemented.
                int ne = type_next[next];
                
                orig = log(b[ne])-b0*(t_next[next]-t);
                is = -a0*(t_next[next]-t);
                L_imp = L_imp + orig - is;
                
                t = t_next[next];
                Z[ne]++;
                next--;
            }

            
            
        }
        
        
        
        while(t<1){
            
            // calculate the rates of the original process.
            int E = (Z[0]-Z[1]-Z[2]);
            
            b[0] = (N-Z[0])*(betaP*(Z[1]-Z[4])+betaS*(Z[4]-Z[3]))*(Z[1]-Z[3]);
            b[1] = sigma*p*E;
            b[2] = sigma*(1-p)*E;
            b[3] = gamma*(Z[4]-Z[3]);
            b[4] = gamma*(Z[1]-Z[4]);
            
            b0 = b[0]+b[1]+b[2]+b[3]+b[4];
            
            
            // set the modified rates.
            
            a[0] = b[0];
            a[1] = b[1];
            a[2] = b[2];
            a[3] = b[3];
            a[4] = 0.0;
            
           
            
            if(Z[2] == N - NF){ a[2]=0.0;}
            
            
            if(Z[0]-Z[2]-Z[3] == 1){
                a[3] = 0.0;
                a[2] = 0.0;
            }
            
            if(Z[1] == NF){a[1]=0.0;}
            
            if(p==1 && Z[0]== NF){
                a[0]=0.0;
            }
            
            a0 = a[0]+a[1]+a[2]+a[3];
            
            
            t_dash = gsl_ran_exponential(r,1/a0);
            
            if (t+t_dash > 1)
            {
                orig = -b0*(1-t);
                is = -a0*(1-t);
                
                L_imp = L_imp + orig - is;
                break;
            }
            else{
                
                r2 = gsl_rng_uniform(r);
                int index = 0;
                double sum = a[0];
                while(sum <= r2*a0)
                {
                    index ++;
                    sum += a[index];
                }
                
                t = t+t_dash;
                Z[index]++;
                
                is = log(a[index]) - a0*t_dash;
                orig = log(b[index]) - b0*t_dash;
                
                L_imp = L_imp + orig - is;
                
            }
   
          
        }

        
        // update the original array
        X->Z[i*5]   = Z[0];
        X->Z[i*5+1] = Z[1];
        X->Z[i*5+2] = Z[2];
        X->Z[i*5+3] = Z[3];
        X->Z[i*5+4] = Z[4];
        
        X->w[i] = exp(L_imp);
        
    }
    
}


double SEIAR_alive(gsl_rng *r, PartStruct *X, int NR, int NF){
    // Alive sampler.
    
    double b[5];
    double b0;
    double t;
   	int Z[5];
    int Z_obs;
    
    
   	// dereference pointers
   	const int N = X->N;
    const double gamma = X->par->data[2];
    const double sigma = X->par->data[1];
    const double p = X->par->data[3];
    const double betaP = (X->par->data[0])/(N-1);
    const double betaS = (X->par->data[4])/(N-1);
    const int part = X->part;
    const int Z2_max = N-NF;
    
    // max number of trials before we stop because prob is too low.
    int K = 100000;

    int j = 0;
    int n = 0;
    
   	while(j < part +1)
    {
        
        // sample a random particle
        unsigned long int i = gsl_rng_uniform_int(r,part);
        
        n++;
        
        Z[0] = X->Z[i*5];
        Z[1] = X->Z[i*5+1];
        Z[2] = X->Z[i*5+2];
        Z[3] = X->Z[i*5+3];
        Z[4] = X->Z[i*5+4];
        
        Z_obs = Z[4] + NR;

        
        t = 0;
        
        while( Z[0]-Z[2]-Z[3] > 0 && Z[4] <= Z_obs && Z[2] <= Z2_max && Z[1] <= NF){
            
            int E = (Z[0]-Z[1]-Z[2]);
            
            b[0] = (N-Z[0])*(betaP*(Z[1]-Z[4])+betaS*(Z[4]-Z[3]))*(Z[1]-Z[3]);
            b[1] = sigma*p*E;
            b[2] = sigma*(1-p)*E;
            b[3] = gamma*(Z[4]-Z[3]);
            b[4] = gamma*(Z[1]-Z[4]);
            
            b0 = b[0]+b[1]+b[2]+b[3]+b[4];
            
            // exponential time to next event.
            t = t + gsl_ran_exponential(r,1/b0);
            
            if (t > 1) break;
            
            double r2 = gsl_rng_uniform(r);
            int index = 0;
            double sum = b[0];
            while(sum <= r2*b0)
            {
                index ++;
                sum += b[index];
            }
            
            Z[index] ++;
            
        }
        
        
        if (Z[4] == Z_obs &&  Z[0]-Z[2]-Z[3] > 0 && Z[2] <= Z2_max && Z[1] <= NF){
            
            // add found particle to the tmp array.
            if(j < part){
                
                // store the particle
                X->Z_tmp[j*5]   = Z[0];
                X->Z_tmp[j*5+1] = Z[1];
                X->Z_tmp[j*5+2] = Z[2];
                X->Z_tmp[j*5+3] = Z[3];
                X->Z_tmp[j*5+4] = Z[4];
           
            }
            j++;
            
        }
        
        if (n > K){
            return 0;
        }
        
    }
    
    // put the tmp particles back into the original array
    for (int i = 0; i < part; ++i)
    {
        
        X->Z[i*5]   = X->Z_tmp[i*5];
        X->Z[i*5+1] = X->Z_tmp[i*5+1];
        X->Z[i*5+2] = X->Z_tmp[i*5+2];
        X->Z[i*5+3] = X->Z_tmp[i*5+3];
        X->Z[i*5+4] = X->Z_tmp[i*5+4];
    }
    
    return (double)part / (n-1);
    
}




double SEIAR_likelihood_alive(gsl_rng *r, PartStruct *X, gsl_vector *theta, data_s *ts){
    // calculate the marginal likelihood using the alive filter.

    
    // calculate the parameters      
    X->par->data[0] = theta->data[0]*theta->data[4] / theta->data[2];   //betaS
    X->par->data[1] = 1.0/theta->data[1];                               // 1/sigma
    X->par->data[2] = 1.0/theta->data[2];                               // 1/gamma
    X->par->data[3] = theta->data[3];                                   // p
    X->par->data[4] = theta->data[0]*(1-theta->data[4]) / theta->data[2]; //betaP
    
    size_t length = ts->vec->size;
    
    double LL[length];

    for (int j=0; j<length; j++) {
        LL[j] = 0.0;
    }

    double LL_tot;
    
    int part = X->part;
    
    // set initial conditions.
    for (int i = 0; i < part; i++){
        
        X->Z[i*5]   = 1;
        X->Z[i*5+1] = 1;
        X->Z[i*5+2] = 0;
        X->Z[i*5+3] = 0;
        X->Z[i*5+4] = 0;
    }
    
    
    for (int i = 0; i < length; i++)
    {
        // simulate forward a day.
        double like = SEIAR_alive(r, X, ts->vec->data[i], ts->NF);
        
        if(like > 0){
            
            LL[i] = log(like);
        }
        
        else{
            return -INFINITY;
        }
        
       
    }
    
    LL_tot = 0;
    for (int j=0; j<length; j++) {
        LL_tot = LL_tot + LL[j];
    }
    
    return LL_tot;
}





double SEIAR_likelihood_is(gsl_rng *r, PartStruct *X, gsl_vector *theta, data_s *ts){
    // calculate the marginal likelihood using importance sampling particle filter
    
    X->par->data[0] = theta->data[0]*theta->data[4] / theta->data[2];     //betaS
    X->par->data[1] = 1.0/theta->data[1];                                 // 1/sigma
    X->par->data[2] = 1.0/theta->data[2];                                 // 1/gamma
    X->par->data[3] = theta->data[3];                                     // p
    X->par->data[4] = theta->data[0]*(1-theta->data[4]) / theta->data[2]; //betaP

    size_t length = ts->vec->size;
    
    double LL[length];

    for (int j=0; j<length; j++) {
        LL[j] = 0.0;
    }

    int part = X->part;
    
    // set initial conditions.
    for (int i = 0; i < part; i++){
        
        X->Z[i*5]   = 1;
        X->Z[i*5+1] = 1;
        X->Z[i*5+2] = 0;
        X->Z[i*5+3] = 0;
        X->Z[i*5+4] = 0;
    }
    
    for (int i = 0; i < length-1; i++)
    {
        // simulate forward a day.
        SEIAR_is(r, X, ts->vec->data[i],ts->NF);
        
        // calculate the log-like contribution.
        LL[i] = log(gsl_stats_mean(X->w,1,part));

        resample_part(r,X);
        
    }
    
    SEIAR_is(r, X, ts->vec->data[length-1],ts->NF);

    LL[length-1] = log(gsl_stats_mean(X->w,1,part));

    
    double LL_tot = 0;
    for (int j=0; j<length; j++) {
        LL_tot = LL_tot + LL[j];
    }
    
    return LL_tot;
}


