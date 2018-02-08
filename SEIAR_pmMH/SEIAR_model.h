//
//  SEIAR_model.h
//  part_test
//
//  Created by Andrew Black on 25/10/17.
//  Copyright © 2017 Andrew Black. All rights reserved.
//

#ifndef SEIAR_model_h
#define SEIAR_model_h

#include "particle_mcmc.h"

double SEIAR_likelihood_alive(gsl_rng*, PartStruct*, gsl_vector*, data_s*);
double SEIAR_likelihood_is(gsl_rng*, PartStruct*, gsl_vector*, data_s*);


#endif /* SEIAR_model_h */
