/*
 *  rand.c
 *  Random number generation
 *  Reference: 
 *    W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
 *    "Numerical Receipes", 3rd Edition, 
 *    Cambridge University Press, 2006
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include <math.h>
#include <time.h>
#include "rand.h"

gsl_rng *rng_ = NULL;    // rundom number generator

void initRan(){
    if( rng_ == NULL ){

        // choose one generatr from the following:
        // three types of fast & simulation quality generators
        rng_ = gsl_rng_alloc( gsl_rng_taus );
//        rng_ = gsl_rng_alloc( gsl_rng_gfsr4 );
//        rng_ = gsl_rng_alloc( gsl_rng_mt19937 );
        // strongest generators
//        rng_ = gsl_rng_alloc( gsl_rng_ranlxs0 );
//        rng_ = gsl_rng_alloc( gsl_rng_ranlxs1 );
//        rng_ = gsl_rng_alloc( gsl_rng_ranlxs2 );

        if( rng_ == NULL ){
            fprintf( stderr, "Rondome Number Generator couldn't be generated.\n" );
            exit(1);
        }
        gsl_rng_set( rng_, time(NULL) );
    }
    
}

void freeRan(){
    gsl_rng_free( rng_ );
    rng_ = NULL;
}

double ranDoub(){
    // returns value in range of [0.0, 1.0).
    return gsl_rng_uniform( rng_ );
}

double ranInRange( double x, double y ){
    // returns value in range of (x, y).
    return gsl_ran_flat( rng_, x, y );
}

double enoise( double x ){
    // returns value in range of (-x, x).
    return gsl_ran_flat( rng_, -x, x );
}

double interval( double I ){
    // returns value from distribution p(x) = I exp(- I x).
    return gsl_ran_exponential( rng_, 1.0/I );
}

//
