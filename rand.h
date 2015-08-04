/*
 *  rand.h
 *  Random number generation
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#ifndef DEF_RAND
#define DEF_RAND

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void initRan();
void freeRan();
double ranInRange( double, double );
double enoise( double );
double interval( double );

#endif

//
