/*
 *  rand.h
 *  Random number generation
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
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
