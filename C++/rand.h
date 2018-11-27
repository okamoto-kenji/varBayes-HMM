/*
 *  rand.h
 *  for vbHMM
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Last modified on 2018.11.26
 */

#ifndef __vbHmm__rand__
#define __vbHmm__rand__

#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


void initRan();
void freeRan();
double ranInRange(double, double);
double enoise(double);
double gnoise(double);
int randomInteger(int, int);

#endif /* defined(__vbHmm__rand__) */

//

