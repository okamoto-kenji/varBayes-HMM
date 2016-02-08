/*
 *  vbHmmGaussIntensityDataHandler.h
 *  File loader for VB-HMM-GAUSS-INTENSITY data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

//// Assumed data format:
// A text file containing time series of real values.
// A data point per line.

#include "vbHmmGaussIntensity.h"

xnDataSet *readGaussIntText( char*, FILE* );
double *readDoubleArrayFromFile( char*, size_t*, FILE* );

//
