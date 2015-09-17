/*
 *  vbHmmGaussDataHandler.h
 *  File loader for VB-HMM-GAUSS data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

//// Assumed data format:
// A text file containing time series of real values.
// A data point per line.

#include "vbHmmGauss.h"

xnDataSet *readGaussText( char*, FILE* );
double *readDoubleArrayFromFile( char*, size_t*, FILE* );

//
