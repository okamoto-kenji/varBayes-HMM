/*
 *  vbHmmTsDataHandler.h
 *  File loader for VB-HMM-TS data.
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
// File contents are binary array of unsigned long integer, big endian.
// 4 * (N+1) bytes, where N is number of data points
// Data point is lapsed time before the previous data point in unit of sampling ticks
// 1st 4bytes is the sampling frequency.
// data[1-] = time interval (ticks) from previous data point
// delta_t(n), in unit of time, is given as
//  delta_t(n) = data[n] / data[0]
// where n = 1...N.

#include "vbHmmTs.h"

xnDataSet *readTimeStampBinary( char*, FILE* );
void tsDataSetFromLongArray( xnDataSet*, unsigned long*, size_t, FILE* );

//
