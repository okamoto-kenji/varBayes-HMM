/*
 *  vbHmmPcDataHandler.h
 *  File loader for VB-HMM-PC data.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

//// Assumed data format:
// File contents are binary array of unsigned short integer, big endian.
// 2 * (N+2) bytes, where N is number of data points
// Data point represents the photon counts.
// 1st 4bytes represents the size of time bin in unit of ns.

#include "vbHmmPc.h"

xnDataSet *readPcBinary( char*, FILE* );

//
