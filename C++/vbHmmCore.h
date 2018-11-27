/*
 *  vbHmmCore.h
 *  Common VB-HMM engine.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#ifndef __vbHmmCore__
#define __vbHmmCore__

#include "vbHmm_Class.h"
#include "vbHmmModel_Common.h"

using namespace std;

// manages the VB-HMM engine to choose likeliest model
int modelComparison(vbHmmCond *c, vbHmmData *d, vbHmmModel *m, fstream *logFS);
// main routine of VB-HMM
int vbHmm_Main( vbHmmCond *c, vbHmmData *d, vbHmmModel *m, fstream *logFS );
// Baum-Welch algorithm
void forwardBackward( vbHmmCond *c, vbHmmData *d, vbHmmModel *m );
// pick up state trajectory by maximum gamma
#ifdef OUTPUT_MAX_GAMMA
void maxGamma( vbHmmCond *c, vbHmmData *d, vbHmmModel *m );
#endif
// Viterbi algorithm
void maxSum( vbHmmCond *c, vbHmmData *d, vbHmmModel *m );

#endif /* defined(__vbHmmCore__) */

//
