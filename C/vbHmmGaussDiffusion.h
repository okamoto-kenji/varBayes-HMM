/*
 *  vbHmmGaussDiffusion.h
 *  Model-specific core functions for VB-HMM-GAUSS-DIFFUSION.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

#ifndef VBHMMGAUSSDIFF_DEF
#define VBHMMGAUSSDIFF_DEF

#include "gVbHmm_Common.h"

// model-specific data
typedef struct _gaussDiffData {
    double *v;       // value
} gaussDiffData;

xnDataSet *newXnDataSet_gaussDiff( const char* );
void freeXnDataSet_gaussDiff( xnDataSet** );

// model-specific parameters
typedef struct gaussDiffParameters_ {
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double *avgDlt, *avgLnDlt;
    double *uAArr, *uBArr;
    double *aDlt, *bDlt;
} gaussDiffParameters;

void *newModelParameters_gaussDiff( xnDataSet*, int );
void freeModelParameters_gaussDiff( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _gaussDiffStats {
    double *Ni, *Ri, *Nii, **Nij;
} gaussDiffStats;

void *newModelStats_gaussDiff( xnDataSet*, globalVars*, indVars* );
void freeModelStats_gaussDiff( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
typedef struct _gaussDiffGlobalStats {
    double *NiR, *RiR, *NiiR, **NijR, *z1iR;
} gaussDiffGlobalStats;

void *newModelStatsG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
void freeModelStatsG_gaussDiff( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_gaussDiff( xnDataSet*, globalVars*, indVars* );
void initializeVbHmmG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_gaussDiff( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputGaussDiffResults( xnDataSet*, globalVars*, indVars*, FILE* );
void outputGaussDiffResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gaussDiff();
void setGFunctions_gaussDiff();

double pTilde_z1_gaussDiff( int, void* );
double pTilde_zn_zn1_gaussDiff( int, int, void* );
double pTilde_xn_zn_gaussDiff( xnDataSet*, size_t, int, void* );

void calcStatsVars_gaussDiff( xnDataSet*, globalVars*, indVars* );
void calcStatsVarsG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_gaussDiff( xnDataSet*, globalVars*, indVars* );
void maximizationG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_gaussDiff( xnDataSet*, globalVars*, indVars* );
double varLowerBoundG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_gaussDiff( xnDataSet*, globalVars*, indVars* );
void reorderParametersG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_gaussDiff( xnDataSet*, globalVars*, indVars*, FILE* );
void outputResultsG_gaussDiff( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//
