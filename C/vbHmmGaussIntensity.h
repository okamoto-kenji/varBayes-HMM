/*
 *  vbHmmGaussIntensity.h
 *  Model-specific core functions for VB-HMM-GAUSS-INTENSITY.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

#ifndef VBHMMGAUSSINT_DEF
#define VBHMMGAUSSINT_DEF

#include "gVbHmm_Common.h"

// model-specific data
typedef struct _gaussIntData {
    double *v;       // value
} gaussIntData;

xnDataSet *newXnDataSet_gaussInt( const char* );
void freeXnDataSet_gaussInt( xnDataSet** );

// model-specific parameters
typedef struct gaussIntParameters_ {
//    int sNo;
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double avgMu, avgLm, avgLnLm;
    double uBt, uMu, uA, uB;
    double btMu, aLm, bLm, mu0;
} gaussIntParameters;

void *newModelParameters_gaussInt( xnDataSet*, int );
void freeModelParameters_gaussInt( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _gaussIntStats {
    double N0, N0i, Nx;
    double **Nij, *Nii, *Ni;
} gaussIntStats;

void *newModelStats_gaussInt( xnDataSet*, globalVars*, indVars* );
void freeModelStats_gaussInt( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
typedef struct _gaussIntGlobalStats {
    double N0R, N0iR, NxR;
    double **NijR, *NiiR, *NiR, *z1iR;
} gaussIntGlobalStats;

void *newModelStatsG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
void freeModelStatsG_gaussInt( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_gaussInt( xnDataSet*, globalVars*, indVars* );
void initializeVbHmmG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_gaussInt( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputGaussIntResults( xnDataSet*, globalVars*, indVars*, FILE* );
void outputGaussIntResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gaussInt();
void setGFunctions_gaussInt();


double pTilde_z1_gaussInt( int, void* );
double pTilde_zn_zn1_gaussInt( int, int, void* );
double pTilde_xn_zn_gaussInt( xnDataSet*, size_t, int, void* );

void calcStatsVars_gaussInt( xnDataSet*, globalVars*, indVars* );
void calcStatsVarsG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_gaussInt( xnDataSet*, globalVars*, indVars* );
void maximizationG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_gaussInt( xnDataSet*, globalVars*, indVars* );
double varLowerBoundG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_gaussInt( xnDataSet*, globalVars*, indVars* );
void reorderParametersG_gaussInt( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_gaussInt( xnDataSet*, globalVars*, indVars*, FILE* );
void outputResultsG_gaussInt( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//

