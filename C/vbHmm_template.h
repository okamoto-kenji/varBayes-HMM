/*
 *  vbHmm_template.h
 *  Template for original model to be analyzed by VB-HMM.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

// Template for arbitraty model to be analyzed by VB-HMM.

#ifndef VBHMMTEMPLATE_DEF
#define VBHMMTEMPLATE_DEF

// Include common global VB-HMM engine.
#include "gVbHmm_Common.h"

// Define variables to describe model-specific data.
typedef struct _tempData {
} tempData;

xnDataSet *newXnDataSet__( const char* );
void freeXnDataSet__( xnDataSet** );

// model-specific parameters
typedef struct _tempParameters {
    // Define model specific parameters here.
    // For example, average and log average of parameters,
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    //   hyperparameter for prior distributions.
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
} tempParameters;

void *newModelParameters__( xnDataSet*, int );
void freeModelParameters__( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _tempStats {
    double **Nij, *Nii, *Ni;
} tempStats;

void *newModelStats__( xnDataSet*, globalVars*, indVars* );
void freeModelStats__( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
typedef struct _tempGlobalStats {
    double **NijR, *NiiR, *NiR, *z1iR;
} tempGlobalStats;

void *newModelStatsG__( xnDataBundle*, globalVars*, indVarBundle* );
void freeModelStatsG__( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm__( xnDataSet*, globalVars*, indVars* );
void initializeVbHmmG__( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars__( xnDataSet*, globalVars*, indVars* );

// model-specific output
void output__Results( xnDataSet*, globalVars*, indVars*, FILE* );
void output__ResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions__();
void setGFunctions__();

double pTilde_z1__( int, void* );
double pTilde_zn_zn1__( int, int, void* );
double pTilde_xn_zn__( xnDataSet*, size_t, int, void* );

void calcStatsVars__( xnDataSet*, globalVars*, indVars* );
void calcStatsVarsG__( xnDataBundle*, globalVars*, indVarBundle* );
void maximization__( xnDataSet*, globalVars*, indVars* );
void maximizationG__( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound__( xnDataSet*, globalVars*, indVars* );
double varLowerBoundG__( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters__( xnDataSet*, globalVars*, indVars* );
void reorderParametersG__( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults__( xnDataSet*, globalVars*, indVars*, FILE* );
void outputResultsG__( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//
