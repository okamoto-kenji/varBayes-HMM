/*
 *  vbHmmGauss.h
 *  Model-specific core functions for VB-HMM-GAUSS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#ifndef VBHMMGAUSS_DEF
#define VBHMMGAUSS_DEF

#include "gVbHmm_Common.h"

// model-specific data
typedef struct _gaussData {
    double *v;       // value
} gaussData;

xnDataSet *newXnDataSet_gauss( const char* );
void freeXnDataSet_gauss( xnDataSet** );

// model-specific parameters
typedef struct _gaussParameters {
//    int sNo;
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double *avgMu, *avgLm, *avgLnLm;
    double *uBtArr, *uMuArr, *uAArr, *uBArr;
    double *btMu, *aLm, *bLm, *mu0;
} gaussParameters;

void *newModelParameters_gauss( xnDataSet*, int );
void freeModelParameters_gauss( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _gaussStats {
    double **Nij, *Nii, *Ni, *barX, *NiSi;
} gaussStats;

void *newModelStats_gauss( xnDataSet*, globalVars*, indVars* );
void freeModelStats_gauss( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
typedef struct _gaussGlobalStats {
    double **NijR, *NiiR, *NiR, *barXR, *NiSiR, *z1iR;
} gaussGlobalStats;

void *newModelStatsG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
void freeModelStatsG_gauss( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_gauss( xnDataSet*, globalVars*, indVars* );
void initializeVbHmmG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_gauss( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputGaussResults( xnDataSet*, globalVars*, indVars*, FILE* );
void outputGaussResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gauss();
void setGFunctions_gauss();

double pTilde_z1_gauss( int, void* );
double pTilde_zn_zn1_gauss( int, int, void* );
double pTilde_xn_zn_gauss( xnDataSet*, size_t, int, void* );

void calcStatsVars_gauss( xnDataSet*, globalVars*, indVars* );
void calcStatsVarsG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_gauss( xnDataSet*, globalVars*, indVars* );
void maximizationG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_gauss( xnDataSet*, globalVars*, indVars* );
double varLowerBoundG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_gauss( xnDataSet*, globalVars*, indVars* );
void reorderParametersG_gauss( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_gauss( xnDataSet*, globalVars*, indVars*, FILE* );
void outputResultsG_gauss( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//

