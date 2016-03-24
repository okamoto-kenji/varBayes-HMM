/*
 *  vbHmmTs.h
 *  Model-specific core functions for VB-HMM-TS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */


#ifndef VBHMMTS_DEF
#define VBHMMTS_DEF

#include "vbHmm_Common.h"
//#include "gVbHmm_Common.h"

// model-specific data
typedef struct _tsData {
    double T;                       // total length of data (may not be always necessary)
    double *dt;
    double *time;
} tsData;

xnDataSet *newXnDataSet_ts( const char* );
void freeXnDataSet_ts( xnDataSet** );

// model-specific parameters
typedef struct tsParameters_ {
    double *uPiArr, sumUPi;
    double **uKMat, *sumUKArr;
    double *aIArr, *bIArr;
    double *avgPi, *avgLnPi, **avgK, **avgLnK, *avgLnKI, *avgI, *avgLnI;
} tsParameters;

void *newModelParameters_ts( xnDataSet*, int );
void freeModelParameters_ts( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _tsStats {
    double *Ni, *Nii, *Nij, **Mij, *Ti;
} tsStats;

void *newModelStats_ts( xnDataSet*, globalVars*, indVars* );
void freeModelStats_ts( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
//typedef struct _tsGlobalStats {
//    double **NiR, *NiiR, *NijR, **MijR, *TiR, *z1iR;
//} tsGlobalStats;

//void *newModelStatsG_ts( xnDataBundle*, globalVars*, indVarBundle* );
//void freeModelStatsG_ts( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_ts( xnDataSet*, globalVars*, indVars* );
//void initializeVbHmmG_ts( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_ts( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputTsResults( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputTsResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_ts();
//void setGFunctions_ts();

double pTilde_z1_ts( int, void* );
double pTilde_zn_zn1_ts( int, int, void* );
double pTilde_xn_zn_ts( xnDataSet*, size_t, int, void* );

void calcStatsVars_ts( xnDataSet*, globalVars*, indVars* );
//void calcStatsVarsG_ts( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_ts( xnDataSet*, globalVars*, indVars* );
//void maximizationG_ts( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_ts( xnDataSet*, globalVars*, indVars* );
//double varLowerBoundG_ts( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_ts( xnDataSet*, globalVars*, indVars* );
//void reorderParametersG_ts( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_ts( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputResultsG_ts( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//
