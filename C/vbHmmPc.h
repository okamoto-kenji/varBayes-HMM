/*
 *  vbHmmPc.h
 *  Model-specific core functions for VB-HMM-PC.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */


#ifndef VBHMMPC_DEF
#define VBHMMPC_DEF

#include "vbHmm_Common.h"
//#include "gVbHmm_Common.h"

// model-specific data
typedef struct _pcData {
    double binSize;
    unsigned int *counts;
} pcData;

xnDataSet *newXnDataSet_pc( const char* );
void freeXnDataSet_pc( xnDataSet** );

// model-specific parameters
typedef struct pcParameters_ {
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *aIArr, *bIArr;
    double *avgPi, *avgLnPi, **avgA, **avgLnA, *avgI, *avgLnI;
} pcParameters;

void *newModelParameters_pc( xnDataSet*, int );
void freeModelParameters_pc( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _pcStats {
    double *Ni, **Nij, *Ci, *Mi;
} pcStats;

void *newModelStats_pc( xnDataSet*, globalVars*, indVars* );
void freeModelStats_pc( void**, xnDataSet*, globalVars*, indVars* );

//// model-specific stats variables
//typedef struct _pcGlobalStats {
//    double *NiR, *CiR, *MiR, **NijR, *z1iR;
//} pcGlobalStats;
//void *newModelStatsG_pc( xnDataBundle*, globalVars*, indVarBundle* );
//void freeModelStatsG_pc( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_pc( xnDataSet*, globalVars*, indVars* );
//void initializeVbHmmG_pc( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_pc( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputPcResults( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputPcResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_pc();
//void setGFunctions_pc();

double pTilde_z1_pc( int, void* );
double pTilde_zn_zn1_pc( int, int, void* );
double pTilde_xn_zn_pc( xnDataSet*, size_t, int, void* );

void calcStatsVars_pc( xnDataSet*, globalVars*, indVars* );
//void calcStatsVarsG_pc( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_pc( xnDataSet*, globalVars*, indVars* );
//void maximizationG_pc( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_pc( xnDataSet*, globalVars*, indVars* );
//double varLowerBoundG_pc( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_pc( xnDataSet*, globalVars*, indVars* );
//void reorderParametersG_pc( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_pc( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputResultsG_pc( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//
