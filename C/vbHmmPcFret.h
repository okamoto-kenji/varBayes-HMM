/*
 *  vbHmmPcFret.h
 *  Model-specific core functions for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#ifndef VBHMMPCFRET_DEF
#define VBHMMPCFRET_DEF

#include "vbHmm_Common.h"
//#include "gVbHmm_Common.h"

// model-specific data
typedef struct _pcFretData {
    double binSize;
    unsigned int *dCounts;
    unsigned int *aCounts;
} pcFretData;

xnDataSet *newXnDataSet_pcFret( const char* );
void freeXnDataSet_pcFret( xnDataSet** );

// model-specific parameters
typedef struct pcFretParameters_ {
    double *uPiArr, sumUPi;
    double *aIArr, *bIArr;
    double **uAMat, *sumUAArr;
    double *uEArr, *vEArr;
    double *avgPi, *avgLnPi, **avgA, **avgLnA, *avgI, *avgLnI;
    double *avgE, **avgLnE;
} pcFretParameters;

void *newModelParameters_pcFret( xnDataSet*, int );
void freeModelParameters_pcFret( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _pcFretStats {
    double *Ni, *Ci, *Di, *Ai, *Mi, **Nij;
} pcFretStats;

void *newModelStats_pcFret( xnDataSet*, globalVars*, indVars* );
void freeModelStats_pcFret( void**, xnDataSet*, globalVars*, indVars* );

//// model-specific stats variables
//typedef struct _pcFretGlobalStats {
//    double *NiR, *CiR, *DiR, *AiR, *MiR, **NijR, *z1iR;
//} pcFretGlobalStats;
//void *newModelStatsG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
//void freeModelStatsG_pcFret( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_pcFret( xnDataSet*, globalVars*, indVars* );
//void initializeVbHmmG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_pcFret( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputPcFretResults( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputPcFretResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_pcFret();
//void setGFunctions_pcFret();

double pTilde_z1_pcFret( int, void* );
double pTilde_zn_zn1_pcFret( int, int, void* );
double pTilde_xn_zn_pcFret( xnDataSet*, size_t, int, void* );

void calcStatsVars_pcFret( xnDataSet*, globalVars*, indVars* );
//void calcStatsVarsG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_pcFret( xnDataSet*, globalVars*, indVars* );
//void maximizationG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_pcFret( xnDataSet*, globalVars*, indVars* );
//double varLowerBoundG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_pcFret( xnDataSet*, globalVars*, indVars* );
//void reorderParametersG_pcFret( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_pcFret( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputResultsG_pcFret( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//

