/*
 *  vbHmmTsFret.h
 *  Model-specific core functions for VB-HMM-TS-FRET.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#ifndef VBHMMTSFRET_DEF
#define VBHMMTSFRET_DEF

#include "vbHmm_Common.h"
//#include "gVbHmm_Common.h"

int modelComparison_tsFret( xnDataSet*, int, int, int, int, double, char*, FILE* );

// model-specific data
typedef struct _tsFretData {
    double T;                       // total length of data (may not be always necessary)
    double *dt;
    double *time;
    int *ch;             //  0:donor,  1:acceptor
} tsFretData;

xnDataSet *newXnDataSet_tsFret( const char* );
void freeXnDataSet_tsFret( xnDataSet** );

// model-specific parameters
typedef struct tsFretParameters_ {
    double *uPiArr, sumUPi;
    double *aIArr, *bIArr;
    double **uKMat, *sumUKArr;
    double *uEArr, *vEArr;
    double *avgPi, *avgLnPi, **avgK, **avgLnK, *avgLnKI, *avgI, *avgLnI;
    double *avgE, **avgLnE;
} tsFretParameters;

void *newModelParameters_tsFret( xnDataSet*, int );
void freeModelParameters_tsFret( void**, xnDataSet*, int );

// model-specific stats variables
typedef struct _tsFretStats {
    double *Ni, *Ti, *eps, *Nii, *Nij, **Mij;
} tsFretStats;

void *newModelStats_tsFret( xnDataSet*, globalVars*, indVars* );
void freeModelStats_tsFret( void**, xnDataSet*, globalVars*, indVars* );

// model-specific stats variables
//typedef struct _tsFretGlobalStats {
//    double *NiR, *TiR, *epsR, *NiiR, *NijR, **MijR, *z1iR;
//} tsFretGlobalStats;

//void *newModelStatsG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
//void freeModelStatsG_tsFret( void**, xnDataBundle*, globalVars*, indVarBundle* );

void initializeVbHmm_tsFret( xnDataSet*, globalVars*, indVars* );
//void initializeVbHmmG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
void initialize_indVars_tsFret( xnDataSet*, globalVars*, indVars* );

// model-specific output
void outputTsFretResults( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputTsFretResultsG( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_tsFret();
//void setGFunctions_tsFret();

double pTilde_z1_tsFret( int, void* );
double pTilde_zn_zn1_tsFret( int, int, void* );
double pTilde_xn_zn_tsFret( xnDataSet*, size_t, int, void* );

void calcStatsVars_tsFret( xnDataSet*, globalVars*, indVars* );
//void calcStatsVarsG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
void maximization_tsFret( xnDataSet*, globalVars*, indVars* );
//void maximizationG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
double varLowerBound_tsFret( xnDataSet*, globalVars*, indVars* );
//double varLowerBoundG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
void reorderParameters_tsFret( xnDataSet*, globalVars*, indVars* );
//void reorderParametersG_tsFret( xnDataBundle*, globalVars*, indVarBundle* );
void outputResults_tsFret( xnDataSet*, globalVars*, indVars*, FILE* );
//void outputResultsG_tsFret( xnDataBundle*, globalVars*, indVarBundle*, FILE* );

#endif

//

