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

// model-specific data
typedef struct _tsData {
    double dt;
    double time;
} tsData;

// model-specific parameters
typedef struct tsParameters_ {
    int sNo;
    double *uPiArr, sumUPi;
    double **uKMat, *sumUKArr;
    double *aIArr, *bIArr;
    double *avgPi, *avgLnPi, **avgK, **avgLnK, *avgLnKI, *avgI, *avgLnI;
    double *Ni, *Ti, *Nii, *Nij, **Mij;
} tsParameters;
tsParameters *blankParameters( int );

// model-specific output
void outputTsResults( vbHmmCommonParameters*, tsParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_ts();
void **mallocParameterArray_ts( size_t );
void *initialize_vbHmm_ts( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_ts( void* );
double pTilde_z1_ts( int, void* );
double pTilde_zn_zn1_ts( int, int, void* );
double pTilde_xn_zn_ts( xnDataSet*, size_t, int, void* );
void calcStatsVars_ts( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_ts( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_ts( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_ts( vbHmmCommonParameters*, void* );
void outputResults_ts( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//
