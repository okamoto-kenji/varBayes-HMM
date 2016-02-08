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

#include "vbHmm_Common.h"

// model-specific data
typedef struct _gaussIntData {
    double *v;       // value
} gaussIntData;

// model-specific parameters
typedef struct gaussIntParameters_ {
    int sNo;
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double avgMu, avgLm, avgLnLm;
    double uBt, uMu, uA, uB;
    double btMu, aLm, bLm, mu0;
    double **Nij, *Nii, *Ni, N0, N0i, Nx;
} gaussIntParameters;
gaussIntParameters *blankParameters_gaussInt( int );
void freegaussIntDataSet( xnDataSet* );

// model-specific output
void outputGaussIntResults( vbHmmCommonParameters*, gaussIntParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gaussInt();
void **mallocParameterArray_gaussInt( size_t );
void *initialize_vbHmm_gaussInt( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_gaussInt( void* );

double pTilde_z1_gaussInt( int, void* );
double pTilde_zn_zn1_gaussInt( int, int, void* );
double pTilde_xn_zn_gaussInt( xnDataSet*, size_t, int, void* );

void calcStatsVars_gaussInt( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_gaussInt( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_gaussInt( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_gaussInt( vbHmmCommonParameters*, void* );
void outputResults_gaussInt( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//

