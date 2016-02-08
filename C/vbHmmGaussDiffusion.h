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

#include "vbHmm_Common.h"

// model-specific data
typedef struct _gaussDiffData {
    double *v;       // value
} gaussDiffData;

// model-specific parameters
typedef struct gaussDiffParameters_ {
    int sNo;
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double *avgDlt, *avgLnDlt;
    double *uAArr, *uBArr;
    double *aDlt, *bDlt;
    double *Ni, *Ri, *Nii, **Nij;
} gaussDiffParameters;
gaussDiffParameters *blankParameters_gaussDiff( int );
void freeGaussDiffDataSet( xnDataSet* );

// model-specific output
void outputGaussDiffResults( vbHmmCommonParameters*, gaussDiffParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gaussDiff();
void **mallocParameterArray_gaussDiff( size_t );
void *initialize_vbHmm_gaussDiff( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_gaussDiff( void* );

double pTilde_z1_gaussDiff( int, void* );
double pTilde_zn_zn1_gaussDiff( int, int, void* );
double pTilde_xn_zn_gaussDiff( xnDataSet*, size_t, int, void* );

void calcStatsVars_gaussDiff( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_gaussDiff( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_gaussDiff( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_gaussDiff( vbHmmCommonParameters*, void* );
void outputResults_gaussDiff( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//

