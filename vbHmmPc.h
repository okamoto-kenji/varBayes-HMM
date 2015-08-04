/*
 *  vbHmmPc.h
 *  Model-specific core functions for VB-HMM-PC.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */


#ifndef VBHMMPC_DEF
#define VBHMMPC_DEF

#include "vbHmm_Common.h"

// model-specific data
typedef struct _pcData {
    double binSize;
    unsigned int *counts;
} pcData;

// model-specific parameters
typedef struct pcParameters_ {
    int sNo;
    double binSize;
    double *uPiArr, sumUPi;
    double **uAMat, *sumUAArr;
    double *aIArr, *bIArr;
    double *avgPi, *avgLnPi, **avgA, **avgLnA, *avgI, *avgLnI;
    double *Ni, *Ci, *Mi, **Nij;
} pcParameters;
pcParameters *blankParameters_pc( int );
void freePcDataSet( xnDataSet* );

// model-specific output
void outputPcResults( vbHmmCommonParameters*, pcParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_pc();
void **mallocParameterArray_pc( size_t );
void *initialize_vbHmm_pc( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_pc( void* );

double pTilde_z1_pc( int, void* );
double pTilde_zn_zn1_pc( int, int, void* );
double pTilde_xn_zn_pc( xnDataSet*, size_t, int, void* );

void calcStatsVars_pc( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_pc( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_pc( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_pc( vbHmmCommonParameters*, void* );
void outputResults_pc( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//
