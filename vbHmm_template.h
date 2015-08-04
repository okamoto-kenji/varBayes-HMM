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

// Include common VB-HMM engine.
#include "vbHmm_Common.h"

#define  SIMU_STATS

int modelComparison__( xnDataSet*, int, int, int, int, double, char*, FILE* );

// Specify variables to describe each data point.
typedef struct _dataStruct {
} dataStruct;

typedef struct _parametersStruct {
    // Assumed number of state is kept for convenience.
    int sNo;
    // Define model specific parameters here.
    // For example, average and log average of parameters,
    //   hyperparameter for prior distributions.
} parametersStruct;

parametersStruct *blankParameters( int );
void output__Results( vbHmmCommonParameters*, parametersStruct*, vbHmmResults*, int, char*, FILE* );


// functions required to work with vbHmm_Common
void setFunctions__();
void **mallocParameterArray__( size_t );
void *initialize_vbHmm__( xnDataSet*, vbHmmCommonParameters* );
void freeParameters__( void* );

double p_z1__( int, void* );
double pTilde_z1__( int, void* );
double p_zn_zn1__( int, int, void* );
double pTilde_zn_zn1__( int, int, void* );
double p_xn_zn__( xnDataSet*, size_t, int, void* );
double pTilde_xn_zn__( xnDataSet*, size_t, int, void* );

void calcStatsVars__( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization__( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound__( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters__( vbHmmCommonParameters*, void* );
void outputResults__( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//
