/*
 *  vbHmm_Common.h
 *  Common VB-HMM engine.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#ifndef VBHMMCOMMON_DEF
#define VBHMMCOMMON_DEF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _individualDataVariables{    // variables for individual trajectories
    // stats
    double **gmMat;                     // gamma matrix for E-step
    double ***xiMat;                    // xi matrix for E-step
    double **aMat;                      // alpha matrix for Baum-Welch
    double **bMat;                      // beta matrix for Baum-Welch
    double *cn;                         // scaling factor for Baum-Welch
    double **valpZnZn1, **valpXnZn;     // temporary storage of calculation to save time
    void *stats;

    // results
    int *stateTraj;
} indVars;

typedef struct _globalVariables{        // global variables
    // parameters
    int sNo;
    void *params;

    // results
    size_t iteration;               // calculation steps
    double maxLq;                   // final lower bound
    double *LqArr;                  // time series of lower bound
} globalVars;


// Structure encapsulating data to analyze
typedef struct _xnDataSet {
    char *name;
    size_t N;                       // number of data points
    void *data;                     // an array conatining data points (arbitrary format)
} xnDataSet;


// manages the VB-HMM engine to choose likeliest model
int modelComparison( xnDataSet*, int, int, int, int, double, FILE* );


//// Core Engines for VB-HMM
// common functions to manage peripheral data
globalVars *newGlobalVars( xnDataSet*, int );
void freeGlobalVars( xnDataSet*, globalVars** );
indVars *newIndVars( xnDataSet*, globalVars* );
void freeIndVars( xnDataSet*, globalVars*, indVars** );
// main routine of VB-HMM
double vbHmm_Main( xnDataSet*, globalVars*, indVars* ,int ,double, FILE* );
// Baum-Welch algorithm
void forwardBackward( xnDataSet*, globalVars*, indVars* );
// Viterbi algorithm
int *maxSum( xnDataSet*, globalVars*, indVars* );


//// functions to be implemented in model-specific sources
// to manage data
typedef void *(*new_model_parameters_func)( xnDataSet*, int );
typedef void (*free_model_parameters_func)( void**, xnDataSet*, int );
typedef void *(*new_model_stats_func)( xnDataSet*, globalVars*, indVars* );
typedef void (*free_model_stats_func)( void**, xnDataSet*, globalVars*, indVars* );
typedef void (*initialize_vbHmm_func)( xnDataSet*, globalVars*, indVars* );
// for E-step calculation
typedef double (*pTilde_z1_func)( int, void* );
typedef double (*pTilde_zn_zn1_func)( int, int, void* );
typedef double (*pTilde_xn_zn_func)( xnDataSet*, size_t, int, void* );
// for M-step calculation
typedef void (*calcStatsVars_func)( xnDataSet*, globalVars*, indVars* );
typedef void (*maximization_func)( xnDataSet*, globalVars*, indVars* );
// calculates the Variational Lower Bound
typedef double (*varLowerBound_func)( xnDataSet*, globalVars*, indVars* );
// to manage & output results
typedef void (*reorderParameters_func)( xnDataSet*, globalVars*, indVars* );
typedef void (*outputResults_func)( xnDataSet*, globalVars*, indVars*, FILE* );
// function pointers
typedef struct _commonFunctions{
    new_model_parameters_func  newModelParameters;
    free_model_parameters_func freeModelParameters;
    new_model_stats_func       newModelStats;
    free_model_stats_func      freeModelStats;
    initialize_vbHmm_func      initializeVbHmm;
    pTilde_z1_func             pTilde_z1;
    pTilde_zn_zn1_func         pTilde_zn_zn1;
    pTilde_xn_zn_func          pTilde_xn_zn;
    calcStatsVars_func         calcStatsVars;
    maximization_func          maximization;
    varLowerBound_func         varLowerBound;
    reorderParameters_func     reorderParameters;
    outputResults_func         outputResults;
} commonFunctions;
void setFunctions( commonFunctions );


#endif

//
