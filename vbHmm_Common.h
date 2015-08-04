/*
 *  vbHmm_Common.h
 *  Common VB-HMM engine.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#ifndef VBHMMCOMMON_DEF
#define VBHMMCOMMON_DEF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Structure encapsulating data to analyze
typedef struct _xnDataSet {
    size_t N;                       // number of data points
    double T;                       // total length of data (may not be always necessary)
    void *data;                     // an array conatining data points (arbitrary format)
} xnDataSet;

// Commonly used parameters
typedef struct _vbHmmCommonParameters {
    size_t dLen;                    // number of data points
    int sNo;                        // number of states
    
    double **gmMat;                 // gamma matrix for E-step
    double ***xiMat;                // xi matrix for E-step
    double **aMat;                  // alpha matrix for Baum-Welch
    double **bMat;                  // beta matrix for Baum-Welch
    double *cn;                     // scaling factor for Baum-Welch
    
    double **valpZnZn1, **valpXnZn; // temporary storage of calculation to save time
} vbHmmCommonParameters;

// Results of analysis
typedef struct _vbHmmResults {
    size_t iteration;               // calculation steps
    double maxLq;                   // final lower bound
    double *LqArr;                  // time series of lower bound
    int *maxSumTraj;                // resulting state transition trajectory
} vbHmmResults;


// manages the VB-HMM engine to choose likeliest model
int modelComparison( xnDataSet*, int, int, int, int, double, char*, FILE* );


//// Core Engines for VB-HMM
// common functions to manage peripheral data
vbHmmCommonParameters *initialize_Common( xnDataSet*, int );
void freeCommonParameters( vbHmmCommonParameters* );
vbHmmResults *newResults();
void freeResults( vbHmmResults* );
// main routine of VB-HMM
double vbHmm_Main( xnDataSet*, vbHmmCommonParameters*, void*, vbHmmResults* ,int ,double, FILE* );
// Baum-Welch algorithm
void forwardBackward( xnDataSet*, vbHmmCommonParameters*, void* );
// Viterbi algorithm
int *maxSum( xnDataSet*, vbHmmCommonParameters*, void*, int** );


//// functions to be implemented in model-specific sources
// to manage data
typedef void **(*mallocParameterArray_func)( size_t );
typedef void *(*initialize_vbHmm_func)( xnDataSet*, vbHmmCommonParameters* );
typedef void (*freeParameters_func)( void* );
// for E-step calculation
typedef double (*pTilde_z1_func)( int, void* );
typedef double (*pTilde_zn_zn1_func)( int, int, void* );
typedef double (*pTilde_xn_zn_func)( xnDataSet*, size_t, int, void* );
// for M-step calculation
typedef void (*calcStatsVars_func)( xnDataSet*, vbHmmCommonParameters*, void* );
typedef void (*maximization_func)( xnDataSet*, vbHmmCommonParameters*, void* );
// calculates the Variational Lower Bound
typedef double (*varLowerBound_func)( xnDataSet*, vbHmmCommonParameters*, void* );
// to manage & output results
typedef void (*reorderParameters_func)( vbHmmCommonParameters*, void* );
typedef void (*outputResults_func)( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );
// for stats simulation
typedef void *(*duplicateParameters_func)( void * );
typedef double (*calcAccuracyFromResult_func)( xnDataSet*, vbHmmResults* );
// function pointers
typedef struct _commonFunctions{
    mallocParameterArray_func  mallocParameterArray;
    initialize_vbHmm_func      initialize_vbHmm;
    freeParameters_func        freeParameters;
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
