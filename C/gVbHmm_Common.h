/*
 *  gVbHmm_Common.h
 *  Common VB-HMM engine for global analysis.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#ifndef GVBHMMCOMMON_DEF
#define GVBHMMCOMMON_DEF

#include "vbHmm_Common.h"

// Structure encapsulating data to analyze
typedef struct _xnDataBundle {
    int R;                          // number of trajectories
    xnDataSet **xn;                       // an array conatining data points (arbitrary format)
} xnDataBundle;

// Structure encapsulating data to analyze
typedef struct _indVarBundle {
    void *stats;
    indVars **indVars;                     // an array conatining data points (arbitrary format)
} indVarBundle;

// manages the VB-HMM engine to choose likeliest model
int gModelComparison( xnDataBundle*, int, int, int, int, double, FILE* );

globalVars *newGlobalVarsG( xnDataBundle*, int );
void freeGlobalVarsG( xnDataBundle*, globalVars** );

typedef void *(*new_model_statsG_func)( xnDataBundle*, globalVars*, indVarBundle* );
typedef void (*free_model_statsG_func)( void**, xnDataBundle*, globalVars*, indVarBundle* );

//// Core Engines for VB-HMM
// main routine of VB-HMM
double gVbHmm_Main( xnDataBundle*, globalVars*, indVarBundle*, int, double, FILE* );


//// functions to be implemented in model-specific sources
// to initialize parameters
typedef void (*initialize_vbHmmG_func)( xnDataBundle*, globalVars*, indVarBundle* );
// for M-step calculation
typedef void (*calcStatsVarsG_func)( xnDataBundle*, globalVars*, indVarBundle* );
typedef void (*maximizationG_func)( xnDataBundle*, globalVars*, indVarBundle* );
// calculates the Variational Lower Bound
typedef double (*varLowerBoundG_func)( xnDataBundle*, globalVars*, indVarBundle* );
// to output results
typedef void (*reorderParametersG_func)( xnDataBundle*, globalVars*, indVarBundle* );
typedef void (*outputResultsG_func)( xnDataBundle*, globalVars*, indVarBundle*, FILE* );
// function pointers
typedef struct _gCommonFunctions{
    new_model_parameters_func   newModelParameters;
    free_model_parameters_func  freeModelParameters;
    new_model_stats_func        newModelStats;
    free_model_stats_func       freeModelStats;
    new_model_statsG_func       newModelStatsG;
    free_model_statsG_func      freeModelStatsG;
    initialize_vbHmmG_func      initializeVbHmmG;
    pTilde_z1_func              pTilde_z1;
    pTilde_zn_zn1_func          pTilde_zn_zn1;
    pTilde_xn_zn_func           pTilde_xn_zn;
    calcStatsVarsG_func         calcStatsVarsG;
    maximizationG_func          maximizationG;
    varLowerBoundG_func         varLowerBoundG;
    reorderParametersG_func     reorderParametersG;
    outputResultsG_func         outputResultsG;
} gCommonFunctions;
void setGFunctions( gCommonFunctions );


#endif

//
