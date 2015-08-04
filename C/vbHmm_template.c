/*
 *  vbHmm_template.c
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

#include <math.h>
#include <string.h>

// Special functions is probably necessary.
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "vbHmm_template.h"

// For random number generation.
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))


void setFunctions__(){
    // Execute this to register functions before call common VB-HMM engine.
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray__;
    funcs.initialize_vbHmm = initialize_vbHmm__;
    funcs.freeParameters = freeParameters__;
    funcs.pTilde_z1 = pTilde_z1__;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1__;
    funcs.pTilde_xn_zn = pTilde_xn_zn__;
    funcs.calcStatsVars = calcStatsVars__;
    funcs.maximization = maximization__;
    funcs.varLowerBound = varLowerBound__;
    funcs.reorderParameters = reorderParameters__;
    funcs.outputResults = outputResults__;
    setFunctions( funcs );
}

void **mallocParameterArray__( n )
size_t n;
{
    return (void**)malloc( n * sizeof(parametersStruct*) );
}


// Redirect output function to that for model-specific data.
void outputResults__( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    output__Results( cParams, (parametersStruct*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm__( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    int sNo = cParams->sNo;
    parametersStruct *params = blankParameters( sNo );
    params->sNo = sNo;

    // Initialize parameters or any necessary for analysis.

    return params;
}


parametersStruct *blankParameters( sNo )
int sNo;
{
    // Provide a blank set of a model-specific parameters.
    parametersStruct *params = (parametersStruct*)malloc( sizeof(parametersStruct) );
    return params;
}

void freeParameters__( params )
void *params;
{
    // Responsible to free memory for a set of model-specific parameters.
    parametersStruct *p = (parametersStruct*)params;
    free( p );
    p = NULL;
}


double pTilde_z1__( i, params )
int i;
void *params;
{
    // Return the value for p(z_{1i}|pi_i) for given i.
    parametersStruct *p = (parametersStruct*)params;
    return 0.0;
}

double pTilde_zn_zn1__( i, j, params )
int i, j;
void *params;
{
    // Return the value for p(z_{nj}|z_{n-1,i},theta) for given i and j.
    parametersStruct *p = (parametersStruct*)params;
    return 0.0;
}

double pTilde_xn_zn__( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    // Return the value for p(x_{mi}|z_{ni},phi) for given i.
    parametersStruct *p = (parametersStruct*)params;
    return 0.0;
}


void calcStatsVars__( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    // Execute calculation statistics, which may be necessary for M-step.
    parametersStruct *p = (parametersStruct*)params;
}


void maximization__( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    // Calculation of M-step, depending on the model.
    parametersStruct *p = (parametersStruct*)params;
}


double varLowerBound__( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    // Calculation of M-step, depending on the model.
    parametersStruct *p = (parametersStruct*)params;
    return 0.0;
}


void reorderParameters__( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    // Reorder states in the final result, for example, in the order of intensity.
    parametersStruct *p = (parametersStruct*)params;
}


void output__Results( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
parametersStruct *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    // Output the final result (after reordering) in any format, e.g. stdout, stderr or files.
}

//
