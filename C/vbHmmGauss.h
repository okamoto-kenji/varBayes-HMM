/*
 *  vbHmmGauss.h
 *  Model-specific core functions for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.xx
 */

#ifndef VBHMMGAUSS_DEF
#define VBHMMGAUSS_DEF

#include "vbHmm_Common.h"

// model-specific data
typedef struct _gaussData {
//    double binSize;
    double *v;       // value
//    unsigned int *dCounts;
//    unsigned int *aCounts;
} gaussData;

// model-specific parameters
typedef struct gaussParameters_ {
    int sNo;
//    double binSize;
    double *uPiArr, sumUPi;
//    double *aIArr, *bIArr;
    double **uAMat, *sumUAArr;
//    double *uEArr, *vEArr;
    double *avgPi, *avgLnPi;
    double **avgA, **avgLnA;
    double *avgMu, *avgLm, *avgLnLm;
    double *uBtArr, *uMuArr, *uAArr, *uBArr;
    double *btMu, *aLm, *bLm, *mu0;
//    double *avgI, *avgLnI;
//    double *avgE, **avgLnE;
    double **Nij, *Nii, *Ni, *barX, *NiSi;
//    double *Ci, *Di, *Ai, *Mi, **Nij;

//    int sNo;
//    double binSize;
//    double *uPiArr, sumUPi;
//    double *aIArr, *bIArr;
//    double **uAMat, *sumUAArr;
//    double *uEArr, *vEArr;
//    double *avgPi, *avgLnPi, **avgA, **avgLnA, *avgI, *avgLnI;
//    double *avgE, **avgLnE;
//    double *Ni, *Ci, *Di, *Ai, *Mi, **Nij;
} gaussParameters;
gaussParameters *blankParameters_gauss( int );
void freeGaussDataSet( xnDataSet* );

// model-specific output
void outputGaussResults( vbHmmCommonParameters*, gaussParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_gauss();
void **mallocParameterArray_gauss( size_t );
void *initialize_vbHmm_gauss( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_gauss( void* );

double pTilde_z1_gauss( int, void* );
double pTilde_zn_zn1_gauss( int, int, void* );
double pTilde_xn_zn_gauss( xnDataSet*, size_t, int, void* );

void calcStatsVars_gauss( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_gauss( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_gauss( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_gauss( vbHmmCommonParameters*, void* );
void outputResults_gauss( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//

