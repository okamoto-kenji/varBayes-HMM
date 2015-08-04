/*
 *  vbHmmTsFret.h
 *  Model-specific core functions for VB-HMM-TS-FRET.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#ifndef VBHMMTSFRET_DEF
#define VBHMMTSFRET_DEF

#include "vbHmm_Common.h"

int modelComparison_tsFret( xnDataSet*, int, int, int, int, double, char*, FILE* );

// model-specific data
typedef struct _tsFretData {
    double dt;
    double time;
    int ch;             //  0:donor,  1:acceptor
} tsFretData;

// model-specific parameters
typedef struct tsFretParameters_ {
    int sNo;
    double *uPiArr, sumUPi;
    double *aIArr, *bIArr;
    double **uKMat, *sumUKArr;
    double *uEArr, *vEArr;
    double *avgPi, *avgLnPi, **avgK, **avgLnK, *avgLnKI, *avgI, *avgLnI;
    double *avgE, **avgLnE;
    double *Ni, *Ti, *eps, *Nii, *Nij, **Mij;
} tsFretParameters;
tsFretParameters *blankParameters_tsFret( int );

// model-specific output
void outputTsFretResults( vbHmmCommonParameters*, tsFretParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_tsFret();
void **mallocParameterArray_tsFret( size_t );
void *initialize_vbHmm_tsFret( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_tsFret( void* );

double pTilde_z1_tsFret( int, void* );
double pTilde_zn_zn1_tsFret( int, int, void* );
double pTilde_xn_zn_tsFret( xnDataSet*, size_t, int, void* );

void calcStatsVars_tsFret( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_tsFret( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_tsFret( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_tsFret( vbHmmCommonParameters*, void* );
void outputResults_tsFret( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//

