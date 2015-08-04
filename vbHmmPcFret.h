/*
 *  vbHmmPcFret.h
 *  Model-specific core functions for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#ifndef VBHMMPCFRET_DEF
#define VBHMMPCFRET_DEF

#include "vbHmm_Common.h"

// model-specific data
typedef struct _pcFretData {
    double binSize;
    unsigned int *dCounts;
    unsigned int *aCounts;
} pcFretData;

// model-specific parameters
typedef struct pcFretParameters_ {
    int sNo;
    double binSize;
    double *uPiArr, sumUPi;
    double *aIArr, *bIArr;
    double **uAMat, *sumUAArr;
    double *uEArr, *vEArr;
    double *avgPi, *avgLnPi, **avgA, **avgLnA, *avgI, *avgLnI;
    double *avgE, **avgLnE;
    double *Ni, *Ci, *Di, *Ai, *Mi, **Nij;
} pcFretParameters;
pcFretParameters *blankParameters_pcFret( int );
void freePcFretDataSet( xnDataSet* );

// model-specific output
void outputPcFretResults( vbHmmCommonParameters*, pcFretParameters*, vbHmmResults*, int, char*, FILE* );

// functions required to work with vbHmm_Common
void setFunctions_pcFret();
void **mallocParameterArray_pcFret( size_t );
void *initialize_vbHmm_pcFret( xnDataSet*, vbHmmCommonParameters* );
void freeParameters_pcFret( void* );

double pTilde_z1_pcFret( int, void* );
double pTilde_zn_zn1_pcFret( int, int, void* );
double pTilde_xn_zn_pcFret( xnDataSet*, size_t, int, void* );

void calcStatsVars_pcFret( xnDataSet*, vbHmmCommonParameters*, void* );
void maximization_pcFret( xnDataSet*, vbHmmCommonParameters*, void* );
double varLowerBound_pcFret( xnDataSet*, vbHmmCommonParameters*, void* );
void reorderParameters_pcFret( vbHmmCommonParameters*, void* );
void outputResults_pcFret( vbHmmCommonParameters*, void*, vbHmmResults*, int, char*, FILE* );

#endif

//

