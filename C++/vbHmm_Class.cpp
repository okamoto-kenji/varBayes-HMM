/*
 *  vbHmm_Class.cpp
 *  Common VB-HMM classes.
 *  Reference: 
 *    Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#include "vbHmm_Class.h"


//////////  conditions
vbHmmCond::vbHmmCond(){
}

vbHmmCond::vbHmmCond(const vbHmmCond &other){  // do not copy
    sFrom = other.sFrom;
    sTo = other.sTo;
    trials = other.trials;
    maxIteration = other.maxIteration;
    lbTh = other.lbTh;
}

vbHmmCond::~vbHmmCond(){
}


//////////  trace
vbHmmTrace::vbHmmTrace(){
    name = "";
    N = 0;                       // number of data points
}

vbHmmTrace::vbHmmTrace(const vbHmmTrace &other){  // do not copy
    name = "";
    N = 0;
}

vbHmmTrace::~vbHmmTrace(){
}

double vbHmmTrace::kMeans_x( int n ){
    return 0.0;
}


//////////  data
vbHmmData::vbHmmData(){
    name = "";
    R = 0;                       // number of traces
    dR = (double)R;
    totalN = 0;
    traces = NULL;
}

vbHmmData::vbHmmData(const vbHmmData &other){  // do not copy
    name = "";
    R = 0;
    dR = (double)R;
    totalN = 0;
    traces = NULL;
}

vbHmmData::~vbHmmData(){
    for( int i = 0 ; i < R ; i++ ){
        delete traces[i];
        traces[i] = NULL;
    }
    free(traces);
}

void vbHmmData::appendTrace( vbHmmTrace *t ){
    if( t != NULL ){
        if( R == 0 ){
            traces = (vbHmmTrace**)malloc( sizeof(vbHmmTrace*) );
            name = t->name;
        } else {
            traces =(vbHmmTrace**)realloc( traces, (R + 1) * sizeof(vbHmmTrace*) );
        }
        traces[R] = t;
        R++;
        dR = (double)R;
        totalN += t->N;
    }
}


//////////  individual variables
vbHmmIndVars::vbHmmIndVars( int _N, int _sNo ){
    N = _N;                       // number of data points

    gmMat = new Mat2D(_N, _sNo);
    xiMat = new Ten3D(_N, _sNo, _sNo);
    aMat = new Mat2D(_N, _sNo);
    bMat = new Mat2D(_N, _sNo);
    cn = new Vec1D(_N);
    valpZnZn1 = new Mat2D(_N, _sNo);
    valpXnZn = new Mat2D(_N, _sNo);
    
#ifdef OUTPUT_MAX_GAMMA
    gammaTraj = new Vec1I(_N);
#endif
    stateTraj = new Vec1I(_N);
}

vbHmmIndVars::vbHmmIndVars(const vbHmmIndVars &other){  // do not copy
    gmMat = NULL;
    xiMat = NULL;
    aMat = NULL;
    bMat = NULL;
    cn = NULL;
    valpZnZn1 = NULL;
    valpXnZn = NULL;
    
#ifdef OUTPUT_MAX_GAMMA
    gammaTraj = NULL;
#endif
    stateTraj = NULL;
}

vbHmmIndVars::~vbHmmIndVars(){
    delete gmMat;
    delete xiMat;
    delete aMat;
    delete bMat;
    delete cn;
    delete valpZnZn1;
    delete valpXnZn;
    
#ifdef OUTPUT_MAX_GAMMA
    delete gammaTraj;
#endif
    delete stateTraj;
}


//
