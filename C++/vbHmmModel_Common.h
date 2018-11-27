/*
 *  vbHmmModel_Common.h
 *  Common VB-HMM model class.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#ifndef __vbHmmModel__Common__
#define __vbHmmModel__Common__

#include <stdio.h>
#include <string>
#include "vbHmm_Class.h"
#include "mathUtils.h"
#include "rand.h"

using namespace std;

class vbHmmModel{
public:
    int R;                      //  number of traces
    double dR;
    int sNo;                    // number of States

    // parameters
    Vec1D *uPiArr;
    double sumUPi;
    Mat2D *uAMat;
    Vec1D *sumUAArr;
    Vec1D *avgPi, *avgLnPi;
    Mat2D *avgA, *avgLnA;

    Vec1D *NiiR;
    Mat2D *NijR;
    Vec1D *z1iR;

    int iteration;              // calculation steps
    double maxLq;               // final lower bound
    Vec1D *LqArr;               // time series of lower bound

    vbHmmIndVars **ivs;

    vbHmmModel();
    vbHmmModel(const vbHmmModel&);
    ~vbHmmModel();

    virtual vbHmmModel *newInstance();
    
    virtual void setSNo(int _sNo);
    virtual void appendIndVars(int _N);

    virtual void initialize(vbHmmCond *c, vbHmmData *d);
    virtual void initialize_gmMat_random( vbHmmCond *c, vbHmmData *d );
    virtual void initialize_gmMat_kMeans( vbHmmCond *c, vbHmmData *d );
    virtual void initialize_modelParams(vbHmmCond *c, vbHmmData *d);
    virtual double pTilde_z1(int i);
    virtual double pTilde_zn_zn1(int i, int j);
    virtual double pTilde_xn_zn(vbHmmTrace *t, int n, int i);
    virtual void calcStatsVars(vbHmmData *d);
    virtual int maximization();
    virtual int maximization_modelVars(int s);
    virtual double varLowerBound(vbHmmData *d);
    virtual double varLowerBound_modelTerm(vbHmmData *d);
    virtual void reorderParameters(vbHmmData *d);
    virtual void outputResults(vbHmmData *d, fstream *logFS );
};


#endif /* defined(__vbHmmModel__Common__) */

//
