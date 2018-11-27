/*
 *  vbHmm_GaussDiffusion.h
 *  Model-specific class for VB-HMM-GAUSS-DIFFUSION.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

//#include "vbHmmGaussDiffusion.h"
#ifndef __vbHmm__GaussDiffusion__
#define __vbHmm__GaussDiffusion__

#include <stdio.h>
#include <string>
#include "vbHmm_Class.h"
#include "vbHmmModel_Common.h"

using namespace std;


class vbHmmCond_GaussDiff : public vbHmmCond {
public:
    
    vbHmmCond_GaussDiff();
    vbHmmCond_GaussDiff(const vbHmmCond_GaussDiff&);
    ~vbHmmCond_GaussDiff();
};

class vbHmmTrace_GaussDiff : public vbHmmTrace {
public:
    Vec1D *xn;

    vbHmmTrace_GaussDiff();
    vbHmmTrace_GaussDiff(const vbHmmTrace_GaussDiff&);
    ~vbHmmTrace_GaussDiff();

    virtual double kMeans_x(int n);

    void readGaussDiffText(const char *filename, fstream *logFS );
private:
    Vec1D *readDoubleArrayFromFile(const char *filename, int *dlen, fstream *logFS);
};

class vbHmmModel_GaussDiff : public vbHmmModel {
public:
    Vec1D *avgDlt, *avgLnDlt;
    Vec1D *uAArr, *uBArr;
    Vec1D *aDlt, *bDlt;

    Vec1D *NiR;
    Vec1D *RiR;

    vbHmmModel_GaussDiff();
    vbHmmModel_GaussDiff(const vbHmmModel_GaussDiff&);
    ~vbHmmModel_GaussDiff();

    virtual vbHmmModel *newInstance();

    virtual void setSNo(int _sNo);
//    virtual void appendIndVars(int _N);
    
//    virtual void initialize(vbHmmCond *c, vbHmmData *d);
    virtual void initialize_modelParams(vbHmmCond *c, vbHmmData *d);
//    virtual double pTilde_z1(int i);
//    virtual double pTilde_zn_zn1(int i, int j);
    virtual double pTilde_xn_zn(vbHmmTrace *t, int n, int i);
    virtual void calcStatsVars(vbHmmData *d);
//    virtual int maximization();
    virtual int maximization_modelVars(int s);
//    virtual double varLowerBound(vbHmmData *d);
    virtual double varLowerBound_modelTerm(vbHmmData *d);
    virtual void reorderParameters(vbHmmData *d);
    virtual void outputResults(vbHmmData *d, fstream *logFS );
};

#endif /* defined(__vbHmm__GaussDiffusion__) */

//
