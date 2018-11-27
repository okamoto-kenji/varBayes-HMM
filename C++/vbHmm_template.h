/*
 *  vbHmm_template.h
 *  Model-specific class for VB-HMM-TEMPLATE.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#ifndef __vbHmm__Template__
#define __vbHmm__Template__


#include <stdio.h>
#include <string>
#include "vbHmm_Class.h"
#include "vbHmmModel_Common.h"

using namespace std;


class vbHmmCond_template : public vbHmmCond {
public:
    
    vbHmmCond_template();
    vbHmmCond_template(const vbHmmCond_template&);
    ~vbHmmCond_template();
};

class vbHmmTrace_template : public vbHmmTrace {
public:
    // define model-specific data variables

    vbHmmTrace_template();
    vbHmmTrace_template(const vbHmmTrace_template&);
    ~vbHmmTrace_template();

    virtual double kMeans_x(int n);         // return 1D variable for k-Means clustering at initialization

    // may define function to import data here
    void importDataFromFile(const char*, fstream*);
};

class vbHmmModel_template : public vbHmmModel {
public:
    // define any model-specific variables

    vbHmmModel_template();
    vbHmmModel_template(const vbHmmModel_template&);
    ~vbHmmModel_template();

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

#endif /* defined(__vbHmm__Template__) */

//
