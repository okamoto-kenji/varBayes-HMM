/*
 *  vbHmm_Class.h
 *  Common VB-HMM classes.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#ifndef __vbHmm__Class__
#define __vbHmm__Class__

#include <stdio.h>
#include <string>
#include "mathUtils.h"

// uncomment to enable calculation & output of max gamma trajectories
//#define OUTPUT_MAX_GAMMA

using namespace std;

enum {
    initType_random = 0,
    initType_kMeans = 1
};

class vbHmmCond{
public:
    int sFrom, sTo;
    int trials;
    int maxIteration;
    double lbTh;                        // threshold for lower bound
    int initType;

    vbHmmCond();
    vbHmmCond(const vbHmmCond&);
    ~vbHmmCond();
};

class vbHmmTrace{
public:
    string name;
    int N;                              // number of data points
    
    vbHmmTrace();
    vbHmmTrace(const vbHmmTrace&);
    ~vbHmmTrace();

    virtual double kMeans_x(int n);
};

class vbHmmData{
public:
    string name;
    int R;                              // number of traces
    double dR;
    int totalN;

    vbHmmTrace **traces;

    vbHmmData();
    vbHmmData(const vbHmmData&);
    ~vbHmmData();

    virtual void appendTrace(vbHmmTrace*t);
};

class vbHmmIndVars{                     // individual variables
public:
    int N;

    Mat2D *gmMat;                       // gamma matrix for E-step
    Ten3D *xiMat;                       // xi matrix for E-step
    Mat2D *aMat;                        // alpha matrix for Baum-Welch
    Mat2D *bMat;                        // beta matrix for Baum-Welch
    Vec1D *cn;                          // scaling factor for Baum-Welch
    Mat2D *valpZnZn1, *valpXnZn;        // temporary storage of calculation to save time

    vbHmmTrace **traces;

    // results
#ifdef OUTPUT_MAX_GAMMA
    Vec1I *gammaTraj;                   // state trajectory by max gamma
#endif
    Vec1I *stateTraj;                   // state trajectory by max sum

    vbHmmIndVars(int _N, int _sNo);
    vbHmmIndVars(const vbHmmIndVars&);
    ~vbHmmIndVars();
};


#endif /* defined(__vbHmm__Class__) */

//
