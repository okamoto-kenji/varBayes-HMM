/*
 *  mathUtils.h
 *  for vbHMM
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Last modified on 2018.11.26
 */

#ifndef __vbHMM__mathUtils__
#define __vbHMM__mathUtils__

#include "rand.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>

#define  VBHMM_GSL_ERROR_ON_LnGAMMA  999
#define  VBHMM_GSL_ERROR_ON_FACT     888
#define  VBHMM_GSL_ERROR_ON_PSI      777
#define  VBHMM_GSL_NAN               666

void disableGslErrorHandler();

double lnGamma(double x);
double fact(double x);
double psi(double x);

class Vec1D {
public:
    int d1;
    double *p;

    Vec1D(int _d1=0);
    Vec1D(const Vec1D&);
    ~Vec1D();

    void resize(int);
    void allSet(double);
    void append(double);
    double var();
    double sdev();
};

class Vec1I {
public:
    int d1;
    int *p;
    
    Vec1I(int _d1=0);
    Vec1I(const Vec1I&);
    ~Vec1I();
    
    void resize(int);
    void allSet(int);
    double var();
    double sdev();
};

class Mat2D {
public:
    int d1, d2;
    double **p;

    Mat2D(int _d1=0, int _d2=0);
    Mat2D(const Mat2D&);
    ~Mat2D();
//
//    void allSet(double);
};

class Ten3D {
public:
    int d1, d2, d3;
    double ***p;
    
    Ten3D(int _d1=0, int _d2=0, int _d3=0);
    Ten3D(const Ten3D&);
    ~Ten3D();
};


#endif /* defined(__vbHmm__mathUtils__) */

//

