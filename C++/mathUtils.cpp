/*
 *  mathUtils.cpp
 *  for vbHMM
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Last modified on 2018.11.26
 */

#include <iostream>
#include "mathUtils.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_errno.h>


void disableGslErrorHandler(){
    gsl_set_error_handler_off();
}

double lnGamma( double x ){
    gsl_sf_result r;
    int s = gsl_sf_lngamma_e( x, &r );
    if( s == GSL_SUCCESS){
        return r.val;
    } else {
        throw VBHMM_GSL_ERROR_ON_LnGAMMA;
        return 0.0;
    }
}

double fact( double x ){
    gsl_sf_result r;
    int s = gsl_sf_fact_e( x, &r );
    if( s == GSL_SUCCESS){
        return r.val;
    } else {
        throw VBHMM_GSL_ERROR_ON_FACT;
        return 0.0;
    }
}

double psi( double x ){
    gsl_sf_result r;
    int s = gsl_sf_psi_e( x, &r );
    if( s == GSL_SUCCESS){
        return r.val;
    } else {
        throw VBHMM_GSL_ERROR_ON_PSI;
        return 0.0;
    }
}


Vec1D::Vec1D(int _d1){
    d1 = _d1;
    if( d1 > 0 ){
        p = (double*)malloc(sizeof(double) * d1);
    } else {
        p = NULL;
    }
}

Vec1D::Vec1D(const Vec1D &other){
    d1 = other.d1;
    if( d1 > 0 ){
        p = (double*)malloc(sizeof(double) * d1);
        int i;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = other.p[i];
        }
    } else {
        p = NULL;
    }
}

Vec1D::~Vec1D(){
    free(p);
}

void Vec1D::resize(int d){
    if( d != d1 ){
        if( d == 0 ){
            if( p != NULL ){
                free(p);
                p = NULL;
            }
        } else {
            p = (double*)realloc( p, d * sizeof(double) );
            if( d > d1 ){
                for( int i = d1 ; i < d ; i++ ){
                    p[i] = 0.0;
                }
            }
        }
        d1 = d;
    }
}

void Vec1D::allSet(double v){
    int i;
    for( i = 0 ; i < d1 ; i++ ){
        p[i] = v;
    }
}

void Vec1D::append(double v){
    resize( d1 + 1 );
    if( d1 > 0 )  p[d1-1] = v;
}

double Vec1D::var(){
    double n = (double)d1, mean = 0.0, ms = 0.0;
    int i;
    for( i = 0 ; i < d1 ; i++ ){
        mean += p[i];
        ms += p[i] * p[i];
    }
    mean /= n;
    ms /= n;
    return (ms - mean * mean) * n / (n - 1.0);
}


double Vec1D::sdev(){
    return sqrt(var());
}


Vec1I::Vec1I(int _d1){
    d1 = _d1;
    if( d1 > 0 ){
        p = (int*)malloc(sizeof(int) * d1);
    } else {
        p = NULL;
    }
}

Vec1I::Vec1I(const Vec1I &other){
    d1 = other.d1;
    if( d1 > 0 ){
        p = (int*)malloc(sizeof(int) * d1);
        int i;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = other.p[i];
        }
    } else {
        p = NULL;
    }
}

Vec1I::~Vec1I(){
    free(p);
}

void Vec1I::resize(int d){
    if( d != d1 ){
        if( d == 0 ){
            if( p != NULL ){
                free(p);
                p = NULL;
            }
        } else {
            p = (int*)realloc( p, d * sizeof(int) );
            if( d > d1 ){
                for( int i = d1 ; i < d ; i++ ){
                    p[i] = 0;
                }
            }
        }
        d1 = d;
    }
}

void Vec1I::allSet(int v){
    int i;
    for( i = 0 ; i < d1 ; i++ ){
        p[i] = v;
    }
}


Mat2D::Mat2D(int _d1, int _d2){
    d1 = _d1;
    d2 = _d2;
    if( (d1 > 0) && (d2 > 0) ){
        p = (double **)malloc(sizeof(double*) * d1);
        int i;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = (double*)malloc(sizeof(double) * d2);
        }
    } else {
        p = NULL;
    }
}

Mat2D::Mat2D(const Mat2D &other){
    d1 = other.d1;
    d2 = other.d2;
    if( (d1 > 0) && (d2 > 0) ){
        p = (double **)malloc(sizeof(double*) * d1);
        int i, j;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = (double*)malloc(sizeof(double) * d2);
            for( j = 0 ; j < d2 ; j++ ){
                p[i][j] = other.p[i][j];
            }
        }
    } else {
        p = NULL;
    }
}

Mat2D::~Mat2D(){
    if( p != NULL ){
        int i;
        for( i = 0 ; i < d1 ; i++ ){
            free(p[i]);
        }
        free(p);
    }
}


Ten3D::Ten3D(int _d1, int _d2, int _d3){
    d1 = _d1;
    d2 = _d2;
    d3 = _d3;
    if( (d1 > 0) && (d2 > 0) && (d3 > 0) ){
        p = (double ***)malloc(sizeof(double**) * d1);
        int i, j;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = (double**)malloc(sizeof(double*) * d2);
            for( j = 0 ; j < d2 ; j++ ){
                p[i][j] = (double*)malloc(sizeof(double) * d3);
            }
        }
    } else {
        p = NULL;
    }
}

Ten3D::Ten3D(const Ten3D &other){
    d1 = other.d1;
    d2 = other.d2;
    d3 = other.d3;
    if( (d1 > 0) && (d2 > 0) && (d3 > 0) ){
        p = (double ***)malloc(sizeof(double**) * d1);
        int i, j, k;
        for( i = 0 ; i < d1 ; i++ ){
            p[i] = (double**)malloc(sizeof(double*) * d2);
            for( j = 0 ; j < d2 ; j++ ){
                p[i][j] = (double*)malloc(sizeof(double) * d3);
                for( k = 0 ; k < d3 ; k++ ){
                    p[i][j][k] = other.p[i][j][k];
                }
            }
        }
    } else {
        p = NULL;
    }
}

Ten3D::~Ten3D(){
    if( p != NULL ){
        int i, j;
        for( i = 0 ; i < d1 ; i++ ){
            for( j = 0 ; j < d2 ; j++ ){
                free(p[i][j]);
            }
            free(p[i]);
        }
        free(p);
    }
}


//
