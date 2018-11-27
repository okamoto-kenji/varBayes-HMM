/*
 *  vbHmmModel_Common.cpp
 *  Common VB-HMM model class.
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

#include <iostream>
#include <fstream>
#include "vbHmmModel_Common.h"


//////////  model
vbHmmModel::vbHmmModel(){
    R = 0;
    dR = (double)R;
    sNo = 0;
    
    uPiArr = NULL;
    sumUPi = 0.0;
    uAMat = NULL;
    sumUAArr = NULL;
    avgPi = NULL;
    avgLnPi = NULL;
    avgA = NULL;
    avgLnA = NULL;
    
    NiiR = NULL;
    NijR = NULL;
    z1iR = NULL;
    
    iteration = 0;
    maxLq = 0.0;
    LqArr = NULL;

    ivs = NULL;
}

vbHmmModel::vbHmmModel(const vbHmmModel &other){  // do not copy
    R = 0;
    dR = (double)R;
    sNo = 0;

    uPiArr = NULL;
    sumUPi = 0.0;
    uAMat = NULL;
    sumUAArr = NULL;
    avgPi = NULL;
    avgLnPi = NULL;
    avgA = NULL;
    avgLnA = NULL;

    NiiR = NULL;
    NijR = NULL;
    z1iR = NULL;

    iteration = 0;
    maxLq = 0.0;
    LqArr = NULL;

    ivs = NULL;
}

vbHmmModel::~vbHmmModel(){
    delete uPiArr;
    delete uAMat;
    delete sumUAArr;
    delete avgPi;
    delete avgLnPi;
    delete avgA;
    delete avgLnA;

    delete NiiR;
    delete NijR;
    delete z1iR;

    delete LqArr;

    for( int i = 0 ; i < R ; i++ ){
        delete ivs[i];
        ivs[i] = NULL;
    }
    free(ivs);
}

vbHmmModel *vbHmmModel::newInstance(){
    return new vbHmmModel();
}


void vbHmmModel::setSNo( int _sNo ){
    sNo = _sNo;

    uPiArr = new Vec1D(_sNo);
    sumUPi = 0.0;
    uAMat = new Mat2D(_sNo, _sNo);
    sumUAArr = new Vec1D(_sNo);
    avgPi = new Vec1D(_sNo);
    avgLnPi = new Vec1D(_sNo);
    avgA = new Mat2D(_sNo, _sNo);
    avgLnA = new Mat2D(_sNo, _sNo);

    NiiR = new Vec1D(_sNo);
    NijR = new Mat2D(_sNo, _sNo);
    z1iR = new Vec1D(_sNo);

    LqArr = new Vec1D(0);                  // time series of lower bound
}

void vbHmmModel::appendIndVars( int _N ){
    if( R == 0 ){
        ivs = (vbHmmIndVars**)malloc( sizeof(vbHmmIndVars*) );
    } else {
        ivs =(vbHmmIndVars**)realloc( ivs, (R + 1) * sizeof(vbHmmIndVars*) );
    }
    ivs[R] = new vbHmmIndVars( _N, sNo );
    R++;
    dR = (double)R;
}


void vbHmmModel::initialize( vbHmmCond *c, vbHmmData *d ){
    int i, j;
        
    // hyper parameter for p( pi(i) )
    sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        uPiArr->p[i] = 1.0;
        sumUPi += uPiArr->p[i];
    }
    
    // hyper parameter for p( A(i,j) )
    for( i = 0 ; i < sNo ; i++ ){
        sumUAArr->p[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                uAMat->p[i][j] = 30.0;
            } else {
                uAMat->p[i][j] = 1.0;
            }
            sumUAArr->p[i] += uAMat->p[i][j];
        }
    }

    initialize_modelParams( c, d );

    switch( c->initType ){
        case initType_random:
        default:
            initialize_gmMat_random( c, d );
            break;
        case initType_kMeans:
            initialize_gmMat_kMeans( c, d );
            break;
    }

    calcStatsVars( d );
    maximization();
}

void vbHmmModel::initialize_gmMat_random( vbHmmCond *c, vbHmmData *d ){
    int i, n, r;
    double sumPar;
    for( r = 0 ; r < d->R ; r++ ){
        for( n = 0 ; n < d->traces[r]->N ; n++ ){
            sumPar = 0.0;
            for( i = 0 ; i < sNo ; i++ ){
                ivs[r]->gmMat->p[n][i] = enoise(1.0) + 1.0;
                sumPar += ivs[r]->gmMat->p[n][i];
            }
            for( i = 0 ; i < sNo ; i++ ){
                ivs[r]->gmMat->p[n][i] /= sumPar;
            }  // i
        }  // n
    }  // r
}

void vbHmmModel::initialize_gmMat_kMeans( vbHmmCond *c, vbHmmData *d ){
    int r, n, s, i;
    
    Vec1I index( d->totalN );
    Vec1I lastIndex( d->totalN );
    Vec1D sAvg( sNo );
    for( n = 0 ; n < d->totalN ; n++ ){
        index.p[n] = randomInteger(0, sNo - 1);
    }  // n
    
    double ck, dis, minDis;
    do{
        for( n = 0 ; n < d->totalN ; n++ ){
            lastIndex.p[n] = index.p[n];
        }  // n
        
        // calc avg.
        for( s = 0 ; s < sNo ; s++ ){
            ck = 0.0;
            sAvg.p[s] = 0.0;
            for( r = 0, n = 0 ; r < d->R ; r++ ){
                for( i = 0 ; i < d->traces[r]->N ; i++, n++ ){
                    if( index.p[n] == s ){
                        sAvg.p[s] += d->traces[r]->kMeans_x(i);
                        ck += 1.0;
                    }
                }  //i
            }  //r
            if( ck > 0.0 )  sAvg.p[s] /= ck;
        }  //s
        
        // assign index
        for( r = 0, n = 0 ; r < d->R ; r++ ){
            for( i = 0 ; i < d->traces[r]->N ; i++, n++ ){
                minDis = INFINITY;
                for( s = 0 ; s < sNo ; s++ ){
                    dis = d->traces[r]->kMeans_x(i) - sAvg.p[s];
                    if( dis < minDis ){
                        index.p[n] = s;
                        minDis = dis;
                    }
                }  //s
            }  //i
        }  //r
        
        // check
        for( n = 0 ; n < d->totalN ; n++ ){
            if( index.p[n] != lastIndex.p[n] ){
                break;
            }
        }  //n
    }while( n < d->totalN );
    
    for( r = 0, n = 0 ; r < d->R ; r++ ){
        for( i = 0 ; i < d->traces[r]->N ; i++, n++ ){
            double sumPar = 0.0;
            for( s = 0 ; s < sNo ; s++ ){
                if( s == index.p[n] ){
                    ivs[r]->gmMat->p[i][s] = 9.0 + enoise(1.0);
                } else {
                    ivs[r]->gmMat->p[i][s] = 1.1 + enoise(1.0);
                }
                sumPar += ivs[r]->gmMat->p[i][s];
            }  //s
            for( s = 0 ; s < sNo ; s++ ){
                ivs[r]->gmMat->p[i][s] /= sumPar;
            }  //s
        }  //i
    }  //r
    
}

void vbHmmModel::initialize_modelParams( vbHmmCond *c, vbHmmData *d ){
}

double vbHmmModel::pTilde_z1( int i ){
    return exp( avgLnPi->p[i] );
}

double vbHmmModel::pTilde_zn_zn1( int i, int j ){
    return exp( avgLnA->p[i][j] );
}

double vbHmmModel::pTilde_xn_zn( vbHmmTrace *t, int n, int i ){
    return 1.0;
}


void vbHmmModel::calcStatsVars( vbHmmData *d ){
    int i, j, n, r, N;

    double *NiiRp = NiiR->p, **NijRp = NijR->p, *z1iRp = z1iR->p;
    double **gmMat, ***xiMat;
    for( i = 0 ; i < sNo ; i++ ){
        NiiRp[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijRp[i][j] = 1e-10;
        }
        z1iRp[i] = 1e-10;
    }  // i
    for( r = 0 ; r < R ; r++ ){
        for( i = 0 ; i < sNo ; i++ ){
            N = ivs[r]->N;
            gmMat = ivs[r]->gmMat->p;
            xiMat = ivs[r]->xiMat->p;

            z1iRp[i] += gmMat[0][i];
            for( n = 0 ; n < N ; n++ ){
                for( j = 0 ; j < sNo ; j++ ){
                    NiiRp[i]    += xiMat[n][i][j];
                    NijRp[i][j] += xiMat[n][i][j];
                }  // j
            }  // n
        }  // i
    }  // r
}

int vbHmmModel::maximization(){
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi->p[i] = ( uPiArr->p[i] + z1iR->p[i] ) / ( sumUPi + dR );
        try{
            avgLnPi->p[i] = psi( uPiArr->p[i] + z1iR->p[i] ) - psi( sumUPi + dR );
        }catch (int status){
            return status;

        }

        for( j = 0 ; j < sNo ; j++ ){
            avgA->p[i][j] = ( uAMat->p[i][j] + NijR->p[i][j] ) / ( sumUAArr->p[i] + NiiR->p[i] );
            try{
                avgLnA->p[i][j] = psi( uAMat->p[i][j] + NijR->p[i][j] ) - psi( sumUAArr->p[i] + NiiR->p[i] );
            }catch (int status){
                return status;
            }
        }

        int status = maximization_modelVars(i);
        if( status != 0 ){
            return status;
        }
    }
    return 0;
}

int vbHmmModel::maximization_modelVars( int s ){
    return 0;
}


double vbHmmModel::varLowerBound( vbHmmData *d ){
    int i, j, n, r;

    double lnpPi, lnpA, lnqPi, lnqA;
    try{
        lnpPi = lnGamma(sumUPi);
        lnpA = 0.0;
        lnqPi = lnGamma(sumUPi + dR);
        lnqA = 0.0;
    }catch (int status){
        return numeric_limits<float>::quiet_NaN();
    }
    for( i = 0 ; i < sNo ; i++ ){
        try{
            lnpPi += (uPiArr->p[i] - 1.0) * avgLnPi->p[i] - lnGamma( uPiArr->p[i] );
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }

        try{
            lnqPi += (uPiArr->p[i] + z1iR->p[i] - 1.0) * ( psi( uPiArr->p[i] + z1iR->p[i] ) - psi( sumUPi + dR ));
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }
        try{
            lnqPi -= lnGamma( uPiArr->p[i] + z1iR->p[i] );
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }

        try{
            lnpA += lnGamma( sumUAArr->p[i] );
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }
        try{
            lnqA += lnGamma( sumUAArr->p[i] + NiiR->p[i] );
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }
        for( j = 0 ; j < sNo ; j++ ){
            try{
                lnpA += (uAMat->p[i][j] - 1.0) * avgLnA->p[i][j] - lnGamma( uAMat->p[i][j] );
            }catch (int status){
                return numeric_limits<float>::quiet_NaN();
            }

            try{
                lnqA += (uAMat->p[i][j] + NijR->p[i][j] - 1.0) * (psi( uAMat->p[i][j] + NijR->p[i][j] ) - psi( sumUAArr->p[i] + NiiR->p[i] ));
            }catch (int status){
                return numeric_limits<float>::quiet_NaN();
            }
            try{
                lnqA -= lnGamma( uAMat->p[i][j] + NijR->p[i][j] );
            }catch (int status){
                return numeric_limits<float>::quiet_NaN();
            }
        }  // j
    }  // i
    
    double lnpX = 0.0;
    for( r = 0 ; r < R ; r++ ){
        for( n = 0 ; n < ivs[r]->N ; n++ ){
            lnpX += log( ivs[r]->cn->p[n] );
        }
    }
    
    double val;
    val  = lnpPi + lnpA;
    val -= lnqPi + lnqA;
    double modelTerm = varLowerBound_modelTerm( d );
    if( !isfinite( modelTerm ) ){
        return numeric_limits<float>::quiet_NaN();
    }
    val += modelTerm;
    val += lnpX;
    try{
        val += log(fact(sNo));
    }catch (int status){
        return numeric_limits<float>::quiet_NaN();
   }

    return val;
}

double vbHmmModel::varLowerBound_modelTerm( vbHmmData *d ){
    return 0.0;
}


void vbHmmModel::reorderParameters( vbHmmData *d ){
}


void vbHmmModel::outputResults( vbHmmData *d, fstream *logFS ){
}


//
