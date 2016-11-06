/*
 *  vbHmmGauss.c
 *  Model-specific core functions for VB-HMM-GAUSS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2016.11.04
 */

#include "vbHmmGauss.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <string.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

static int isGlobalAnalysis = 0;

void setFunctions_gauss(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_gauss;
    funcs.freeModelParameters   = freeModelParameters_gauss;
    funcs.newModelStats         = newModelStats_gauss;
    funcs.freeModelStats        = freeModelStats_gauss;
    funcs.initializeVbHmm       = initializeVbHmm_gauss;
    funcs.pTilde_z1             = pTilde_z1_gauss;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_gauss;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_gauss;
    funcs.calcStatsVars         = calcStatsVars_gauss;
    funcs.maximization          = maximization_gauss;
    funcs.varLowerBound         = varLowerBound_gauss;
    funcs.reorderParameters     = reorderParameters_gauss;
    funcs.outputResults         = outputResults_gauss;
    setFunctions( funcs );
}

void setGFunctions_gauss(){
    gCommonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_gauss;
    funcs.freeModelParameters   = freeModelParameters_gauss;
    funcs.newModelStats         = newModelStats_gauss;
    funcs.freeModelStats        = freeModelStats_gauss;
    funcs.newModelStatsG        = newModelStatsG_gauss;
    funcs.freeModelStatsG       = freeModelStatsG_gauss;
    funcs.initializeVbHmmG      = initializeVbHmmG_gauss;
    funcs.pTilde_z1             = pTilde_z1_gauss;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_gauss;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_gauss;
    funcs.calcStatsVarsG        = calcStatsVarsG_gauss;
    funcs.maximizationG         = maximizationG_gauss;
    funcs.varLowerBoundG        = varLowerBoundG_gauss;
    funcs.reorderParametersG    = reorderParametersG_gauss;
    funcs.outputResultsG        = outputResultsG_gauss;
    setGFunctions( funcs );
    isGlobalAnalysis = 1;
}


void outputResults_gauss( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputGaussResults( xn, gv, iv, logFP );
}

void outputResultsG_gauss( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    outputGaussResultsG( xns, gv, ivs, logFP );
}


void *newModelParameters_gauss( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    gaussParameters *p = (void*)malloc( sizeof(gaussParameters) );

    p->uPiArr = (double*)malloc( sNo * sizeof(double) );
    p->sumUPi = 0.0;
    p->uAMat = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->uAMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    p->sumUAArr = (double*)malloc( sNo * sizeof(double) );
    
    p->avgPi = (double *)malloc( sNo * sizeof(double) );
    p->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    p->avgA = (double **)malloc( sNo * sizeof(double*) );
    p->avgLnA = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->avgA[i] = (double *)malloc( sNo * sizeof(double) );
        p->avgLnA[i] = (double *)malloc( sNo * sizeof(double) );
    }
    p->avgMu = (double *)malloc( sNo * sizeof(double) );
    p->avgLm = (double *)malloc( sNo * sizeof(double) );
    p->avgLnLm = (double *)malloc( sNo * sizeof(double) );
    
    p->uBtArr = (double *)malloc( sNo * sizeof(double) );
    p->uMuArr = (double *)malloc( sNo * sizeof(double) );
    p->uAArr = (double *)malloc( sNo * sizeof(double) );
    p->uBArr = (double *)malloc( sNo * sizeof(double) );
    p->btMu = (double *)malloc( sNo * sizeof(double) );
    p->aLm = (double *)malloc( sNo * sizeof(double) );
    p->bLm = (double *)malloc( sNo * sizeof(double) );
    p->mu0 = (double *)malloc( sNo * sizeof(double) );

    return p;
}

void freeModelParameters_gauss( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    gaussParameters *gp = *p;
    int i;
    
    free( gp->uPiArr );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->uAMat[i] );
    }
    free( gp->uAMat );
    free( gp->sumUAArr );
    
    free( gp->avgPi );
    free( gp->avgLnPi );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgA[i] );
        free( gp->avgLnA[i] );
    }
    free( gp->avgA );
    free( gp->avgLnA );
    free( gp->avgMu );
    free( gp->avgLm );
    free( gp->avgLnLm );
    
    free( gp->uBtArr );
    free( gp->uMuArr );
    free( gp->uAArr );
    free( gp->uBArr );
    free( gp->btMu );
    free( gp->aLm );
    free( gp->bLm );
    free( gp->mu0 );

    free( *p );
    *p = NULL;
}


void *newModelStats_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussStats *s = (gaussStats*)malloc( sizeof(gaussStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
        s->Nii = (double *)malloc( sNo * sizeof(double) );
        s->barX = (double *)malloc( sNo * sizeof(double) );
        s->NiSi = (double *)malloc( sNo * sizeof(double) );
        
        return s;

    } else {

        return NULL;

    }
}

void freeModelStats_gauss( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussStats *gs = *s;
        int i;
        free( gs->Ni );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Nij[i] );   }
        free( gs->Nij );
        free( gs->Nii );
        free( gs->barX );
        free( gs->NiSi );

        free( gs );
        *s = NULL;
    }
}

void *newModelStatsG_gauss( xns, gv, ivs)
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussGlobalStats *gs = (gaussGlobalStats*)malloc( sizeof(gaussGlobalStats) );

    int i;
    gs->NiR = (double *)malloc( sNo * sizeof(double) );
    gs->NijR = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   gs->NijR[i] = (double *)malloc( sNo * sizeof(double) );   }
    gs->NiiR = (double *)malloc( sNo * sizeof(double) );
    gs->barXR = (double *)malloc( sNo * sizeof(double) );
    gs->NiSiR = (double *)malloc( sNo * sizeof(double) );
    gs->z1iR = (double *)malloc( sNo * sizeof(double) );

    return gs;
}

void freeModelStatsG_gauss( gs, xns, gv, ivs )
void **gs;
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussGlobalStats *ggs = *gs;
    int i;
    free( ggs->NiR );
    for( i = 0 ; i < sNo ; i++ )
    {   free( ggs->NijR[i] );   }
    free( ggs->NijR );
    free( ggs->NiiR );
    free( ggs->barXR );
    free( ggs->NiSiR );
    free( ggs->z1iR );
    
    free( *gs );
    *gs = NULL;
}


void initializeVbHmm_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussData *d = xn->data;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    gaussParameters *p = gv->params;

    int i, j;
    double totalX = 0.0, precX, varX = 0.0;
    for( i = 0 ; i < dLen ; i++ ){
        totalX += d->v[i];
    }
    double meanX = totalX / (double)dLen;
    for( i = 0 ; i < dLen ; i++ ){
        varX += pow(d->v[i] - meanX, 2.0);
    }
    varX /= (double)(dLen - 1);
    precX = 1.0 / varX;

    // hyper parameter for p( pi(i) )
    p->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        p->uPiArr[i] = 1.0;
        p->sumUPi += p->uPiArr[i];
    }

    // hyper parameter for p( A(i,j) )
    for( i = 0 ; i < sNo ; i++ ){
        p->sumUAArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                p->uAMat[i][j] = 5.0;
            } else {
                p->uAMat[i][j] = 1.0;
            }
            p->sumUAArr[i] += p->uAMat[i][j];
        }
    }
    
    // hyper parameter for p( mu(k), lm(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uBtArr[i] = 0.25;
        p->uMuArr[i] = meanX;
        p->uAArr[i] = 0.01 * varX;
        p->uBArr[i] = 0.01;
    }
    
    initialize_indVars_gauss( xn, gv, iv );

    calcStatsVars_gauss( xn, gv, iv );
    maximization_gauss( xn, gv, iv );
}

void initializeVbHmmG_gauss( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussData *d;
    size_t dLen, totalN;
    int sNo = gv->sNo, rNo = xns->R;
    gaussParameters *p = gv->params;
    
    int i, j, r;
    double totalX, precX, varX, meanX;

    totalX = 0.0;
    totalN = 0;
    for( r = 0 ; r < rNo ; r++ ){
        d = xns->xn[r]->data;
        dLen = xns->xn[r]->N;
        for( i = 0 ; i < dLen ; i++ ){
            totalX += d->v[i];
        }
        totalN += dLen;
    }
    meanX = totalX / (double)totalN;
    varX = 0.0;
    for( r = 0 ; r < rNo ; r++ ){
        d = xns->xn[r]->data;
        dLen = xns->xn[r]->N;
        for( i = 0 ; i < dLen ; i++ ){
            varX += pow(d->v[i] - meanX, 2.0);
        }
    }
    varX /= (double)(totalN - 1);
    precX = 1.0 / varX;

    p->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        p->uPiArr[i] = 1.0;
        p->sumUPi += p->uPiArr[i];
    }
    
    for( i = 0 ; i < sNo ; i++ ){
        p->sumUAArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                p->uAMat[i][j] = 5.0;
            } else {
                p->uAMat[i][j] = 1.0;
            }
            p->sumUAArr[i] += p->uAMat[i][j];
        }
    }
    
    // hyper parameter for p( mu(k), lm(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uBtArr[i] = 0.25;
        p->uMuArr[i] = meanX;
        p->uAArr[i] = 0.01 * varX;
        p->uBArr[i] = 0.01;
    }

    for( r = 0 ; r < rNo ; r++ ){
        initialize_indVars_gauss( xns->xn[r], gv, ivs->indVars[r] );
    }

    calcStatsVarsG_gauss( xns, gv, ivs );
    maximizationG_gauss( xns, gv, ivs );
}


void initialize_indVars_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    
    int i;
    size_t n;
    double sumPar;
    for( n = 0 ; n < dLen ; n++ ){
        sumPar = 0.0;
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] = enoise(1.0) + 1.0;
            sumPar += gmMat[n][i];
        }
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] /= sumPar;
        }
    }
}


xnDataSet *newXnDataSet_gauss( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (gaussData*)malloc( sizeof(gaussData) );
    gaussData *d = (gaussData*)xn->data;
    d->v = NULL;
    return xn;
}

void freeXnDataSet_gauss( xn )
xnDataSet **xn;
{
    gaussData *d = (gaussData*)(*xn)->data;
    free( d->v );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_gauss( i, params )
int i;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_gauss( i, j, params )
int i, j;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_gauss( xn, n, i, params )
xnDataSet *xn;
size_t n;
int i;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    gaussData *d = (gaussData*)xn->data;
    double val;
    val  = p->avgLnLm[i] - log(2.0 * M_PI);
    val -= 1.0/p->btMu[i] + p->aLm[i] / p->bLm[i] * pow( d->v[n] - p->avgMu[i], 2.0);
    return exp(val / 2.0);
}


void calcStatsVars_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussData *d = (gaussData*)xn->data;
    gaussStats *s = (gaussStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Nii = s->Nii, **Nij = s->Nij, *Ni = s->Ni;
    double *barX = s->barX, *NiSi = s->NiSi;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i]   = 1e-10;
        Nii[i]  = 1e-10;
        barX[i] = 1e-10;
        NiSi[i] = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }

        for( n = 0 ; n < dLen ; n++ ){
            Ni[i]   += gmMat[n][i];
            barX[i] += gmMat[n][i] * d->v[n];
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        barX[i] /= Ni[i];
        for( n = 0 ; n < dLen ; n++ ){
            NiSi[i] += gmMat[n][i] * pow( d->v[n] - barX[i], 2.0);
        }
    }
}

void calcStatsVarsG_gauss( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussData *d;
    gaussGlobalStats *gs = (gaussGlobalStats*)ivs->stats;
    int sNo = gv->sNo, rNo = xns->R;
    double **gmMat, ***xiMat;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *NiR = gs->NiR;
    double *barXR = gs->barXR, *NiSiR = gs->NiSiR, *z1iR = gs->z1iR;
    size_t dLen, n;
    int i, j, r;

    for( i = 0 ; i < sNo ; i++ ){
        NiR[i]   = 1e-10;
        NiiR[i]  = 1e-10;
        barXR[i] = 1e-10;
        NiSiR[i] = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijR[i][j] = 1e-10;
        }
        z1iR[i] = 1e-10;
    }
    for( r = 0 ; r < rNo ; r++ ){
        d = (gaussData*)xns->xn[r]->data;
        dLen = xns->xn[r]->N;
        gmMat = ivs->indVars[r]->gmMat;
        xiMat = ivs->indVars[r]->xiMat;
        for( i = 0 ; i < sNo ; i++ ){
            z1iR[i] += gmMat[0][i];
            for( n = 0 ; n < dLen ; n++ ){
                NiR[i]   += gmMat[n][i];
                barXR[i] += gmMat[n][i] * d->v[n];
                for( j = 0 ; j < sNo ; j++ ){
                    NiiR[i]  += xiMat[n][i][j];
                    NijR[i][j] += xiMat[n][i][j];
                }
            }
        }
    }
    for( i = 0 ; i < sNo ; i++ ){
        barXR[i] /= NiR[i];
        for( r = 0 ; r < rNo ; r++ ){
            d = (gaussData*)xns->xn[r]->data;
            dLen = xns->xn[r]->N;
            gmMat = ivs->indVars[r]->gmMat;
            for( n = 0 ; n < dLen ; n++ ){
                NiSiR[i] += gmMat[n][i] * pow( d->v[n] - barXR[i], 2.0);
            }
        }
    }
}


void maximization_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussParameters *p = (gaussParameters*)gv->params;
    gaussStats *s = (gaussStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    double *Ni = s->Ni, *Nii = s->Nii, **Nij = s->Nij, *barX = s->barX, *NiSi = s->NiSi;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );

        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Nii[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Nii[i] );
        }

        btMu[i] = uBtArr[i] + Ni[i];
        mu0[i]  = (uBtArr[i] * uMuArr[i] + Ni[i] * barX[i]) / btMu[i];
        aLm[i]  = uAArr[i] + Ni[i] / 2.0;
        bLm[i]  = uBArr[i] + (NiSi[i] / 2.0);
        bLm[i] += uBtArr[i] * Ni[i] * pow( barX[i] - uMuArr[i], 2.0) / 2.0 / (uBtArr[i] + Ni[i]);
        
        avgMu[i]    = mu0[i];
        avgLm[i]    = aLm[i] / bLm[i];
        avgLnLm[i]  = gsl_sf_psi( aLm[i] ) - log( bLm[i] );
    }
}

void maximizationG_gauss( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussParameters *p = (gaussParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    gaussGlobalStats *gs = (gaussGlobalStats*)ivs->stats;
    double *NiR = gs->NiR, *NiiR = gs->NiiR, **NijR = gs->NijR, *barXR = gs->barXR, *NiSiR = gs->NiSiR;
    double *z1iR = gs->z1iR, dR = (double)(xns->R);
    int sNo = gv->sNo;
    int i, j;
    
    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + z1iR[i] ) / ( sumUPi + dR );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + z1iR[i] ) - gsl_sf_psi( sumUPi + dR );
        
        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + NijR[i][j] ) / ( sumUAArr[i] + NiiR[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + NijR[i][j] ) - gsl_sf_psi( sumUAArr[i] + NiiR[i] );
        }

        btMu[i] = uBtArr[i] + NiR[i];
        mu0[i]  = (uBtArr[i] * uMuArr[i] + NiR[i] * barXR[i]) / btMu[i];
        aLm[i]  = uAArr[i] + NiR[i] / 2.0;
        bLm[i]  = uBArr[i] + (NiSiR[i] / 2.0);
        bLm[i] += uBtArr[i] * NiR[i] * pow( barXR[i] - uMuArr[i], 2.0) / 2.0 / (uBtArr[i] + NiR[i]);

        avgMu[i]    = mu0[i];
        avgLm[i]    = aLm[i] / bLm[i];
        avgLnLm[i]  = gsl_sf_psi( aLm[i] ) - log( bLm[i] );
    }
}


double varLowerBound_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussParameters *p = (gaussParameters*)gv->params;
    gaussStats *s = (gaussStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    double *Nii = s->Nii, **Nij = s->Nij;
    size_t n;
    int i, j;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpMuLm = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqMuLm = - sNo / 2.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpMuLm += log(uBtArr[i]) / 2.0;
        lnpMuLm -= uBtArr[i] / 2.0 *( 1.0 / btMu[i] + aLm[i] / bLm[i] * pow(mu0[i] - uMuArr[i], 2.0) );
        lnpMuLm += - gsl_sf_lngamma(uAArr[i]) + uAArr[i] * log(uBArr[i]);
        lnpMuLm += (uAArr[i] - 0.5) * avgLnLm[i] - uBArr[i] * avgLm[i];
        
        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Nii[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi(sumUAArr[i]+Nii[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }

        lnqMuLm += log(btMu[i]) / 2.0 - gsl_sf_lngamma(aLm[i]) + aLm[i] * log(bLm[i]);
        lnqMuLm += (aLm[i] - 0.5) * avgLnLm[i] - aLm[i];
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }

    double val;
    val  = lnpPi + lnpA + lnpMuLm;
    val -= lnqPi + lnqA + lnqMuLm;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));

    return val;
}    

double varLowerBoundG_gauss( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo, rNo = xns->R;
    gaussParameters *p = (gaussParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    gaussGlobalStats *gs = (gaussGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *z1iR = gs->z1iR, dR = (double)xns->R;
    size_t n;
    int i, j, r;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpMuLm = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + dR);
    double lnqA = 0.0;
    double lnqMuLm = - sNo / 2.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);
        
        lnpMuLm += log(uBtArr[i]) / 2.0;
        lnpMuLm -= uBtArr[i] / 2.0 *( 1.0 / btMu[i] + aLm[i] / bLm[i] * pow(mu0[i] - uMuArr[i], 2.0) );
        lnpMuLm += - gsl_sf_lngamma(uAArr[i]) + uAArr[i] * log(uBArr[i]);
        lnpMuLm += (uAArr[i] - 0.5) * avgLnLm[i] - uBArr[i] * avgLm[i];
        
        lnqPi += (uPiArr[i]+z1iR[i]-1.0) * (gsl_sf_psi(uPiArr[i]+z1iR[i]) - gsl_sf_psi(sumUPi+dR));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + z1iR[i]);
        
        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + NiiR[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);
            
            lnqA += (uAMat[i][j] + NijR[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+NijR[i][j]) - gsl_sf_psi(sumUAArr[i]+NiiR[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + NijR[i][j] );
        }
        
        lnqMuLm += log(btMu[i]) / 2.0 - gsl_sf_lngamma(aLm[i]) + aLm[i] * log(bLm[i]);
        lnqMuLm += (aLm[i] - 0.5) * avgLnLm[i] - aLm[i];
    }

    double lnpX = 0.0;
    for( r = 0 ; r < rNo ; r++ ){
        size_t dLen = xns->xn[r]->N;
        for( n = 0 ; n < dLen ; n++ ){
            lnpX += log( ivs->indVars[r]->cn[n] );
        }
    }
    
    double val;
    val  = lnpPi + lnpA + lnpMuLm;
    val -= lnqPi + lnqA + lnqMuLm;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));
    
    return val;
}    


void reorderParameters_gauss( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussParameters *p = (gaussParameters*)gv->params;
    gaussStats *s = (gaussStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    double *Ni = s->Ni;
    size_t n;
    int i, j;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }

    // index indicates order of avgMu values (0=biggest avgMu -- sNo=smallest avgMu).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgMu[i] < avgMu[j] ){
                    index[i]--;
                } else if( avgMu[i] == avgMu[j] ){
                    if( j > i )
                    {   index[i]--;   }
                }
            }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgMu[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgMu[i] = store[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLm[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnLm[i] = store[i];   }

    //
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = btMu[i];   }
    for( i = 0 ; i < sNo ; i++ ){   btMu[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = aLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   aLm[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = bLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   bLm[i] = store[i];   }
    //
    
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgA[i][j] = s2D[i][j];   }
    }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgLnA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgLnA[i][j] = s2D[i][j];   }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ni[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ni[i] = store[i];   }

    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = gmMat[n][i];   }
        for( i = 0 ; i < sNo ; i++ ){   gmMat[n][i] = store[i];   }
    }

    for( n = 0 ; n < dLen ; n++ ){
        for( j = 0 ; j < sNo ; j++ ){
            for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = xiMat[n][i][j];   }
        }
        for( j = 0 ; j < sNo ; j++ ){
            for( i = 0 ; i < sNo ; i++ ){   xiMat[n][i][j] = s2D[i][j];   }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   free( s2D[i] );   }
    free( s2D );
    free( store );
    free( index );
}

void reorderParametersG_gauss( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    size_t dLen;
    int sNo = gv->sNo, rNo = xns->R;
    gaussParameters *p = (gaussParameters*)gv->params;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
    size_t n;
    int i, j, r;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }
    
    // index indicates order of avgMu values (0=biggest avgMu -- sNo=smallest avgMu).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgMu[i] < avgMu[j] ){
                    index[i]--;
                } else if( avgMu[i] == avgMu[j] ){
                    if( j > i )
                    {   index[i]--;   }
                }
            }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgMu[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgMu[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLm[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnLm[i] = store[i];   }

    //
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = btMu[i];   }
    for( i = 0 ; i < sNo ; i++ ){   btMu[i] = store[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = aLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   aLm[i] = store[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = bLm[i];   }
    for( i = 0 ; i < sNo ; i++ ){   bLm[i] = store[i];   }
    //
    
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgA[i][j] = s2D[i][j];   }
    }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgLnA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgLnA[i][j] = s2D[i][j];   }
    }

    double *NiR = ((gaussGlobalStats*)ivs->stats)->NiR;
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = NiR[i];   }
    for( i = 0 ; i < sNo ; i++ ){   NiR[i] = store[i];   }

    for( r = 0 ; r < rNo ; r++ ){
        double **gmMat = ivs->indVars[r]->gmMat;
        dLen = xns->xn[r]->N;
        for( n = 0 ; n < dLen ; n++ ){
            for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = gmMat[n][i];   }
            for( i = 0 ; i < sNo ; i++ ){   gmMat[n][i] = store[i];   }
        }
    }

    for( r = 0 ; r < rNo ; r++ ){
        double ***xiMat = ivs->indVars[r]->xiMat;
        dLen = xns->xn[r]->N;
        for( n = 0 ; n < dLen ; n++ ){
            for( j = 0 ; j < sNo ; j++ ){
                for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = xiMat[n][i][j];   }
            }
            for( j = 0 ; j < sNo ; j++ ){
                for( i = 0 ; i < sNo ; i++ ){   xiMat[n][i][j] = s2D[i][j];   }
            }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   free( s2D[i] );   }
    free( s2D );
    free( store );
    free( index );
}


void outputGaussResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    gaussParameters *p = (gaussParameters*)gv->params;
    int sNo = gv->sNo;

    int i, j;
    fprintf(logFP, "  results: K = %d \n", sNo);

    fprintf(logFP, "   means: ( %g", p->avgMu[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgMu[i]);   }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   lambda: ( %g", p->avgLm[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgLm[i]);
    }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   pi: ( %g", p->avgPi[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgPi[i]);
    }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   A_matrix: [");
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, " ( %g", p->avgA[i][0]);
        for( j = 1 ; j < sNo ; j++ )
        {   fprintf(logFP, ", %g", p->avgA[i][j]);   }
        fprintf(logFP, ")");
    }
    fprintf(logFP, " ] \n\n");

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.param%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "mu, lambda, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g, %g", p->avgMu[i], p->avgLm[i], p->avgPi[i]);
            for( j = 0 ; j < sNo ; j++ )
            {   fprintf(fp, ", %g", p->avgA[j][i]);   }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    sprintf( fn, "%s.Lq%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < gv->iteration ; n++ ){
            fprintf( fp, "%24.20e\n", gv->LqArr[n] );
        }
        fclose(fp);
    }

#ifdef OUTPUT_MAX_GAMMA
    sprintf( fn, "%s.maxG%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < xn->N ; n++ ){
            fprintf( fp, "%d\n", iv->gammaTraj[n] );
        }
        fclose(fp);
    }
#endif

    sprintf( fn, "%s.maxS%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < xn->N ; n++ ){
            fprintf( fp, "%d\n", iv->stateTraj[n] );
        }
        fclose(fp);
    }

}

void outputGaussResultsG( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    int sNo = gv->sNo, rNo = xns->R;
    gaussParameters *p = (gaussParameters*)gv->params;
    int i, j, r;

    fprintf(logFP, "  results: K = %d \n", sNo);
    
    fprintf(logFP, "   means: ( %g", p->avgMu[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgMu[i]);   }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   lambda: ( %g", p->avgLm[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgLm[i]);
    }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   pi: ( %g", p->avgPi[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgPi[i]);
    }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   A_matrix: [");
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, " ( %g", p->avgA[i][0]);
        for( j = 1 ; j < sNo ; j++ )
        {   fprintf(logFP, ", %g", p->avgA[i][j]);   }
        fprintf(logFP, ")");
    }
    fprintf(logFP, " ] \n\n");

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.g.param%03d", xns->xn[0]->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "mu, lambda, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g, %g", p->avgMu[i], p->avgLm[i], p->avgPi[i]);
            for( j = 0 ; j < sNo ; j++ )
            {   fprintf(fp, ", %g", p->avgA[j][i]);   }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    
    sprintf( fn, "%s.g.Lq%03d", xns->xn[0]->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < gv->iteration ; n++ ){
            fprintf( fp, "%24.20e\n", gv->LqArr[n] );
        }
        fclose(fp);
    }

    int flag;
#ifdef OUTPUT_MAX_GAMMA
    sprintf( fn, "%s.g.maxG%03d", xns->xn[0]->name, sNo );
    flag = 0;
    if( (fp = fopen( fn, "w")) != NULL ){
        n = 0;
        do{
            flag = 1;
            for( r = 0 ; r < rNo ; r++ ){
                xnDataSet *xn = xns->xn[r];
                indVars *iv = ivs->indVars[r];
                
                if( r > 0 ){
                    fprintf( fp, "," );
                }
                if( n < xn->N ){
                    fprintf( fp, "%d", iv->gammaTraj[n] );
                }
                flag &= (n >= (xn->N - 1));
            }
            fprintf( fp, "\n" );
            n++;
        }while( !flag );
        fclose(fp);
    }
#endif

    sprintf( fn, "%s.g.maxS%03d", xns->xn[0]->name, sNo );
    flag = 0;
    if( (fp = fopen( fn, "w")) != NULL ){
        n = 0;
        do{
            flag = 1;
            for( r = 0 ; r < rNo ; r++ ){
                xnDataSet *xn = xns->xn[r];
                indVars *iv = ivs->indVars[r];

                if( r > 0 ){
                    fprintf( fp, "," );
                }
                if( n < xn->N ){
                    fprintf( fp, "%d", iv->stateTraj[n] );
                }
                flag &= (n >= (xn->N - 1));
            }
            fprintf( fp, "\n" );
            n++;
        }while( !flag );
        fclose(fp);
    }
}

//
