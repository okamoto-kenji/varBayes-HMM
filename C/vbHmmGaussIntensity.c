/*
 *  vbHmmGaussIntensity.c
 *  Model-specific core functions for VB-HMM-GAUSS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

#include "vbHmmGaussIntensity.h"
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

void setFunctions_gaussInt(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_gaussInt;
    funcs.freeModelParameters   = freeModelParameters_gaussInt;
    funcs.newModelStats         = newModelStats_gaussInt;
    funcs.freeModelStats        = freeModelStats_gaussInt;
    funcs.initializeVbHmm       = initializeVbHmm_gaussInt;
    funcs.pTilde_z1             = pTilde_z1_gaussInt;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_gaussInt;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_gaussInt;
    funcs.calcStatsVars         = calcStatsVars_gaussInt;
    funcs.maximization          = maximization_gaussInt;
    funcs.varLowerBound         = varLowerBound_gaussInt;
    funcs.reorderParameters     = reorderParameters_gaussInt;
    funcs.outputResults         = outputResults_gaussInt;
    setFunctions( funcs );
}

void setGFunctions_gaussInt(){
    gCommonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_gaussInt;
    funcs.freeModelParameters   = freeModelParameters_gaussInt;
    funcs.newModelStats         = newModelStats_gaussInt;
    funcs.freeModelStats        = freeModelStats_gaussInt;
    funcs.newModelStatsG        = newModelStatsG_gaussInt;
    funcs.freeModelStatsG       = freeModelStatsG_gaussInt;
    funcs.initializeVbHmmG      = initializeVbHmmG_gaussInt;
    funcs.pTilde_z1             = pTilde_z1_gaussInt;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_gaussInt;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_gaussInt;
    funcs.calcStatsVarsG        = calcStatsVarsG_gaussInt;
    funcs.maximizationG         = maximizationG_gaussInt;
    funcs.varLowerBoundG        = varLowerBoundG_gaussInt;
    funcs.reorderParametersG    = reorderParametersG_gaussInt;
    funcs.outputResultsG        = outputResultsG_gaussInt;
    setGFunctions( funcs );
    isGlobalAnalysis = 1;
}



void outputResults_gaussInt( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputGaussIntResults( xn, gv, iv, logFP );
}

void outputResultsG_gaussInt( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    outputGaussIntResultsG( xns, gv, ivs, logFP );
}


void *newModelParameters_gaussInt( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    gaussIntParameters *p = (gaussIntParameters*)malloc( sizeof(gaussIntParameters) );
    
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

    return p;
}



void freeModelParameters_gaussInt( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    gaussIntParameters *gp = *p;
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
    
    free( gp );
    *p = NULL;
}


void *newModelStats_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussIntStats *s = (gaussIntStats*)malloc( sizeof(gaussIntStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
        s->Nii = (double *)malloc( sNo * sizeof(double) );

        return s;
        
    } else {
        
        return NULL;
        
    }
}

void freeModelStats_gaussInt( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussIntStats *gs = *s;
        int i;

        free( gs->Ni );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Nij[i] );   }
        free( gs->Nij );
        free( gs->Nii );

        free( gs );
        *s = NULL;
    }
}

void *newModelStatsG_gaussInt( xns, gv, ivs)
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussIntGlobalStats *gs = (gaussIntGlobalStats*)malloc( sizeof(gaussIntGlobalStats) );
    
    int i;
    gs->NijR = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   gs->NijR[i] = (double *)malloc( sNo * sizeof(double) );   }
    gs->NiiR = (double *)malloc( sNo * sizeof(double) );
    gs->NiR = (double *)malloc( sNo * sizeof(double) );
    gs->z1iR = (double *)malloc( sNo * sizeof(double) );
    
    return gs;
}

void freeModelStatsG_gaussInt( gs, xns, gv, ivs )
void **gs;
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussIntGlobalStats *ggs = *gs;
    int i;
    free( ggs->NiR );
    for( i = 0 ; i < sNo ; i++ )
    {   free( ggs->NijR[i] );   }
    free( ggs->NijR );
    free( ggs->NiiR );
    free( ggs->z1iR );
    
    free( *gs );
    *gs = NULL;
}


void initializeVbHmm_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussIntData *d = xn->data;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    gaussIntParameters *p = gv->params;

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
        p->uBt = 1.0;
        p->uMu = meanX * 2.0 / (double)sNo;
        p->uA  = 1.0;
        p->uB  = p->uA / precX ;
    }
    
    initialize_indVars_gaussInt( xn, gv, iv );
    
    calcStatsVars_gaussInt( xn, gv, iv );
    maximization_gaussInt( xn, gv, iv );

}

void initializeVbHmmG_gaussInt( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussIntData *d;
    size_t dLen, totalN;
    int sNo = gv->sNo, rNo = xns->R;
    gaussIntParameters *p = gv->params;
    
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
        p->uBt = 1.0;
        p->uMu = meanX * 2.0 / (double)sNo;
        p->uA  = 1.0;
        p->uB  = p->uA / precX ;
    }
    
    for( r = 0 ; r < rNo ; r++ ){
        initialize_indVars_gaussInt( xns->xn[r], gv, ivs->indVars[r] );
    }
    
    calcStatsVarsG_gaussInt( xns, gv, ivs );
    maximizationG_gaussInt( xns, gv, ivs );
}


void initialize_indVars_gaussInt( xn, gv, iv )
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


xnDataSet *newXnDataSet_gaussInt( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (gaussIntData*)malloc( sizeof(gaussIntData) );
    gaussIntData *d = (gaussIntData*)xn->data;
    d->v = NULL;
    return xn;
}

void freeXnDataSet_gaussInt( xn )
xnDataSet **xn;
{
    gaussIntData *d = (gaussIntData*)(*xn)->data;
    free( d->v );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_gaussInt( i, params )
int i;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_gaussInt( i, j, params )
int i, j;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_gaussInt( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    gaussIntData *xn = (gaussIntData*)xnWv->data;
    double val, di = (double)(i+1);
    val  = p->avgLnLm - log(2.0 * M_PI * di);
    val -= di/p->btMu + p->aLm * pow(xn->v[n] - di * p->mu0, 2.0) / di / p->bLm;
    return exp(val / 2.0);
}


void calcStatsVars_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussIntData *d = (gaussIntData*)xn->data;
    gaussIntStats *s = (gaussIntStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Nii = s->Nii, **Nij = s->Nij, *Ni = s->Ni;
    size_t n;
    int i, j;

    s->N0 = (double)xn->N;
    s->N0i = 1e-10;
    s->Nx  = 1e-10;
    for( i = 0 ; i < sNo ; i++ ){
        Ni[i]   = 1e-10;
        Nii[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }

        for( n = 0 ; n < dLen ; n++ ){
            Ni[i] += gmMat[n][i];
            s->Nx += gmMat[n][i] * d->v[n];
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        s->N0i += (double)(i+1) * Ni[i];
    }
}

void calcStatsVarsG_gaussInt( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussIntData *d;
    gaussIntGlobalStats *gs = (gaussIntGlobalStats*)ivs->stats;
    int sNo = gv->sNo, rNo = xns->R;
    double **gmMat, ***xiMat;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *NiR = gs->NiR, *z1iR = gs->z1iR;
    size_t dLen, n;
    int i, j, r;

    gs->N0R = 0.0;
    gs->N0iR = 1e-10;
    gs->NxR  = 1e-10;
    for( i = 0 ; i < sNo ; i++ ){
        NiR[i]   = 1e-10;
        NiiR[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijR[i][j] = 1e-10;
        }
        z1iR[i] = 1e-10;
    }
    for( r = 0 ; r < rNo ; r++ ){
        d = (gaussIntData*)xns->xn[r]->data;
        dLen = xns->xn[r]->N;
        gmMat = ivs->indVars[r]->gmMat;
        xiMat = ivs->indVars[r]->xiMat;
        for( i = 0 ; i < sNo ; i++ ){
            z1iR[i] += gmMat[0][i];
            for( n = 0 ; n < dLen ; n++ ){
                NiR[i]  += gmMat[n][i];
                gs->NxR += gmMat[n][i] * d->v[n];
                for( j = 0 ; j < sNo ; j++ ){
                    NiiR[i]    += xiMat[n][i][j];
                    NijR[i][j] += xiMat[n][i][j];
                }
            }
        }
        gs->N0R += (double)xns->xn[r]->N;
    }
    for( i = 0 ; i < sNo ; i++ ){
        gs->N0iR += (double)(i+1) * NiR[i];
    }
}


void maximization_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    gaussIntStats *s = (gaussIntStats*)iv->stats;
    gaussIntData *d = (gaussIntData*)xn->data;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *Nii = s->Nii, **Nij = s->Nij, N0 = s->N0, N0i = s->N0i, Nx = s->Nx;
    int i, j, m;

    double bLm = 0.0, di;
    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );

        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Nii[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Nii[i] );
        }

        di = (double)(i+1);
        for( m = 0 ; m < dLen ; m++ ){
            bLm += gmMat[m][i] * d->v[m] * d->v[m] / di;
        }
    }
    bLm /= 2.0;

    p->btMu = uBt + N0i;
    p->mu0  = (uBt * uMu + Nx) / p->btMu;
    p->aLm  = uA + (N0 / 2.0);
    bLm += uB + (uBt * uMu * uMu / 2.0);
    bLm -= p->btMu * p->mu0 * p->mu0 / 2.0;
    p->bLm = bLm;

    p->avgMu    = p->mu0;
    p->avgLm    = p->aLm / p->bLm;
    p->avgLnLm  = gsl_sf_psi( p->aLm ) - log( p->bLm );
}

void maximizationG_gaussInt( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    gaussIntGlobalStats *gs = (gaussIntGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, N0R = gs->N0R, N0iR = gs->N0iR, NxR = gs->NxR;
    double *z1iR = gs->z1iR, dR = (double)(xns->R), di, bLm = 0.0;
    int sNo = gv->sNo, rNo = xns->R;
    int i, j, m, r;
    
    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + z1iR[i] ) / ( sumUPi + dR );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + z1iR[i] ) - gsl_sf_psi( sumUPi + dR );
        
        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + NijR[i][j] ) / ( sumUAArr[i] + NiiR[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + NijR[i][j] ) - gsl_sf_psi( sumUAArr[i] + NiiR[i] );
        }

        di = (double)(i+1);
        for( r = 0 ; r < rNo ; r++ ){
            size_t dLen = xns->xn[r]->N;
            gaussIntData *d = (gaussIntData*)xns->xn[r]->data;
            double **gmMat = ivs->indVars[r]->gmMat;
            for( m = 0 ; m < dLen ; m++ ){
                bLm += gmMat[m][i] * d->v[m] * d->v[m] / di;
            }
        }

    }
    bLm /= 2.0;

    p->btMu = uBt + N0iR;
    p->mu0  = (uBt * uMu + NxR) / p->btMu;
    p->aLm  = uA + N0R / 2.0;
    bLm += uB + uBt * uMu * uMu / 2.0;
    bLm -= p->btMu * p->mu0 * p->mu0 / 2.0;
    p->bLm = bLm;
    
    p->avgMu    = p->mu0;
    p->avgLm    = p->aLm / p->bLm;
    p->avgLnLm  = gsl_sf_psi( p->aLm ) - log( p->bLm );
}


double varLowerBound_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double avgLm = p->avgLm, avgLnLm = p->avgLnLm;
    double mu0 = p->mu0, btMu = p->btMu, aLm = p->aLm, bLm = p->bLm;
    gaussIntStats *s = (gaussIntStats*)iv->stats;
    double *Nii = s->Nii, **Nij = s->Nij;
    size_t n;
    int i, j;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpMuLm = log(uBt) / 2.0;
    lnpMuLm -= uBt * (1.0 / btMu + aLm * pow(mu0 - uMu, 2.0) / bLm ) / 2.0;
    lnpMuLm += - gsl_sf_lngamma(uA) + uA * log(uB);
    lnpMuLm += (uA - 0.5) * avgLnLm - uB * avgLm;

    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqMuLm = (log(btMu) - 1.0) / 2.0 - gsl_sf_lngamma(aLm) + aLm * log(bLm);
    lnqMuLm += (aLm - 0.5) * avgLnLm - aLm;

    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        
        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Nii[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi(sumUAArr[i]+Nii[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }
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

double varLowerBoundG_gaussInt( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double avgLm = p->avgLm, avgLnLm = p->avgLnLm;
    double mu0 = p->mu0, btMu = p->btMu, aLm = p->aLm, bLm = p->bLm;
    gaussIntGlobalStats *gs = (gaussIntGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *z1iR = gs->z1iR, dR = (double)xns->R;
    size_t n;
    int sNo = gv->sNo, rNo = xns->R;
    int i, j, r;
    
    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpMuLm = log(uBt) / 2.0;
    lnpMuLm -= uBt * (1.0 / btMu + aLm * pow(mu0 - uMu, 2.0) / bLm ) / 2.0;
    lnpMuLm += - gsl_sf_lngamma(uA) + uA * log(uB);
    lnpMuLm += (uA - 0.5) * avgLnLm - uB * avgLm;

    double lnqPi = gsl_sf_lngamma(sumUPi + dR);
    double lnqA = 0.0;
    double lnqMuLm = (log(btMu) - 1.0) / 2.0 - gsl_sf_lngamma(aLm) + aLm * log(bLm);
    lnqMuLm += (aLm - 0.5) * avgLnLm - aLm;

    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnqPi += (uPiArr[i]+z1iR[i]-1.0) * (gsl_sf_psi(uPiArr[i]+z1iR[i]) - gsl_sf_psi(sumUPi+dR));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + z1iR[i]);
        
        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + NiiR[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);
            
            lnqA += (uAMat[i][j] + NijR[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+NijR[i][j]) - gsl_sf_psi(sumUAArr[i]+NiiR[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + NijR[i][j] );
        }
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


void reorderParameters_gaussInt( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
}

void reorderParametersG_gaussInt( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
}


void outputGaussIntResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    int sNo = gv->sNo;

    int i, j;
    fprintf(logFP, "  results: K = %d \n", sNo);

    fprintf(logFP, "   mean: ( %g ) \n", p->avgMu);

    fprintf(logFP, "   lambda: ( %g ) \n", p->avgLm);

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
            fprintf(fp, "%g, %g, %g", (double)(i+1)*p->avgMu, (double)(i+1)*p->avgLm, p->avgPi[i]);
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
    
    sprintf( fn, "%s.maxS%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < xn->N ; n++ ){
            fprintf( fp, "%d\n", iv->stateTraj[n] );
        }
        fclose(fp);
    }
    
}

void outputGaussIntResultsG( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    int sNo = gv->sNo, rNo = xns->R;
    gaussIntParameters *p = (gaussIntParameters*)gv->params;
    int i, j, r;
    
    fprintf(logFP, "  results: K = %d \n", sNo);
    
    fprintf(logFP, "   mean: ( %g ) \n", p->avgMu);
    
    fprintf(logFP, "   lambda: ( %g ) \n", p->avgLm);
    
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
    
    sprintf( fn, "%s.param%03d", xns->xn[0]->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "mu, lambda, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g, %g", (double)(i+1)*p->avgMu, (double)(i+1)*p->avgLm, p->avgPi[i]);
            for( j = 0 ; j < sNo ; j++ )
            {   fprintf(fp, ", %g", p->avgA[j][i]);   }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    
    sprintf( fn, "%s.Lq%03d", xns->xn[0]->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < gv->iteration ; n++ ){
            fprintf( fp, "%24.20e\n", gv->LqArr[n] );
        }
        fclose(fp);
    }
    
    sprintf( fn, "%s.maxS%03d", xns->xn[0]->name, sNo );
    int flag = 0;
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
