/*
 *  vbHmmGaussDiffusion.c
 *  Model-specific core functions for VB-HMM-GAUSS-DIFFUSION.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

#include "vbHmmGaussDiffusion.h"
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

void setFunctions_gaussDiff(){
    commonFunctions funcs;
    funcs.newModelParameters  = newModelParameters_gaussDiff;
    funcs.freeModelParameters = freeModelParameters_gaussDiff;
    funcs.newModelStats       = newModelStats_gaussDiff;
    funcs.freeModelStats      = freeModelStats_gaussDiff;
    funcs.initializeVbHmm     = initializeVbHmm_gaussDiff;
    funcs.pTilde_z1           = pTilde_z1_gaussDiff;
    funcs.pTilde_zn_zn1       = pTilde_zn_zn1_gaussDiff;
    funcs.pTilde_xn_zn        = pTilde_xn_zn_gaussDiff;
    funcs.calcStatsVars       = calcStatsVars_gaussDiff;
    funcs.maximization        = maximization_gaussDiff;
    funcs.varLowerBound       = varLowerBound_gaussDiff;
    funcs.reorderParameters   = reorderParameters_gaussDiff;
    funcs.outputResults       = outputResults_gaussDiff;
    setFunctions( funcs );
}

void setGFunctions_gaussDiff(){
    gCommonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_gaussDiff;
    funcs.freeModelParameters   = freeModelParameters_gaussDiff;
    funcs.newModelStats         = newModelStats_gaussDiff;
    funcs.freeModelStats        = freeModelStats_gaussDiff;
    funcs.newModelStatsG        = newModelStatsG_gaussDiff;
    funcs.freeModelStatsG       = freeModelStatsG_gaussDiff;
    funcs.initializeVbHmmG      = initializeVbHmmG_gaussDiff;
    funcs.pTilde_z1             = pTilde_z1_gaussDiff;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_gaussDiff;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_gaussDiff;
    funcs.calcStatsVarsG        = calcStatsVarsG_gaussDiff;
    funcs.maximizationG         = maximizationG_gaussDiff;
    funcs.varLowerBoundG        = varLowerBoundG_gaussDiff;
    funcs.reorderParametersG    = reorderParametersG_gaussDiff;
    funcs.outputResultsG        = outputResultsG_gaussDiff;
    setGFunctions( funcs );
    isGlobalAnalysis = 1;
}


void outputResults_gaussDiff( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputGaussDiffResults( xn, gv, iv, logFP );
}

void outputResultsG_gaussDiff( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    outputGaussDiffResultsG( xns, gv, ivs, logFP );
}


void *newModelParameters_gaussDiff( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    gaussDiffParameters *p = (gaussDiffParameters*)malloc( sizeof(gaussDiffParameters) );
    
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
    p->avgDlt = (double *)malloc( sNo * sizeof(double) );
    p->avgLnDlt = (double *)malloc( sNo * sizeof(double) );

    p->uAArr = (double *)malloc( sNo * sizeof(double) );
    p->uBArr = (double *)malloc( sNo * sizeof(double) );
    p->aDlt = (double *)malloc( sNo * sizeof(double) );
    p->bDlt = (double *)malloc( sNo * sizeof(double) );

    return p;
}

void freeModelParameters_gaussDiff( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    gaussDiffParameters *gp = *p;
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
    free( gp->avgDlt );
    free( gp->avgLnDlt );

    free( gp->uAArr );
    free( gp->uBArr );
    free( gp->aDlt );
    free( gp->bDlt );

    free( *p );
    *p = NULL;
}


void *newModelStats_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussDiffStats *s = (gaussDiffStats*)malloc( sizeof(gaussDiffStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Ri = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
        s->Nii = (double *)malloc( sNo * sizeof(double) );

        return s;

    } else {

        return NULL;

    }
}

void freeModelStats_gaussDiff( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        gaussDiffStats *gs = *s;
        int i;

        free( gs->Ni );
        free( gs->Ri );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Nij[i] );   }
        free( gs->Nij );
        free( gs->Nii );

        free( gs );
        *s = NULL;
    }
}

void *newModelStatsG_gaussDiff( xns, gv, ivs)
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussDiffGlobalStats *gs = (gaussDiffGlobalStats*)malloc( sizeof(gaussDiffGlobalStats) );
    
    int i;
    gs->NiR = (double *)malloc( sNo * sizeof(double) );
    gs->RiR = (double *)malloc( sNo * sizeof(double) );
    gs->NijR = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   gs->NijR[i] = (double *)malloc( sNo * sizeof(double) );   }
    gs->NiiR = (double *)malloc( sNo * sizeof(double) );
    gs->z1iR = (double *)malloc( sNo * sizeof(double) );

    return gs;
}

void freeModelStatsG_gaussDiff( gs, xns, gv, ivs )
void **gs;
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    gaussDiffGlobalStats *ggs = *gs;
    int i;
    free( ggs->NiR );
    for( i = 0 ; i < sNo ; i++ )
    {   free( ggs->NijR[i] );   }
    free( ggs->NijR );
    free( ggs->NiiR );
    free( ggs->RiR );
    free( ggs->z1iR );

    free( *gs );
    *gs = NULL;
}


void initializeVbHmm_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    int sNo = gv->sNo;
    gaussDiffParameters *p = gv->params;
    int i, j;
    
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
    
    // hyper parameter for p( delta(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uAArr[i] = 1.0;
        p->uBArr[i] = 0.0005;
    }
    
    initialize_indVars_gaussDiff( xn, gv, iv );
    
    calcStatsVars_gaussDiff( xn, gv, iv );
    maximization_gaussDiff( xn, gv, iv );
}

void initializeVbHmmG_gaussDiff( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo, rNo = xns->R;
    gaussDiffParameters *p = gv->params;
    int i, j, r;

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
    
    // hyper parameter for p( delta(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uAArr[i] = 1.0;
        p->uBArr[i] = 0.0005;
    }
    
    for( r = 0 ; r < rNo ; r++ ){
        initialize_indVars_gaussDiff( xns->xn[r], gv, ivs->indVars[r] );
    }
    
    calcStatsVarsG_gaussDiff( xns, gv, ivs );
    maximizationG_gaussDiff( xns, gv, ivs );
}


void initialize_indVars_gaussDiff( xn, gv, iv )
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


xnDataSet *newXnDataSet_gaussDiff( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (gaussDiffData*)malloc( sizeof(gaussDiffData) );
    gaussDiffData *d = (gaussDiffData*)xn->data;
    d->v = NULL;
    return xn;
}

void freeXnDataSet_gaussDiff( xn )
xnDataSet **xn;
{
    gaussDiffData *d = (gaussDiffData*)(*xn)->data;
    free( d->v );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_gaussDiff( i, params )
int i;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_gaussDiff( i, j, params )
int i, j;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_gaussDiff( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    gaussDiffData *xn = (gaussDiffData*)xnWv->data;
    double val;
    val  = p->avgLnDlt[i] - log(2.0);
    val -= log(xn->v[n]) + p->avgDlt[i] * pow( xn->v[n], 2.0) / 4.0;
    return exp(val);
}


void calcStatsVars_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussDiffData *d = (gaussDiffData*)xn->data;
    gaussDiffStats *s = (gaussDiffStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Ni = s->Ni, *Ri = s->Ri, *Nii = s->Nii, **Nij = s->Nij;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i]   = 1e-10;
        Ri[i]   = 1e-10;
        Nii[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }
        for( n = 0 ; n < dLen ; n++ ){
            Ni[i] += gmMat[n][i];
            Ri[i] += gmMat[n][i] * pow( d->v[n], 2.0 );
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
    }
}

void calcStatsVarsG_gaussDiff( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussDiffData *d;
    gaussDiffGlobalStats *gs = (gaussDiffGlobalStats*)ivs->stats;
    int sNo = gv->sNo, rNo = xns->R;
    double **gmMat, ***xiMat;
    double *NiR = gs->NiR, *RiR = gs->RiR, *NiiR = gs->NiiR;
    double **NijR = gs->NijR, *z1iR = gs->z1iR;
    size_t dLen, n;
    int i, j, r;
    
    for( i = 0 ; i < sNo ; i++ ){
        NiR[i]   = 1e-10;
        RiR[i]  = 1e-10;
        NiiR[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijR[i][j] = 1e-10;
        }
        z1iR[i] = 1e-10;
    }
    for( r = 0 ; r < rNo ; r++ ){
        d = (gaussDiffData*)xns->xn[r]->data;
        dLen = xns->xn[r]->N;
        gmMat = ivs->indVars[r]->gmMat;
        xiMat = ivs->indVars[r]->xiMat;
        for( i = 0 ; i < sNo ; i++ ){
            z1iR[i] += gmMat[0][i];
            for( n = 0 ; n < dLen ; n++ ){
                NiR[i] += gmMat[n][i];
                RiR[i] += gmMat[n][i] * pow( d->v[n], 2.0 );
                for( j = 0 ; j < sNo ; j++ ){
                    NiiR[i]    += xiMat[n][i][j];
                    NijR[i][j] += xiMat[n][i][j];
                }
            }
        }
    }
}


void maximization_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    gaussDiffStats *s = (gaussDiffStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr, *aDlt = p->aDlt, *bDlt = p->bDlt;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Ni = s->Ni, *Ri = s->Ri, *Nii = s->Nii, **Nij = s->Nij;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );

        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Nii[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Nii[i] );
        }

        aDlt[i] = uAArr[i] + Ni[i];
        bDlt[i] = uBArr[i] + Ri[i] / 4.0;

        avgDlt[i]   = aDlt[i] / bDlt[i];
        avgLnDlt[i] = gsl_sf_psi( aDlt[i] ) - log( bDlt[i] );
    }
}

void maximizationG_gaussDiff( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr, *aDlt = p->aDlt, *bDlt = p->bDlt;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    gaussDiffGlobalStats *gs = (gaussDiffGlobalStats*)ivs->stats;
    double *NiR = gs->NiR, *NiiR = gs->NiiR, **NijR = gs->NijR, *RiR = gs->RiR;
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
        
        aDlt[i] = uAArr[i] + NiR[i];
        bDlt[i] = uBArr[i] + RiR[i] / 4.0;
        
        avgDlt[i]   = aDlt[i] / bDlt[i];
        avgLnDlt[i] = gsl_sf_psi( aDlt[i] ) - log( bDlt[i] );
    }
}


double varLowerBound_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    gaussDiffStats *s = (gaussDiffStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Nii = s->Nii, **Nij = s->Nij;
    double *aDlt = p->aDlt, *bDlt = p->bDlt;
    size_t n;
    int i, j;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpDlt = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqDlt = - sNo / 2.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpDlt += - gsl_sf_lngamma(uAArr[i]) + uAArr[i] * log(uBArr[i]);
        lnpDlt += (uAArr[i] - 1.0) * avgLnDlt[i] - uBArr[i] * avgDlt[i];
        
        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Nii[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi(sumUAArr[i]+Nii[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }

        lnqDlt += - gsl_sf_lngamma(aDlt[i]) + aDlt[i] * log(bDlt[i]);
        lnqDlt += (aDlt[i] - 1.0) * avgLnDlt[i] - aDlt[i];
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }

    double val;
    val  = lnpPi + lnpA + lnpDlt;
    val -= lnqPi + lnqA + lnqDlt;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));

    return val;
}

double varLowerBoundG_gaussDiff( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    int sNo = gv->sNo, rNo = xns->R;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *aDlt = p->aDlt, *bDlt = p->bDlt;
    gaussDiffGlobalStats *gs = (gaussDiffGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *z1iR = gs->z1iR, dR = (double)xns->R;
    size_t n;
    int i, j, r;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpDlt = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + dR);
    double lnqA = 0.0;
    double lnqDlt = - sNo / 2.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);
        
        lnpDlt += - gsl_sf_lngamma(uAArr[i]) + uAArr[i] * log(uBArr[i]);
        lnpDlt += (uAArr[i] - 1.0) * avgLnDlt[i] - uBArr[i] * avgDlt[i];
        
        lnqPi += (uPiArr[i]+z1iR[i]-1.0) * (gsl_sf_psi(uPiArr[i]+z1iR[i]) - gsl_sf_psi(sumUPi+dR));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + z1iR[i]);
        
        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + NiiR[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);
            
            lnqA += (uAMat[i][j] + NijR[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+NijR[i][j]) - gsl_sf_psi(sumUAArr[i]+NiiR[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + NijR[i][j] );
        }
        
        lnqDlt += - gsl_sf_lngamma(aDlt[i]) + aDlt[i] * log(bDlt[i]);
        lnqDlt += (aDlt[i] - 1.0) * avgLnDlt[i] - aDlt[i];
    }
    
    double lnpX = 0.0;
    for( r = 0 ; r < rNo ; r++ ){
        size_t dLen = xns->xn[r]->N;
        for( n = 0 ; n < dLen ; n++ ){
            lnpX += log( ivs->indVars[r]->cn[n] );
        }
    }
    
    double val;
    val  = lnpPi + lnpA + lnpDlt;
    val -= lnqPi + lnqA + lnqDlt;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));
    
    return val;
}    


void reorderParameters_gaussDiff( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    gaussDiffStats *s = (gaussDiffStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Ni = s->Ni;
    size_t n;
    int i, j;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }

    // index indicates order of avgDlt values (0=biggest avgDlt -- sNo=smallest avgDlt).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgDlt[i] < avgDlt[j] ){
                    index[i]--;
                } else if( avgDlt[i] == avgDlt[j] ){
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

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgDlt[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgDlt[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnDlt[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnDlt[i] = store[i];   }

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

void reorderParametersG_gaussDiff( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    size_t dLen;
    int sNo = gv->sNo, rNo = xns->R;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
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
                if( avgDlt[i] < avgDlt[j] ){
                    index[i]--;
                } else if( avgDlt[i] == avgDlt[j] ){
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
    
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgDlt[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgDlt[i] = store[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnDlt[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnDlt[i] = store[i];   }
    
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
    
    double *NiR = ((gaussDiffGlobalStats*)ivs->stats)->NiR;
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


void outputGaussDiffResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    int sNo = gv->sNo;
    int i, j;
    fprintf(logFP, "  results: K = %d \n", sNo);

    fprintf(logFP, "   delta: ( %g", p->avgDlt[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgDlt[i]);
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
        fprintf(fp, "delta, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g", p->avgDlt[i], p->avgPi[i]);
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

void outputGaussDiffResultsG( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    int sNo = gv->sNo, rNo = xns->R;
    gaussDiffParameters *p = (gaussDiffParameters*)gv->params;
    int i, j, r;
    
    fprintf(logFP, "  results: K = %d \n", sNo);
    
    fprintf(logFP, "   delta: ( %g", p->avgDlt[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgDlt[i]);   }
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
    
    sprintf( fn, "%s.param%03d", xns->xn[0]->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "delta, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g", p->avgDlt[i], p->avgPi[i]);
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
