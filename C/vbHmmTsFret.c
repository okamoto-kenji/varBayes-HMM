/*
 *  vbHmmTsFret.c
 *  Model-specific core functions for VB-HMM-TS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmTsFret.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <string.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

//// Uncomment one/both of the following defenitions to activate constraint on I or/and K.
//#define  INTENSITY_CAP
//#define  TRANSITION_RATE_CAP
#ifdef INTENSITY_CAP
#define  maxIntensityRatio  10.0
#endif
#ifdef TRANSITION_RATE_CAP
#define  minPhotonNumPerState  10.0
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))


//static int isGlobalAnalysis = 0;

void setFunctions_tsFret(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_tsFret;
    funcs.freeModelParameters   = freeModelParameters_tsFret;
    funcs.newModelStats         = newModelStats_tsFret;
    funcs.freeModelStats        = freeModelStats_tsFret;
    funcs.initializeVbHmm       = initializeVbHmm_tsFret;
    funcs.pTilde_z1             = pTilde_z1_tsFret;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_tsFret;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_tsFret;
    funcs.calcStatsVars         = calcStatsVars_tsFret;
    funcs.maximization          = maximization_tsFret;
    funcs.varLowerBound         = varLowerBound_tsFret;
    funcs.reorderParameters     = reorderParameters_tsFret;
    funcs.outputResults         = outputResults_tsFret;
    setFunctions( funcs );
}

//void setGFunctions_tsFret(){
//    gCommonFunctions funcs;
//    funcs.newModelParameters    = newModelParameters_tsFret;
//    funcs.freeModelParameters   = freeModelParameters_tsFret;
//    funcs.newModelStats         = newModelStats_tsFret;
//    funcs.freeModelStats        = freeModelStats_tsFret;
//    funcs.newModelStatsG        = newModelStatsG_tsFret;
//    funcs.freeModelStatsG       = freeModelStatsG_tsFret;
//    funcs.initializeVbHmmG      = initializeVbHmmG_tsFret;
//    funcs.pTilde_z1             = pTilde_z1_tsFret;
//    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_tsFret;
//    funcs.pTilde_xn_zn          = pTilde_xn_zn_tsFret;
//    funcs.calcStatsVarsG        = calcStatsVarsG_tsFret;
//    funcs.maximizationG         = maximizationG_tsFret;
//    funcs.varLowerBoundG        = varLowerBoundG_tsFret;
//    funcs.reorderParametersG    = reorderParametersG_tsFret;
//    funcs.outputResultsG        = outputResultsG_tsFret;
//    setGFunctions( funcs );
//    isGlobalAnalysis = 1;
//}


void outputResults_tsFret( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputTsFretResults( xn, gv, iv, logFP );
}

//void outputResultsG_tsFret( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//    outputTsFretResultsG( xns, gv, ivs, logFP );
//}


void *newModelParameters_tsFret( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    tsFretParameters *p = (void*)malloc( sizeof(tsFretParameters) );
    
    p->uPiArr = (double*)malloc( sNo * sizeof(double) );
    p->sumUPi = 0.0;
    // hyperparameter for p( k(i,j) ) (i != j)
    p->uKMat = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->uKMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    p->sumUKArr = (double*)malloc( sNo * sizeof(double) );
    // hyperparameter for p( I(k) )
    p->aIArr = (double*)malloc( sNo * sizeof(double) );
    p->bIArr = (double*)malloc( sNo * sizeof(double) );
    // hyperparameter for p( E(i) )
    p->uEArr = (double*)malloc( sNo * sizeof(double) );
    p->vEArr = (double*)malloc( sNo * sizeof(double) );
    
    // parameters
    p->avgPi = (double *)malloc( sNo * sizeof(double) );
    p->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    p->avgI = (double *)malloc( sNo * sizeof(double) );
    p->avgLnI = (double *)malloc( sNo * sizeof(double) );
    p->avgK = (double **)malloc( sNo * sizeof(double*) );
    p->avgLnK = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->avgK[i] = (double *)malloc( sNo * sizeof(double) );
        p->avgLnK[i] = (double *)malloc( sNo * sizeof(double) );
    }
    p->avgLnKI = (double *)malloc( sNo * sizeof(double) );
    p->avgE = (double *)malloc( sNo * sizeof(double) );
    p->avgLnE = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   p->avgLnE[i] = (double *)malloc( 2 * sizeof(double) );   }
    
    return p;
}

void freeModelParameters_tsFret( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    tsFretParameters *gp = *p;
    int i;
    
    free( gp->uPiArr );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->uKMat[i] );
    }
    free( gp->uKMat );
    free( gp->sumUKArr );
    free( gp->aIArr );
    free( gp->bIArr );
    free( gp->uEArr );
    free( gp->vEArr );
    
    free( gp->avgPi );
    free( gp->avgLnPi );
    free( gp->avgI );
    free( gp->avgLnI );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgK[i] );
        free( gp->avgLnK[i] );
    }
    free( gp->avgK );
    free( gp->avgLnK );
    free( gp->avgLnKI );
    free( gp->avgE );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgLnE[i] );
    }
    free( gp->avgLnE );

    free( *p );
    *p = NULL;
}


void *newModelStats_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tsFretStats *s = (tsFretStats*)malloc( sizeof(tsFretStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Ti = (double *)malloc( sNo * sizeof(double) );
        s->Nii = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double *)malloc( sNo * sizeof(double) );
        s->Mij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Mij[i] = (double *)malloc( sNo * sizeof(double) );   }
        s->eps = (double *)malloc( sNo * sizeof(double) );

    return s;
    
//    } else {
//
//        return NULL;
//
//    }
}

void freeModelStats_tsFret( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tsFretStats *gs = *s;
        int i;
        free( gs->Ni );
        free( gs->Ti );
        free( gs->Nii );
        free( gs->Nij );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Mij[i] );   }
        free( gs->Mij );
        free( gs->eps );

        free( gs );
        *s = NULL;
//    }
}

//void *newModelStatsG_tsFret( xns, gv, ivs)
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    tsGlobalStats *gs = (tsGlobalStats*)malloc( sizeof(tsGlobalStats) );
//
//    return gs;
//}

//void freeModelStatsG_tsFret( gs, xns, gv, ivs )
//void **gs;
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//
//    free( *gs );
//    *gs = NULL;
//}


void initializeVbHmm_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsFretData *d = xn->data;
    int sNo = gv->sNo;
    tsFretParameters *p = gv->params;
    int i, j;
    

    // hyperparameter for p( pi(i) )
    p->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        p->uPiArr[i] = 1.0;
        p->sumUPi += p->uPiArr[i];
    }

    // hyperparameter for p( k(i,j) ) (i != j)
    for( i = 0 ; i < sNo ; i++ ){
        p->uKMat[i] = (double*)malloc( sNo * sizeof(double) );
        p->sumUKArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            p->uKMat[i][j] = 1.0;
            if( j != i ){
                p->sumUKArr[i] += p->uKMat[i][j];
            }
        }
    }

    double meanI = (double)xn->N / d->T;

    // hyperparameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->aIArr[i] = 1.0;
        p->bIArr[i] = 1.0 / meanI;
    }

    // hyperparameter for p( E(i) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uEArr[i] = 1.0;
        p->vEArr[i] = 1.0;
    }

    initialize_indVars_tsFret( xn, gv, iv );
    
    calcStatsVars_tsFret( xn, gv, iv );
    maximization_tsFret( xn, gv, iv );
}

//void initializeVbHmmG_tsFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void initialize_indVars_tsFret( xn, gv, iv )
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


xnDataSet *newXnDataSet_tsFret( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (tsFretData*)malloc( sizeof(tsFretData) );
    tsFretData *d = (tsFretData*)xn->data;
    d->T = 0.0;
    d->dt = NULL;
    d->time = NULL;
    d->ch = NULL;
    return xn;
}

void freeXnDataSet_tsFret( xn )
xnDataSet **xn;
{
    tsFretData *d = (tsFretData*)(*xn)->data;
    free( d->dt );
    free( d->time );
    free( d->ch );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_tsFret( i, params )
int i;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_tsFret( i, j, params )
int i, j;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    if( i == j ){
        return exp( p->avgLnI[i] - p->avgLnKI[i] );
    } else {
        return exp( p->avgLnK[i][i] + p->avgLnK[i][j] - p->avgLnKI[i] );
    }
}

double pTilde_xn_zn_tsFret( xn, n, i, params )
xnDataSet *xn;
size_t n;
int i;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    tsFretData *d = (tsFretData*)xn->data;
    return exp( p->avgLnI[i] - (p->avgI[i] * d->dt[n]) + p->avgLnE[i][d->ch[n]] );
}



void calcStatsVars_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsFretData *d = (tsFretData*)xn->data;
    tsFretStats *s = (tsFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Ni = s->Ni, *Ti = s->Ti, *eps = s->eps;
    double *Nii = s->Nii, *Nij = s->Nij, **Mij = s->Mij;
    size_t n;
    int i, j;
    
    for( i = 0 ; i < sNo ; i++ ){
        Ni[i] = 1e-10;
        Ti[i] = 1e-10;
        eps[i] = 1e-10;
        Nii[i] = 0.0;
        Nij[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            Mij[i][j] = 1e-10;
        }
    }
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            Ni[i]  += gmMat[n][i];
            Ti[i]  += gmMat[n][i] * d->dt[n];
            eps[i] += gmMat[n][i] * d->ch[n];
            Nii[i] += xiMat[n][i][i];
            for( j = 0 ; j < sNo ; j++ ){
                if( j != i ){
                    Mij[i][j] += xiMat[n][i][j];
                    Nij[i]    += xiMat[n][i][j];
                }
            }
        }
    }
    for( i = 0 ; i < sNo ; i++ ){
        Nii[i] = MAX( Nii[i], 1.0 );
        Nij[i] = MAX( Nij[i], 1.0 );
    }
}

//void calcStatsVarsG_tsFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}




void maximization_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsFretParameters *p = (tsFretParameters*)gv->params;
    tsFretStats *s = (tsFretStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *Ni = s->Ni, *Ti = s->Ti, *eps = s->eps;
    double *Nii = s->Nii, *Nij = s->Nij, **Mij = s->Mij;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );

        avgK[i][i] = (Ni[i] + aIArr[i]) * Nij[i] / (Nii[i] - 1.0) / (Ti[i] + bIArr[i]);
        avgLnK[i][i]  = gsl_sf_psi(Ni[i] + aIArr[i]) + gsl_sf_psi(Nij[i]);
        avgLnK[i][i] += -gsl_sf_psi(Nii[i]) - log(Ti[i] + bIArr[i]);
        
        avgI[i] = (Ni[i] + aIArr[i]) / (Ti[i] + bIArr[i]);
        avgLnI[i] = gsl_sf_psi(Ni[i] + aIArr[i]) - log(Ti[i] + bIArr[i]);
        
        avgLnKI[i]  = gsl_sf_psi(Ni[i] + aIArr[i]) + gsl_sf_psi(Nii[i] + Nij[i]);
        avgLnKI[i] += -gsl_sf_psi(Nii[i]) - log(Ti[i] + bIArr[i]);

#ifdef INTENSITY_CAP
        double meanI = (double)xnWv->N / xnWv->T;
        avgI[i]  = MIN( avgI[i], maxIntensityRatio * meanI );
        avgLnI[i]  = MIN( avgLnI[i], log(maxIntensityRatio * meanI) );
#endif
#ifdef TRANSITION_RATE_CAP
        avgK[i][i] = MIN( avgK[i][i], avgI[i]/minPhotonNumPerState );
        avgLnK[i][i]  = MIN( avgLnK[i][i], avgLnI[i] - log(minPhotonNumPerState) );
        avgLnKI[i]  = MIN( avgLnKI[i], avgLnI[i] + log(1.0 + 1.0/minPhotonNumPerState) );
#endif
        
        for( j = 0 ; j < sNo ; j++ ){
            if( i != j ){
                avgK[i][j] = ( uKMat[i][j] + Mij[i][j] ) / ( sumUKArr[i] + Nij[i] );
                avgLnK[i][j] = gsl_sf_psi( uKMat[i][j] + Mij[i][j] ) - gsl_sf_psi( sumUKArr[i] + Nij[i] );
            }
        }

        avgE[i] = ( uEArr[i] + eps[i] ) / ( uEArr[i] + vEArr[i] + Ni[i] );
        avgLnE[i][0]  = gsl_sf_psi( vEArr[i] + Ni[i] - eps[i] );    // ln(1-E) for donor
        avgLnE[i][0] -= gsl_sf_psi( uEArr[i] + vEArr[i] + Ni[i] );
        avgLnE[i][1]  = gsl_sf_psi( uEArr[i] + eps[i] );            // ln(E) for acceptor
        avgLnE[i][1] -= gsl_sf_psi( uEArr[i] + vEArr[i] + Ni[i] );
    }
}

//void maximizationG_tsFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}




double varLowerBound_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsFretParameters *p = (tsFretParameters*)gv->params;
    tsFretStats *s = (tsFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *avgLnPi = p->avgLnPi, **avgLnK = p->avgLnK, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *avgLnKI = p->avgLnKI, **avgLnE = p->avgLnE;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *Ni = s->Ni, *Ti = s->Ti, *eps = s->eps;
    double *Nii = s->Nii, *Nij = s->Nij, **Mij = s->Mij;
    size_t n;
    int i, j;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double Ck = 1.0, lnpKii = sNo * log(Ck);
    double lnpKij = 0.0;
    double lnpI = 0.0;
    double lnpE = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1);
    double lnqKiiI = 0.0;
    double lnqKij = 0.0;
    double lnqE = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpKii -= avgLnK[i][i];

        lnpI += aIArr[i] * log(bIArr[i]) - gsl_sf_lngamma(aIArr[i]);
        lnpI += (aIArr[i]-1.0)*avgLnI[i] - bIArr[i]*avgI[i];

        lnpE += gsl_sf_lngamma(uEArr[i]+vEArr[i]) - gsl_sf_lngamma(uEArr[i]);
        lnpE += -gsl_sf_lngamma(vEArr[i]);
        lnpE += (uEArr[i]-1.0)*avgLnE[i][1] + (vEArr[i]-1.0)*avgLnE[i][0];

        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnqKiiI += gsl_sf_lngamma(Nii[i] + Nij[i]) + (Ni[i] + aIArr[i]) * log(Ti[i] + bIArr[i]);
        lnqKiiI += -gsl_sf_lngamma(Ni[i] + aIArr[i]) -gsl_sf_lngamma(Nii[i]);
        lnqKiiI += -gsl_sf_lngamma(Nij[i]) + (Nij[i]-1.0)*avgLnK[i][i];
        lnqKiiI += (Ni[i] + Nii[i] + aIArr[i] - 1.0)*avgLnI[i];
        lnqKiiI += -(Nii[i] + Nij[i])*avgLnKI[i] - (Ti[i] + bIArr[i])*avgI[i];
        
        lnqE += gsl_sf_lngamma(uEArr[i] + vEArr[i] + Ni[i]);
        lnqE -= gsl_sf_lngamma(uEArr[i] + eps[i]);
        lnqE -= gsl_sf_lngamma(vEArr[i] + Ni[i] - eps[i]);
        lnqE += (uEArr[i] + eps[i] - 1.0)*avgLnE[i][1];
        lnqE += (vEArr[i] + Ni[i] - eps[i] - 1.0)*avgLnE[i][0];

        lnpKij += gsl_sf_lngamma(sumUKArr[i]);
        lnqKij += gsl_sf_lngamma(sumUKArr[i] + Nij[i]);
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                lnpKij += (uKMat[i][j]-1.0)*avgLnK[i][j] - gsl_sf_lngamma(uKMat[i][j]);
                lnqKij += (uKMat[i][j]+Mij[i][j]-1)*(gsl_sf_psi(uKMat[i][j]+Mij[i][j])-gsl_sf_psi(sumUKArr[i]+Nij[i]));
                lnqKij -= gsl_sf_lngamma( uKMat[i][j] + Mij[i][j] );
            }
        }
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }

    double val;
    val  = lnpPi + lnpKii + lnpKij + lnpI + lnpE;
    val -= lnqPi + lnqKiiI + lnqKij + lnqE;
    val += lnpX;
    
    return val;
}

//double varLowerBoundG_tsFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}




void reorderParameters_tsFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsFretParameters *p = (tsFretParameters*)gv->params;
    tsFretStats *s = (tsFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK;
    double **avgLnK = p->avgLnK, *avgLnKI = p->avgLnKI;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *Ni = s->Ni, *Ti = s->Ti, *eps = s->eps;
    size_t n;
    int i, j;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }

    // index indicates order of avgE values (0=biggest avgE -- sNo=smallest avgE).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgE[i] < avgE[j] ){
                    index[i]--;
                } else if( avgE[i] == avgE[j] ){
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

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgI[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgI[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnI[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnI[i] = store[i];   }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgK[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgK[i][j] = s2D[i][j];   }
    }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgLnK[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgLnK[i][j] = s2D[i][j];   }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnKI[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnKI[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgE[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgE[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ )
    {   s2D[index[i]][0] = avgLnE[i][0];
        s2D[index[i]][1] = avgLnE[i][1];   }
    for( i = 0 ; i < sNo ; i++ )
    {   avgLnE[i][0] = s2D[i][0];
        avgLnE[i][1] = s2D[i][1];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ni[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ni[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ti[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ti[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = eps[i];   }
    for( i = 0 ; i < sNo ; i++ ){   eps[i] = store[i];   }

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

//void reorderParametersG_tsFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}




void outputTsFretResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    tsFretParameters *p = (tsFretParameters*)gv->params;
    int sNo = gv->sNo;
    int i, j;

    fprintf(logFP, "  results: K = %d \n", sNo);

    fprintf(logFP, "   intensities: ( %g", p->avgI[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgI[i]);   }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   FRET efficiencies: ( %g", p->avgE[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgE[i]);
    }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   pi: ( %g", p->avgPi[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgPi[i]);
    }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   k_matrix: [");
    for( i = 0 ; i < sNo ; i++ ){
            fprintf(logFP, " ( %g", p->avgK[i][0]);
            for( j = 1 ; j < sNo ; j++ )
            {   fprintf(logFP, ", %g", p->avgK[i][j]);   }
            fprintf(logFP, ")");
    }
    fprintf(logFP, " ] \n\n");

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.param%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "I, E, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", K%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g, %g", p->avgI[i], p->avgE[i], p->avgPi[i]);
            for( j = 0 ; j < sNo ; j++ )
            {   fprintf(fp, ", %g", p->avgK[j][i]);   }
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

//void outputTsFretResultsG( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//}

//
