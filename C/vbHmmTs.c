/*
 *  vbHmmTs.c
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

#include "vbHmmTs.h"
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

void setFunctions_ts(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_ts;
    funcs.freeModelParameters   = freeModelParameters_ts;
    funcs.newModelStats         = newModelStats_ts;
    funcs.freeModelStats        = freeModelStats_ts;
    funcs.initializeVbHmm       = initializeVbHmm_ts;
    funcs.pTilde_z1             = pTilde_z1_ts;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_ts;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_ts;
    funcs.calcStatsVars         = calcStatsVars_ts;
    funcs.maximization          = maximization_ts;
    funcs.varLowerBound         = varLowerBound_ts;
    funcs.reorderParameters     = reorderParameters_ts;
    funcs.outputResults         = outputResults_ts;
    setFunctions( funcs );
}

//void setGFunctions_ts(){
//    gCommonFunctions funcs;
//    funcs.newModelParameters    = newModelParameters_ts;
//    funcs.freeModelParameters   = freeModelParameters_ts;
//    funcs.newModelStats         = newModelStats_ts;
//    funcs.freeModelStats        = freeModelStats_ts;
//    funcs.newModelStatsG        = newModelStatsG_ts;
//    funcs.freeModelStatsG       = freeModelStatsG_ts;
//    funcs.initializeVbHmmG      = initializeVbHmmG_ts;
//    funcs.pTilde_z1             = pTilde_z1_ts;
//    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_ts;
//    funcs.pTilde_xn_zn          = pTilde_xn_zn_ts;
//    funcs.calcStatsVarsG        = calcStatsVarsG_ts;
//    funcs.maximizationG         = maximizationG_ts;
//    funcs.varLowerBoundG        = varLowerBoundG_ts;
//    funcs.reorderParametersG    = reorderParametersG_ts;
//    funcs.outputResultsG        = outputResultsG_ts;
//    setGFunctions( funcs );
//    isGlobalAnalysis = 1;
//}


void outputResults_ts( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputTsResults( xn, gv, iv, logFP );
}

//void outputResultsG_ts( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//    outputTsResultsG( xns, gv, ivs, logFP );
//}


void *newModelParameters_ts( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    tsParameters *p = (void*)malloc( sizeof(tsParameters) );
    
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
    
    return p;
}

void freeModelParameters_ts( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    tsParameters *gp = *p;
    int i;
    
    free( gp->uPiArr );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->uKMat[i] );
    }
    free( gp->uKMat );
    free( gp->sumUKArr );
    free( gp->aIArr );
    free( gp->bIArr );
    
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
    
    free( *p );
    *p = NULL;
}


void *newModelStats_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tsStats *s = (tsStats*)malloc( sizeof(tsStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Ti = (double *)malloc( sNo * sizeof(double) );
        s->Nii = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double *)malloc( sNo * sizeof(double) );
        s->Mij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Mij[i] = (double *)malloc( sNo * sizeof(double) );   }
        
        return s;
        
//    } else {
//        
//        return NULL;
//        
//    }
}

void freeModelStats_ts( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tsStats *gs = *s;
        int i;
        free( gs->Ni );
        free( gs->Ti );
        free( gs->Nii );
        free( gs->Nij );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Mij[i] );   }
        free( gs->Mij );
    
        free( gs );
        *s = NULL;
//    }
}

//void *newModelStatsG_ts( xns, gv, ivs)
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    tsGlobalStats *gs = (tsGlobalStats*)malloc( sizeof(tsGlobalStats) );
//    
//    return gs;
//}

//void freeModelStatsG_ts( gs, xns, gv, ivs )
//void **gs;
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    
//    free( *gs );
//    *gs = NULL;
//}


void initializeVbHmm_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsData *d = xn->data;
    int sNo = gv->sNo;
    tsParameters *p = gv->params;
    
    int i, j;

    // hyperparameter for p( pi(k) )
    p->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        p->uPiArr[i] = 1.0;
        p->sumUPi += p->uPiArr[i];
    }

    // hyperparameter for p( k(i,j) ) (i != j)
    for( i = 0 ; i < sNo ; i++ ){
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
    
    initialize_indVars_ts( xn, gv, iv );
    
    calcStatsVars_ts( xn, gv, iv );
    maximization_ts( xn, gv, iv );
}

//void initializeVbHmmG_ts( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void initialize_indVars_ts( xn, gv, iv )
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


xnDataSet *newXnDataSet_ts( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (tsData*)malloc( sizeof(tsData) );
    tsData *d = (tsData*)xn->data;
    d->T = 0.0;
    d->dt = NULL;
    d->time = NULL;
    return xn;
}

void freeXnDataSet_ts( xn )
xnDataSet **xn;
{
    tsData *d = (tsData*)(*xn)->data;
    free( d->dt );
    free( d->time );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_ts( i, params )
int i;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_ts( i, j, params )
int i, j;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    if( i == j ){
        return exp( p->avgLnI[i] - p->avgLnKI[i] );
    } else {
        return exp( p->avgLnK[i][i] + p->avgLnK[i][j] - p->avgLnKI[i] );
    }
}

double pTilde_xn_zn_ts( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    tsData *d = (tsData*)xnWv->data;
    return exp( p->avgLnI[i] - p->avgI[i] * d->dt[n] );
}


void calcStatsVars_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsData *d = (tsData*)xn->data;
    tsStats *s = (tsStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Ni = s->Ni, *Ti = s->Ti;
    double *Nii = s->Nii, *Nij = s->Nij, **Mij = s->Mij;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i] = 1e-10;
        Ti[i] = 1e-10;
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

//void calcStatsVarsG_ts( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void maximization_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsParameters *p = (tsParameters*)gv->params;
    tsStats *s = (tsStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *Ni = s->Ni, *Ti = s->Ti;
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
    }
}

//void maximizationG_ts( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


double varLowerBound_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsParameters *p = (tsParameters*)gv->params;
    tsStats *s = (tsStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *avgLnPi = p->avgLnPi, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *Ni = s->Ni, *Ti = s->Ti;
    double *Nii = s->Nii, *Nij = s->Nij, **Mij = s->Mij;
    size_t n;
    int i, j;
    
    double lnpPi;
    lnpPi = gsl_sf_lngamma(sumUPi);
    double Ck = 1.0, lnpKii = (double)sNo * log(Ck);
    double lnpKij = 0.0;
    double lnpI = 0.0;
    double lnqPi;
    lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqKiiI = 0.0;
    double lnqKij = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);   

        lnpKii -= avgLnK[i][i];

        lnpI += aIArr[i] * log(bIArr[i]) - gsl_sf_lngamma(aIArr[i]);
        lnpI += (aIArr[i]-1.0)*avgLnI[i] - bIArr[i]*avgI[i];

        lnqPi += (uPiArr[i]+gmMat[0][i]-1) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnqKiiI += gsl_sf_lngamma(Nii[i] + Nij[i]) + (Ni[i] + aIArr[i]) * log(Ti[i] + bIArr[i]);
        lnqKiiI += -gsl_sf_lngamma(Ni[i] + aIArr[i]) -gsl_sf_lngamma(Nii[i]);
        lnqKiiI += -gsl_sf_lngamma(Nij[i]) + (Nij[i]-1.0)*avgLnK[i][i];
        lnqKiiI += (Ni[i] + Nii[i] + aIArr[i] - 1.0)*avgLnI[i];
        lnqKiiI += -(Nii[i] + Nij[i])*avgLnKI[i] - (Ti[i] + bIArr[i])*avgI[i];
        
        lnpKij += gsl_sf_lngamma(sumUKArr[i]);
        lnqKij += gsl_sf_lngamma(sumUKArr[i] + Nij[i]);
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                lnpKij += (uKMat[i][j]-1)*avgLnK[i][j] - gsl_sf_lngamma(uKMat[i][j]);
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
    val  = lnpPi + lnpKii + lnpKij + lnpI;
    val -= lnqPi + lnqKiiI + lnqKij;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));
    
    return val;
}

//double varLowerBoundG_ts( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void reorderParameters_ts( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    tsParameters *p = (tsParameters*)gv->params;
    tsStats *s = (tsStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK;
    double **avgLnK = p->avgLnK, *avgLnKI = p->avgLnKI;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = s->Ni, *Ti = s->Ti;
    size_t n;
    int i, j;
    
    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( sNo * sizeof(double) );   }

    // index indicates order of avgI values (0=biggest avgI -- sNo=smallest avgI).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgI[i] < avgI[j] ){
                    index[i]--;
                } else if( avgI[i] == avgI[j] ){
                    if( i < j )
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

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ni[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ni[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ti[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ti[i] = store[i];   }

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

//void reorderParametersG_ts( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void outputTsResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    tsParameters *p = (tsParameters*)gv->params;
    int sNo = gv->sNo;
    int i, j;

    fprintf(logFP, "  results: K = %d \n", sNo);
        
    fprintf(logFP, "   intensities: ( %g", p->avgI[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgI[i]);   }
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
        fprintf(fp, "I, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", K%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g", p->avgI[i], p->avgPi[i]);
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

//void outputTsResultsG( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//}

//
