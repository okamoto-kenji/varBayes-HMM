/*
 *  vbHmmPcFret.c
 *  Model-specific core functions for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmPcFret.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <string.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

//#define  DEBUG

//// Uncomment one/both of the following defenitions to activate constraint on I and K.
//#define  INTENSITY_CAP
#ifdef INTENSITY_CAP
#define  maxIntensityRatio  10.0
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

//static int isGlobalAnalysis = 0;

void setFunctions_pcFret(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_pcFret;
    funcs.freeModelParameters   = freeModelParameters_pcFret;
    funcs.newModelStats         = newModelStats_pcFret;
    funcs.freeModelStats        = freeModelStats_pcFret;
    funcs.initializeVbHmm       = initializeVbHmm_pcFret;
    funcs.pTilde_z1             = pTilde_z1_pcFret;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_pcFret;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_pcFret;
    funcs.calcStatsVars         = calcStatsVars_pcFret;
    funcs.maximization          = maximization_pcFret;
    funcs.varLowerBound         = varLowerBound_pcFret;
    funcs.reorderParameters     = reorderParameters_pcFret;
    funcs.outputResults         = outputResults_pcFret;
    setFunctions( funcs );
}

//void setGFunctions_pcFret(){
//    gCommonFunctions funcs;
//    funcs.newModelParameters    = newModelParameters_pcFret;
//    funcs.freeModelParameters   = freeModelParameters_pcFret;
//    funcs.newModelStats         = newModelStats_pcFret;
//    funcs.freeModelStats        = freeModelStats_pcFret;
//    funcs.newModelStatsG        = newModelStatsG_pcFret;
//    funcs.freeModelStatsG       = freeModelStatsG_pcFret;
//    funcs.initializeVbHmmG      = initializeVbHmmG_pcFret;
//    funcs.pTilde_z1             = pTilde_z1_pcFret;
//    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_pcFret;
//    funcs.pTilde_xn_zn          = pTilde_xn_zn_pcFret;
//    funcs.calcStatsVarsG        = calcStatsVarsG_pcFret;
//    funcs.maximizationG         = maximizationG_pcFret;
//    funcs.varLowerBoundG        = varLowerBoundG_pcFret;
//    funcs.reorderParametersG    = reorderParametersG_pcFret;
//    funcs.outputResultsG        = outputResultsG_pcFret;
//    setGFunctions( funcs );
//    isGlobalAnalysis = 1;
//}


void outputResults_pcFret( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputPcFretResults( xn, gv, iv, logFP );
}

//void outputResultsG_pcFret( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//    outputPcFretResultsG( xns, gv, ivs, logFP );
//}


void *newModelParameters_pcFret( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    pcFretParameters *p = (void*)malloc( sizeof(pcFretParameters) );
    
    p->uPiArr = (double*)malloc( sNo * sizeof(double) );
    p->sumUPi = 0.0;
    p->uAMat = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->uAMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    p->sumUAArr = (double*)malloc( sNo * sizeof(double) );

    p->aIArr = (double*)malloc( sNo * sizeof(double) );
    p->bIArr = (double*)malloc( sNo * sizeof(double) );
    p->uEArr = (double*)malloc( sNo * sizeof(double) );
    p->vEArr = (double*)malloc( sNo * sizeof(double) );
    
    p->avgPi = (double *)malloc( sNo * sizeof(double) );
    p->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    p->avgA = (double **)malloc( sNo * sizeof(double*) );
    p->avgLnA = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->avgA[i] = (double *)malloc( sNo * sizeof(double) );
        p->avgLnA[i] = (double *)malloc( sNo * sizeof(double) );
    }
    
    p->avgI = (double *)malloc( sNo * sizeof(double) );
    p->avgLnI = (double *)malloc( sNo * sizeof(double) );
    p->avgE = (double *)malloc( sNo * sizeof(double) );
    p->avgLnE = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   p->avgLnE[i] = (double *)malloc( 2 * sizeof(double) );   }
    
    return p;
}

void freeModelParameters_pcFret( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    pcFretParameters *gp = *p;
    int i;
    
    free( gp->uPiArr );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->uAMat[i] );
    }
    free( gp->uAMat );
    free( gp->sumUAArr );

    free( gp->aIArr );
    free( gp->bIArr );
    free( gp->uEArr );
    free( gp->vEArr );

    
    free( gp->avgPi );
    free( gp->avgLnPi );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgA[i] );
        free( gp->avgLnA[i] );
    }
    free( gp->avgA );
    free( gp->avgLnA );
    
    free( gp->avgI );
    free( gp->avgLnI );
    free( gp->avgE );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgLnE[i] );
    }
    free( gp->avgLnE );
    
    free( *p );
    *p = NULL;
}


void *newModelStats_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
    int sNo = gv->sNo;
    pcFretStats *s = (pcFretStats*)malloc( sizeof(pcFretStats) );
    
    int i;
    s->Ni = (double *)malloc( sNo * sizeof(double) );
    s->Ci = (double *)malloc( sNo * sizeof(double) );
    s->Di = (double *)malloc( sNo * sizeof(double) );
    s->Ai = (double *)malloc( sNo * sizeof(double) );
    s->Mi = (double *)malloc( sNo * sizeof(double) );
    s->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    
    return s;
    
//    } else {
//
//        return NULL;
//
//    }
}

void freeModelStats_pcFret( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
    int sNo = gv->sNo;
    pcFretStats *gs = *s;
    int i;
    free( gs->Ni );
    free( gs->Ci );
    free( gs->Di );
    free( gs->Ai );
    free( gs->Mi );
    for( i = 0 ; i < sNo ; i++ )
    {   free( gs->Nij[i] );   }
    free( gs->Nij );

    free( gs );
    *s = NULL;
//    }
}

//void *newModelStatsG_pcFret( xns, gv, ivs)
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    pcFretGlobalStats *gs = (pcFretGlobalStats*)malloc( sizeof(pcFretGlobalStats) );
//
//    return gs;
//}

//void freeModelStatsG_pcFret( gs, xns, gv, ivs )
//void **gs;
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    pcFretGlobalStats *ggs = *gs;
//
//    free( *gs );
//    *gs = NULL;
//}


void initializeVbHmm_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcFretData *d = xn->data;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    pcFretParameters *p = gv->params;
    
    int i, j;
    size_t totalC = 0;
    for( i = 0 ; i < dLen ; i++ ){
        totalC += d->dCounts[i] + d->aCounts[i];
    }
    double meanI = (double)totalC / (double)dLen;

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
                p->uAMat[i][j] = 100.0;
            } else {
                p->uAMat[i][j] = 1.0;
            }
            p->sumUAArr[i] += p->uAMat[i][j];
        }
    }
    
    // hyper parameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->aIArr[i] = 1.0;
        p->bIArr[i] = 1.0 / meanI;
    }
    
    // hyper parameter for p( E(i) )
    for( i = 0 ; i < sNo ; i++ ){
        p->uEArr[i] = 1.0;
        p->vEArr[i] = 1.0;
    }

    initialize_indVars_pcFret( xn, gv, iv );
    
    calcStatsVars_pcFret( xn, gv, iv );
    maximization_pcFret( xn, gv, iv );
}

//void initializeVbHmmG_pcFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void initialize_indVars_pcFret( xn, gv, iv )
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


xnDataSet *newXnDataSet_pcFret( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (pcFretData*)malloc( sizeof(pcFretData) );
    pcFretData *d = (pcFretData*)xn->data;
    d->binSize = 0.0;
    d->dCounts = NULL;
    d->aCounts = NULL;
    return xn;
}

void freeXnDataSet_pcFret( xn )
xnDataSet **xn;
{
    pcFretData *d = (pcFretData*)(*xn)->data;
    free( d->dCounts );
    free( d->aCounts );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}



double pTilde_z1_pcFret( i, params )
int i;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_pcFret( i, j, params )
int i, j;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_pcFret( xn, n, i, params )
xnDataSet *xn;
size_t n;
int i;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    pcFretData *d = (pcFretData*)xn->data;
    return exp( d->dCounts[n]*p->avgLnE[i][0] + d->aCounts[n]*p->avgLnE[i][1] + (d->dCounts[n]+d->aCounts[n])*p->avgLnI[i] - p->avgI[i] ) / gsl_sf_fact( d->dCounts[n] ) / gsl_sf_fact( d->aCounts[n] );
}


void calcStatsVars_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcFretData *d = (pcFretData*)xn->data;
    pcFretStats *s = (pcFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Ni = s->Ni, *Ci = s->Ci, *Di = s->Di, *Ai = s->Ai;
    double *Mi = s->Mi, **Nij = s->Nij;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i] = 1e-10;
        Ci[i] = 1e-10;
        Di[i] = 1e-10;
        Ai[i] = 1e-10;
        Mi[i] = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }

        for( n = 0 ; n < dLen ; n++ ){
            Ni[i]  += gmMat[n][i];
            Di[i]  += gmMat[n][i] * (double)d->dCounts[n];
            Ai[i]  += gmMat[n][i] * (double)d->aCounts[n];
            for( j = 0 ; j < sNo ; j++ ){
                Mi[i]     += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        Ci[i] = Di[i] + Ai[i];
    }

//#ifdef DEBUG
//#pragma omp critical
//{
//    for( n = 0 ; n < 20 ; n++ ){
//        for( i = 0 ; i < sNo ; i++ ){
//            fprintf(logFP, "%g,", gmMat[n][i]);
//        }
//        fprintf(logFP, "; ");
//    }
//    fprintf(logFP, "\n");
//    for( i = 0 ; i < sNo ; i++ ){
//        fprintf(logFP, "Ni(%d)=%g,  ", i, Ni[i]);
//        fprintf(logFP, "Ti(%d)=%g,  ", i, Ti[i]);
//        fprintf(logFP, "Mi(%d)=%g,  ", i, Mi[i]);
//        for( j = 0 ; j < sNo ; j++ ){
//            if( j != i )
//                fprintf(logFP, "Nij(%d,%d)=%g, ", i, j, Nij[i][j]);
//        }
//        fprintf(logFP, "\n");
//    }
//}
//#endif
}

//void calcStatsVarsG_pcFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void maximization_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcFretParameters *p = (pcFretParameters*)gv->params;
    pcFretStats *s = (pcFretStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = s->Ni, *Ci = s->Ci, *Di = s->Di, *Ai = s->Ai;
    double *Mi = s->Mi, **Nij = s->Nij;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );

        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Mi[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Mi[i] );
        }

        avgI[i] = (Ci[i] + aIArr[i]) / (Ni[i] + bIArr[i]);
        avgLnI[i] = gsl_sf_psi( Ci[i] + aIArr[i] ) - log( Ni[i] + bIArr[i] );
#ifdef  INTENSITY_CAP
        size_t n, totalC = 0;
        pcFretData *pc = xnWv->data;
        for( n = 0 ; n < xnWv->N ; n++ ){
            totalC += pc->dCounts[n] + pc->aCounts[n];
        }
        double meanI = (double)totalC / (double)xnWv->N;
        avgI[i] = MIN( avgI[i], maxIntensityRatio * meanI );
        avgLnI[i] = MIN( avgLnI[i], log(maxIntensityRatio * meanI) );
#endif
        
        avgE[i] = ( uEArr[i] + Ai[i] ) / ( uEArr[i] + vEArr[i] + Ci[i] );
        // ln(1-E) for donor
        avgLnE[i][0]  = gsl_sf_psi( vEArr[i] + Di[i] ) - gsl_sf_psi( uEArr[i] + vEArr[i] + Ci[i] );
        // ln(E) for acceptor
        avgLnE[i][1]  = gsl_sf_psi( uEArr[i] + Ai[i] ) - gsl_sf_psi( uEArr[i] + vEArr[i] + Ci[i] );
    }
}

//void maximizationG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


double varLowerBound_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcFretParameters *p = (pcFretParameters*)gv->params;
    pcFretStats *s = (pcFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, **avgLnE = p->avgLnE;
    double *Ni = s->Ni, *Ci = s->Ci, *Di = s->Di, *Ai = s->Ai;
    double *Mi = s->Mi, **Nij = s->Nij;
    size_t n;
    int i, j;

    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpI = 0.0;
    double lnpE = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqI = 0.0;
    double lnqE = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpI += aIArr[i] * log(bIArr[i]) - gsl_sf_lngamma(aIArr[i]);
        lnpI += (aIArr[i] - 1.0) * avgLnI[i] - bIArr[i] * avgI[i];

        lnpE += gsl_sf_lngamma(uEArr[i]+vEArr[i]) - gsl_sf_lngamma(uEArr[i]);
        lnpE += -gsl_sf_lngamma(vEArr[i]);
        lnpE += (uEArr[i]-1.0)*avgLnE[i][1] + (vEArr[i]-1.0)*avgLnE[i][0];

        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnqI += (Ci[i] + aIArr[i]) * log(Ni[i] + bIArr[i]) - gsl_sf_lngamma(Ci[i] + aIArr[i]);
        lnqI += (Ci[i] + aIArr[i] - 1.0) * avgLnI[i] - (Ni[i] + bIArr[i]) * avgI[i];

        lnqE += gsl_sf_lngamma(uEArr[i] + vEArr[i] + Ci[i]);
        lnqE -= gsl_sf_lngamma(uEArr[i] + Ai[i]);
        lnqE -= gsl_sf_lngamma(vEArr[i] + Di[i]);
        lnqE += (uEArr[i] + Ai[i] - 1.0) * avgLnE[i][1];
        lnqE += (vEArr[i] + Di[i] - 1.0) * avgLnE[i][0];

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Mi[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0)*avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi(sumUAArr[i]+Mi[i]));
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }

    double val;
    val  = lnpPi + lnpA + lnpI + lnpE;
    val -= lnqPi + lnqA + lnqI + lnqE;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));

//#ifdef DEBUG
//#pragma omp critical
//{
//    FILE *logFP = stderr;
//    if( val > 100000 ){
//        fprintf(logFP, "  > %g; %g; %g; %g;", lnpPi, lnpA, lnpI, lnpE);
//        fprintf(logFP, " %g; %g; %g; %g; %g\n", lnqPi, lnqA, lnqI, lnqE, lnpX);
//    }
//}
//#endif

    return val;
}

//double varLowerBoundG_pcFret( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void reorderParameters_pcFret( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcFretParameters *p = (pcFretParameters*)gv->params;
    pcFretStats *s = (pcFretStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA;
    double **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *Ni = s->Ni, *Ci = s->Ci, *Di = s->Di, *Ai = s->Ai;
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

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ci[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ci[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Di[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Di[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ai[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ai[i] = store[i];   }

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

//void reorderParametersG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void outputPcFretResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    pcFretParameters *p = (pcFretParameters*)gv->params;
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
        fprintf(fp, "I, E, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g", p->avgI[i], p->avgE[i]);
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

//void outputPcFretResultsG( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//}

//
