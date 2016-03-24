/*
 *  vbHmm_template.c
 *  Template for original model to be analyzed by VB-HMM.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include <math.h>
#include <string.h>

// Special functions is probably necessary.
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "vbHmm_template.h"

// For random number generation.
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))


static int isGlobalAnalysis = 0;

void setFunctions__(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters__;
    funcs.freeModelParameters   = freeModelParameters__;
    funcs.newModelStats         = newModelStats__;
    funcs.freeModelStats        = freeModelStats__;
    funcs.initializeVbHmm       = initializeVbHmm__;
    funcs.pTilde_z1             = pTilde_z1__;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1__;
    funcs.pTilde_xn_zn          = pTilde_xn_zn__;
    funcs.calcStatsVars         = calcStatsVars__;
    funcs.maximization          = maximization__;
    funcs.varLowerBound         = varLowerBound__;
    funcs.reorderParameters     = reorderParameters__;
    funcs.outputResults         = outputResults__;
    setFunctions( funcs );
}

void setGFunctions__(){
    gCommonFunctions funcs;
    funcs.newModelParameters    = newModelParameters__;
    funcs.freeModelParameters   = freeModelParameters__;
    funcs.newModelStats         = newModelStats__;
    funcs.freeModelStats        = freeModelStats__;
    funcs.newModelStatsG        = newModelStatsG__;
    funcs.freeModelStatsG       = freeModelStatsG__;
    funcs.initializeVbHmmG      = initializeVbHmmG__;
    funcs.pTilde_z1             = pTilde_z1__;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1__;
    funcs.pTilde_xn_zn          = pTilde_xn_zn__;
    funcs.calcStatsVarsG        = calcStatsVarsG__;
    funcs.maximizationG         = maximizationG__;
    funcs.varLowerBoundG        = varLowerBoundG__;
    funcs.reorderParametersG    = reorderParametersG__;
    funcs.outputResultsG        = outputResultsG__;
    setGFunctions( funcs );
    isGlobalAnalysis = 1;
}


// Redirect output function to that for model-specific data.
void outputResults__( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    output__Results( xn, gv, iv, logFP );
}

void outputResultsG__( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    output__ResultsG( xns, gv, ivs, logFP );
}


void *newModelParameters__( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    tempParameters *p = (void*)malloc( sizeof(tempParameters) );
    
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

void freeModelParameters__( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    tempParameters *gp = *p;
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
    
    free( *p );
    *p = NULL;
}


void *newModelStats__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tempStats *s = (tempStats*)malloc( sizeof(tempStats) );
        
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

void freeModelStats__( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        tempStats *gs = *s;
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

void *newModelStatsG__( xns, gv, ivs)
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    tempGlobalStats *gs = (tempGlobalStats*)malloc( sizeof(tempGlobalStats) );
    
    int i;
    gs->NiR = (double *)malloc( sNo * sizeof(double) );
    gs->NijR = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   gs->NijR[i] = (double *)malloc( sNo * sizeof(double) );   }
    gs->NiiR = (double *)malloc( sNo * sizeof(double) );
    gs->z1iR = (double *)malloc( sNo * sizeof(double) );
    
    return gs;
}

void freeModelStatsG__( gs, xns, gv, ivs )
void **gs;
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    int sNo = gv->sNo;
    tempGlobalStats *ggs = *gs;
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


void initializeVbHmm__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // an example of hyperparameter initialization

    int sNo = gv->sNo;
    tempParameters *p = gv->params;

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
    
    initialize_indVars__( xn, gv, iv );
    
    calcStatsVars__( xn, gv, iv );
    maximization__( xn, gv, iv );
}

void initializeVbHmmG__( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    // an example of hyperparameter initialization

    int sNo = gv->sNo, rNo = xns->R;
    tempParameters *p = gv->params;

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
    
    for( r = 0 ; r < rNo ; r++ ){
        initialize_indVars__( xns->xn[r], gv, ivs->indVars[r] );
    }
    
    calcStatsVarsG__( xns, gv, ivs );
    maximizationG__( xns, gv, ivs );
}


void initialize_indVars__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // an example of latent variable distribution
    
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


xnDataSet *newXnDataSet__( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (tempData*)malloc( sizeof(tempData) );
    tempData *d = (tempData*)xn->data;
    return xn;
}

void freeXnDataSet__( xn )
xnDataSet **xn;
{
    tempData *d = (tempData*)(*xn)->data;
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1__( i, params )
int i;
void *params;
{
    // Return the value for p(z_{1i}|pi_i) for given i.

    // for common HMM
    tempParameters *p = (tempParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1__( i, j, params )
int i, j;
void *params;
{
    // Return the value for p(z_{nj}|z_{n-1,i},theta) for given i and j.

    // for common HMM
    tempParameters *p = (tempParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn__( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    // Return the value for p(x_{mi}|z_{ni},phi) for given i.

    // return emission probability for the model.
    tempParameters *p = (tempParameters*)params;
    return 0.0;
}


void calcStatsVars__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // Execute calculation of statistics, which may be necessary for M-step.

    // calculation of stats for common HMM
    tempStats *s = (tempStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Nii = s->Nii, **Nij = s->Nij, *Ni = s->Ni;
    size_t n;
    int i, j;
    
    for( i = 0 ; i < sNo ; i++ ){
        Ni[i]   = 1e-10;
        Nii[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }
        
        for( n = 0 ; n < dLen ; n++ ){
            Ni[i]   += gmMat[n][i];
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
    }
}

void calcStatsVarsG__( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    // Execute calculation of statistics for global analysis, which may be necessary for M-step.
    
    // calculation of stats for global analysis of common HMM
    tempGlobalStats *gs = (tempGlobalStats*)ivs->stats;
    int sNo = gv->sNo, rNo = xns->R;
    double **gmMat, ***xiMat;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *NiR = gs->NiR, *z1iR = gs->z1iR;
    size_t dLen, n;
    int i, j, r;
    
    for( i = 0 ; i < sNo ; i++ ){
        NiR[i]   = 1e-10;
        NiiR[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijR[i][j] = 1e-10;
        }
        z1iR[i] = 1e-10;
    }
    for( r = 0 ; r < rNo ; r++ ){
        dLen = xns->xn[r]->N;
        gmMat = ivs->indVars[r]->gmMat;
        xiMat = ivs->indVars[r]->xiMat;
        for( i = 0 ; i < sNo ; i++ ){
            z1iR[i] += gmMat[0][i];
            for( n = 0 ; n < dLen ; n++ ){
                NiR[i]   += gmMat[n][i];
                for( j = 0 ; j < sNo ; j++ ){
                    NiiR[i]  += xiMat[n][i][j];
                    NijR[i][j] += xiMat[n][i][j];
                }
            }
        }
    }
}


void maximization__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // Calculation of M-step, depending on the model.

    // calculation of M-step for common HMM
    tempParameters *p = (tempParameters*)gv->params;
    tempStats *s = (tempStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *Nii = s->Nii, **Nij = s->Nij;
    int i, j;
    
    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );
        
        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Nii[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Nii[i] );
        }
    }
}

void maximizationG__( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    // Calculation of M-step for global analysis, depending on the model.
    
    // calculation of stats for global analysis of common HMM
    tempParameters *p = (tempParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    tempGlobalStats *gs = (tempGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *z1iR = gs->z1iR, dR = (double)(xns->R);
    int sNo = gv->sNo;
    int i, j;
    
    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + z1iR[i] ) / ( sumUPi + dR );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + z1iR[i] ) - gsl_sf_psi( sumUPi + dR );
        
        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + NijR[i][j] ) / ( sumUAArr[i] + NiiR[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + NijR[i][j] ) - gsl_sf_psi( sumUAArr[i] + NiiR[i] );
        }
    }
}


double varLowerBound__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // Calculation of lower bound, depending on the model.

    // calculation of lower bound of common HMM
    tempParameters *p = (tempParameters*)gv->params;
    tempStats *s = (tempStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *Nii = s->Nii, **Nij = s->Nij;
    size_t n;
    int i, j;
    
    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpPhi = 0.0;            // ln p(phi) for model-specific parameters
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqPhi = 0.0;            // ln q(phi) for model-specific parameters
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
    val  = lnpPi + lnpA + lnpPhi;
    val -= lnqPi + lnqA + lnqPhi;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));
    
    return val;
}

double varLowerBoundG__( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    // Calculation of lower bound for global analysis, depending on the model.
    
    // calculation of lower bound for global analysis of common HMM
    int sNo = gv->sNo, rNo = xns->R;
    tempParameters *p = (tempParameters*)gv->params;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    tempGlobalStats *gs = (tempGlobalStats*)ivs->stats;
    double *NiiR = gs->NiiR, **NijR = gs->NijR, *z1iR = gs->z1iR, dR = (double)xns->R;
    size_t n;
    int i, j, r;
    
    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpPhi = 0.0;            // ln p(phi) for model-specific parameters
    double lnqPi = gsl_sf_lngamma(sumUPi + dR);
    double lnqA = 0.0;
    double lnqPhi = 0.0;            // ln q(phi) for model-specific parameters
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
    val  = lnpPi + lnpA + lnpPhi;
    val -= lnqPi + lnqA + lnqPhi;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));
    
    return val;
}    


void reorderParameters__( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    // Reorder states in the final result, for example, in the order of intensity.

    // an example for common HMM.
    tempParameters *p = (tempParameters*)gv->params;
    tempStats *s = (tempStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *Ni = s->Ni;
    size_t n;
    int i, j;
    
    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }
    
    // index indicates order of avgPi values (0=biggest avgMu -- sNo=smallest avgMu).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgPi[i] < avgPi[j] ){
                    index[i]--;
                } else if( avgPi[i] == avgPi[j] ){
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

void reorderParametersG__( xns, gv, ivs )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
{
    // Reorder states in the final result, for example, in the order of intensity.
    
    // an example for global analysis of common HMM.
    size_t dLen;
    int sNo = gv->sNo, rNo = xns->R;
    tempParameters *p = (tempParameters*)gv->params;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    
    size_t n;
    int i, j, r;
    
    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( MAX(sNo,2) * sizeof(double) );   }
    
    // index indicates order of avgPi values (0=biggest avgMu -- sNo=smallest avgMu).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgPi[i] < avgPi[j] ){
                    index[i]--;
                } else if( avgPi[i] == avgPi[j] ){
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
    
    double *NiR = ((tempGlobalStats*)ivs->stats)->NiR;
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


void output__Results( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    // Output the final result (after reordering) in any format, e.g. stdout, stderr or files.

    // an example for common HMM parameters.
    tempParameters *p = (tempParameters*)gv->params;
    int sNo = gv->sNo;
    
    int i, j;
    fprintf(logFP, "  results: K = %d \n", sNo);
    
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
        fprintf(fp, "pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g", p->avgPi[i]);
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

void output__ResultsG( xns, gv, ivs, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
FILE *logFP;
{
    // Output the final result (after reordering) in any format, e.g. stdout, stderr or files.
    
    // an example for parameters of global analysis of common HMM.
    int sNo = gv->sNo, rNo = xns->R;
    tempParameters *p = (tempParameters*)gv->params;
    int i, j, r;
    
    fprintf(logFP, "  results: K = %d \n", sNo);
    
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
        fprintf(fp, "pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g", p->avgPi[i]);
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
    
    for( r = 0 ; r < rNo ; r++ ){
        xnDataSet *xn = xns->xn[r];
        indVars *iv = ivs->indVars[r];
        
        sprintf( fn, "%s.maxS%03d", xn->name, sNo );
        if( (fp = fopen( fn, "w")) != NULL ){
            for( n = 0 ; n < xn->N ; n++ ){
                fprintf( fp, "%d\n", iv->stateTraj[n] );
            }
            fclose(fp);
        }
    }
}

//
