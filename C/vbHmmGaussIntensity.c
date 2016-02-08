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
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

void setFunctions_gaussInt(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_gaussInt;
    funcs.initialize_vbHmm = initialize_vbHmm_gaussInt;
    funcs.freeParameters = freeParameters_gaussInt;
    funcs.pTilde_z1 = pTilde_z1_gaussInt;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_gaussInt;
    funcs.pTilde_xn_zn = pTilde_xn_zn_gaussInt;
    funcs.calcStatsVars = calcStatsVars_gaussInt;
    funcs.maximization = maximization_gaussInt;
    funcs.varLowerBound = varLowerBound_gaussInt;
    funcs.reorderParameters = reorderParameters_gaussInt;
    funcs.outputResults = outputResults_gaussInt;
    setFunctions( funcs );
}


void **mallocParameterArray_gaussInt( n )
size_t n;
{
    return (void**)malloc( n * sizeof(gaussIntParameters*) );
}


void outputResults_gaussInt( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputGaussIntResults( cParams, (gaussIntParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_gaussInt( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    gaussIntData *xn = xnWv->data;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    gaussIntParameters *p = blankParameters_gaussInt( sNo );
    p->sNo = sNo;

    int i, j;
    size_t n;
    double totalX = 0.0, precX, varX = 0.0;
    for( i = 0 ; i < xnWv->N ; i++ ){
        totalX += xn->v[i];
    }
    double meanX = totalX / (double)xnWv->N;
    for( i = 0 ; i < xnWv->N ; i++ ){
        varX += pow(xn->v[i] - meanX, 2.0);
    }
    varX /= (double)(xnWv->N - 1);
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
        p->uBt = precX / 250.0 / (double)sNo;
        p->uMu = meanX / (double)sNo;
        p->uA  = 2.5;
        p->uB  = 0.01;
    }
    
    double sumPar;
    for( n = 0 ; n < xnWv->N ; n++ ){
        sumPar = 0.0;
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] = enoise(1.0) + 1.0;
            sumPar += gmMat[n][i];
        }
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] /= sumPar;
        }
    }
    calcStatsVars_gaussInt( xnWv, cParams, p );
    maximization_gaussInt( xnWv, cParams, p );

    return p;
}

gaussIntParameters *blankParameters_gaussInt( sNo )
int sNo;
{
    int i;
    gaussIntParameters *p = (gaussIntParameters*)malloc( sizeof(gaussIntParameters) );
    
    p->sNo = 0;
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

    p->Ni = (double *)malloc( sNo * sizeof(double) );
    p->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   p->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    p->Nii = (double *)malloc( sNo * sizeof(double) );
    
    return p;
}


void freeParameters_gaussInt( params )
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    int i;

    free( p->uPiArr );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->uAMat[i] );
    }
    free( p->uAMat );
    free( p->sumUAArr );

    free( p->avgPi );
    free( p->avgLnPi );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->avgA[i] );
        free( p->avgLnA[i] );
    }
    free( p->avgA );
    free( p->avgLnA );

    free( p->Ni );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Nij[i] );   }
    free( p->Nij );
    free( p->Nii );
    
    free( p );
    p = NULL;
}

void freeGaussIntDataSet( gaussIntTraj )
xnDataSet *gaussIntTraj;
{
    gaussIntData *xn = gaussIntTraj->data;
    free( xn->v );
    free( gaussIntTraj->data );
    free( gaussIntTraj );
    gaussIntTraj = NULL;
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


void calcStatsVars_gaussInt( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    gaussIntData *xn = (gaussIntData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Nii = p->Nii, **Nij = p->Nij;
    size_t n;
    int i, j;

    p->N0  = 1e-10;
    p->N0i = 1e-10;
    p->Nx  = 1e-10;
    for( i = 0 ; i < sNo ; i++ ){
        Ni[i]   = 1e-10;
        Nii[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }

        for( n = 0 ; n < dLen ; n++ ){
            Ni[i] += gmMat[n][i];
            p->Nx += gmMat[n][i] * xn->v[n];
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        p->N0  += Ni[i];
        p->N0i += (double)(i+1) * Ni[i];
    }
}


void maximization_gaussInt( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    gaussIntData *xn = (gaussIntData*)xnWv->data;
    int sNo = cParams->sNo;
    size_t dLen = cParams->dLen;
    double **gmMat = cParams->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *Ni = p->Ni, *Nii = p->Nii, **Nij = p->Nij, N0 = p->N0, N0i = p->N0i, Nx = p->Nx;
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
            bLm += gmMat[m][i] * xn->v[m] * xn->v[m] / di;
        }
    }
    bLm /= 2.0;

    p->btMu = uBt + N0i;
    p->mu0  = (uBt * uMu + Nx) / p->btMu;
    p->aLm  = uA + N0 / 2.0;
    bLm += uB + uBt * uMu * uMu / 2.0;
    bLm -= (uBt + N0i) * p->mu0 * p->mu0 / 2.0;
    p->bLm = bLm;

    p->avgMu    = p->mu0;
    p->avgLm    = p->aLm / p->bLm;
    p->avgLnLm  = gsl_sf_psi( p->aLm ) - log( p->bLm );
}


double varLowerBound_gaussInt( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussIntParameters *p = (gaussIntParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double uBt = p->uBt, uMu = p->uMu, uA = p->uA, uB = p->uB;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double avgLm = p->avgLm, avgLnLm = p->avgLnLm;
    double *Ni = p->Ni, *Nii = p->Nii, **Nij = p->Nij, N0i = p->N0i, Nx = p->Nx;
    double mu0 = p->mu0, btMu = p->btMu, aLm = p->aLm, bLm = p->bLm;
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


void reorderParameters_gaussInt( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
}


void outputGaussIntResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
gaussIntParameters *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    int i, j;
    fprintf(logFP, "  results: K = %d \n", s);

    fprintf(logFP, "   mean: ( %g ) \n", params->avgMu);

    fprintf(logFP, "   lambda: ( %g ) \n", params->avgLm);

    fprintf(logFP, "   A_matrix: [");
    for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(logFP, " ( %g", params->avgA[i][0]);
            for( j = 1 ; j < cParams->sNo ; j++ )
            {   fprintf(logFP, ", %g", params->avgA[i][j]);   }
            fprintf(logFP, ")");
    }
    fprintf(logFP, " ] \n\n");

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.param%03d", out_name, s );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "mu, lambda");
        for( i = 0 ; i < cParams->sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g, %g", (double)(i+1)*params->avgMu, (double)(i+1)*params->avgLm);
            for( j = 0 ; j < cParams->sNo ; j++ )
            {   fprintf(fp, ", %g", params->avgA[j][i]);   }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    sprintf( fn, "%s.Lq%03d", out_name, s );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < results->iteration ; n++ ){
            fprintf( fp, "%24.20e\n", results->LqArr[n] );
        }
        fclose(fp);
    }

    sprintf( fn, "%s.maxS%03d", out_name, s );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < cParams->dLen ; n++ ){
            fprintf( fp, "%d\n", results->maxSumTraj[n] );
        }
        fclose(fp);
    }

}

//
