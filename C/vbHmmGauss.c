/*
 *  vbHmmGauss.c
 *  Model-specific core functions for VB-HMM-GAUSS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmGauss.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

void setFunctions_gauss(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_gauss;
    funcs.initialize_vbHmm = initialize_vbHmm_gauss;
    funcs.freeParameters = freeParameters_gauss;
    funcs.pTilde_z1 = pTilde_z1_gauss;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_gauss;
    funcs.pTilde_xn_zn = pTilde_xn_zn_gauss;
    funcs.calcStatsVars = calcStatsVars_gauss;
    funcs.maximization = maximization_gauss;
    funcs.varLowerBound = varLowerBound_gauss;
    funcs.reorderParameters = reorderParameters_gauss;
    funcs.outputResults = outputResults_gauss;
    setFunctions( funcs );
}


void **mallocParameterArray_gauss( n )
size_t n;
{
    return (void**)malloc( n * sizeof(gaussParameters*) );
}


void outputResults_gauss( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputGaussResults( cParams, (gaussParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_gauss( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    gaussData *xn = xnWv->data;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    gaussParameters *p = blankParameters_gauss( sNo );
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
        p->uBtArr[i] = precX / 250.0;
        p->uMuArr[i] = meanX;
        p->uAArr[i] = 2.5;
        p->uBArr[i] = 0.01;
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
    calcStatsVars_gauss( xnWv, cParams, p );
    maximization_gauss( xnWv, cParams, p );

    return p;
}

gaussParameters *blankParameters_gauss( sNo )
int sNo;
{
    int i;
    gaussParameters *p = (gaussParameters*)malloc( sizeof(gaussParameters) );
    
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

    p->Ni = (double *)malloc( sNo * sizeof(double) );
    p->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   p->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    p->Nii = (double *)malloc( sNo * sizeof(double) );
    p->barX = (double *)malloc( sNo * sizeof(double) );
    p->NiSi = (double *)malloc( sNo * sizeof(double) );
    
    return p;
}


void freeParameters_gauss( params )
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
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
    free( p->avgMu );
    free( p->avgLm );
    free( p->avgLnLm );

    free( p->uBtArr );
    free( p->uMuArr );
    free( p->uAArr );
    free( p->uBArr );
    free( p->btMu );
    free( p->aLm );
    free( p->bLm );
    free( p->mu0 );

    free( p->Ni );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Nij[i] );   }
    free( p->Nij );
    free( p->Nii );
    free( p->barX );
    free( p->NiSi );
    
    free( p );
    p = NULL;
}

void freeGaussDataSet( gaussTraj )
xnDataSet *gaussTraj;
{
    gaussData *xn = gaussTraj->data;
    free( xn->v );
    free( gaussTraj->data );
    free( gaussTraj );
    gaussTraj = NULL;
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

double pTilde_xn_zn_gauss( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    gaussData *xn = (gaussData*)xnWv->data;
    double val;
    val  = p->avgLnLm[i] - log(2.0 * M_PI);
    val -= 1.0/p->btMu[i] + p->aLm[i] / p->bLm[i] * pow( xn->v[n] - p->mu0[i], 2.0);
    return exp(val / 2.0);
}


void calcStatsVars_gauss( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    gaussData *xn = (gaussData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Nii = p->Nii, **Nij = p->Nij, *Ni = p->Ni;
    double *barX = p->barX, *NiSi = p->NiSi;
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
            barX[i] += gmMat[n][i] * xn->v[n];
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        barX[i] /= Ni[i];
        for( n = 0 ; n < dLen ; n++ ){
            NiSi[i] += gmMat[n][i] * pow( xn->v[n] - barX[i], 2.0);
        }
    }
}


void maximization_gauss( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *Ni = p->Ni, *Nii = p->Nii, **Nij = p->Nij, *barX = p->barX, *NiSi = p->NiSi;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
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
        aLm[i]  = uAArr[i] + (Ni[i] + 1.0) / 2.0;
        bLm[i]  = uBArr[i] + (NiSi[i] / 2.0);
        bLm[i] += uBtArr[i] * Ni[i] * pow( barX[i] - uMuArr[i], 2.0) / 2.0 / (uBtArr[i] + Ni[i]);
        
        avgMu[i]    = mu0[i];
        avgLm[i]    = aLm[i] / bLm[i];
        avgLnLm[i]  = gsl_sf_psi( aLm[i] ) - log( bLm[i] );
    }
}


double varLowerBound_gauss( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uBtArr = p->uBtArr, *uMuArr = p->uMuArr, *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *Nii = p->Nii, **Nij = p->Nij;
    double *mu0 = p->mu0, *btMu = p->btMu, *aLm = p->aLm, *bLm = p->bLm;
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


void reorderParameters_gauss( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    gaussParameters *p = (gaussParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgMu = p->avgMu, *avgLm = p->avgLm, *avgLnLm = p->avgLnLm;
    double *Ni = p->Ni;
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


void outputGaussResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
gaussParameters *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    int i, j;
    fprintf(logFP, "  results: K = %d \n", s);

    fprintf(logFP, "   means: ( %g", params->avgMu[0]);
    for( i = 1 ; i < cParams->sNo ; i++ )
    {   fprintf(logFP, ", %g", params->avgMu[i]);   }
    fprintf(logFP, " ) \n");

    fprintf(logFP, "   lambda: ( %g", params->avgLm[0]);
    for( i = 1 ; i < cParams->sNo ; i++ ){
        fprintf(logFP, ", %g", params->avgLm[i]);
    }
    fprintf(logFP, " ) \n");

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
            fprintf(fp, "%g, %g", params->avgMu[i], params->avgLm[i]);
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
