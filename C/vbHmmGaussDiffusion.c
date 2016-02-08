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
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

void setFunctions_gaussDiff(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_gaussDiff;
    funcs.initialize_vbHmm = initialize_vbHmm_gaussDiff;
    funcs.freeParameters = freeParameters_gaussDiff;
    funcs.pTilde_z1 = pTilde_z1_gaussDiff;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_gaussDiff;
    funcs.pTilde_xn_zn = pTilde_xn_zn_gaussDiff;
    funcs.calcStatsVars = calcStatsVars_gaussDiff;
    funcs.maximization = maximization_gaussDiff;
    funcs.varLowerBound = varLowerBound_gaussDiff;
    funcs.reorderParameters = reorderParameters_gaussDiff;
    funcs.outputResults = outputResults_gaussDiff;
    setFunctions( funcs );
}


void **mallocParameterArray_gaussDiff( n )
size_t n;
{
    return (void**)malloc( n * sizeof(gaussDiffParameters*) );
}


void outputResults_gaussDiff( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputGaussDiffResults( cParams, (gaussDiffParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_gaussDiff( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    gaussDiffData *xn = xnWv->data;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    gaussDiffParameters *p = blankParameters_gaussDiff( sNo );
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
    calcStatsVars_gaussDiff( xnWv, cParams, p );
    maximization_gaussDiff( xnWv, cParams, p );

    return p;
}

gaussDiffParameters *blankParameters_gaussDiff( sNo )
int sNo;
{
    int i;
    gaussDiffParameters *p = (gaussDiffParameters*)malloc( sizeof(gaussDiffParameters) );
    
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
    p->avgDlt = (double *)malloc( sNo * sizeof(double) );
    p->avgLnDlt = (double *)malloc( sNo * sizeof(double) );

    p->uAArr = (double *)malloc( sNo * sizeof(double) );
    p->uBArr = (double *)malloc( sNo * sizeof(double) );
    p->aDlt = (double *)malloc( sNo * sizeof(double) );
    p->bDlt = (double *)malloc( sNo * sizeof(double) );

    p->Ni = (double *)malloc( sNo * sizeof(double) );
    p->Ri = (double *)malloc( sNo * sizeof(double) );
    p->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   p->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    p->Nii = (double *)malloc( sNo * sizeof(double) );
    
    return p;
}


void freeParameters_gaussDiff( params )
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
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
    free( p->avgDlt );
    free( p->avgLnDlt );

    free( p->uAArr );
    free( p->uBArr );

    free( p->Ni );
    free( p->Ri );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Nij[i] );   }
    free( p->Nij );
    free( p->Nii );
    
    free( p );
    p = NULL;
}

void freegaussDiffDataSet( gaussDiffTraj )
xnDataSet *gaussDiffTraj;
{
    gaussDiffData *xn = gaussDiffTraj->data;
    free( xn->v );
    free( gaussDiffTraj->data );
    free( gaussDiffTraj );
    gaussDiffTraj = NULL;
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


void calcStatsVars_gaussDiff( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    gaussDiffData *xn = (gaussDiffData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Ri = p->Ri, *Nii = p->Nii, **Nij = p->Nij;
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
            Ri[i] += gmMat[n][i] * pow( xn->v[n], 2.0 );
            for( j = 0 ; j < sNo ; j++ ){
                Nii[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
    }
}


void maximization_gaussDiff( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Ni = p->Ni, *Ri = p->Ri, *Nii = p->Nii, **Nij = p->Nij;
    double *aDlt = p->aDlt, *bDlt = p->bDlt;
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


double varLowerBound_gaussDiff( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uAArr = p->uAArr, *uBArr = p->uBArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Ni = p->Ni, *Ri = p->Ri, *Nii = p->Nii, **Nij = p->Nij;
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


void reorderParameters_gaussDiff( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    gaussDiffParameters *p = (gaussDiffParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgDlt = p->avgDlt, *avgLnDlt = p->avgLnDlt;
    double *Ni = p->Ni;
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


void outputGaussDiffResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
gaussDiffParameters *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    int i, j;
    fprintf(logFP, "  results: K = %d \n", s);

    fprintf(logFP, "   delta: ( %g", params->avgDlt[0]);
    for( i = 1 ; i < cParams->sNo ; i++ ){
        fprintf(logFP, ", %g", params->avgDlt[i]);
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
        fprintf(fp, "delta");
        for( i = 0 ; i < cParams->sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g", params->avgDlt[i]);
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
