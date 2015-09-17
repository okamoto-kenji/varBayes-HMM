/*
 *  vbHmmPc.c
 *  Model-specific core functions for VB-HMM-PC.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2015.09.17
 */

#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_pow_int.h>
#include "vbHmmPc.h"
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

void setFunctions_pc(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_pc;
    funcs.initialize_vbHmm = initialize_vbHmm_pc;
    funcs.freeParameters = freeParameters_pc;
    funcs.pTilde_z1 = pTilde_z1_pc;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_pc;
    funcs.pTilde_xn_zn = pTilde_xn_zn_pc;
    funcs.calcStatsVars = calcStatsVars_pc;
    funcs.maximization = maximization_pc;
    funcs.varLowerBound = varLowerBound_pc;
    funcs.reorderParameters = reorderParameters_pc;
    funcs.outputResults = outputResults_pc;
    setFunctions( funcs );
}


void **mallocParameterArray_pc( n )
size_t n;
{
    return (void**)malloc( n * sizeof(pcParameters*) );
}


void outputResults_pc( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputPcResults( cParams, (pcParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_pc( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    pcData *pc = xnWv->data;
    int sNo = cParams->sNo;
    pcParameters *params = blankParameters_pc( sNo );
    params->sNo = sNo;
    params->binSize = pc->binSize;

    int i, j;
    size_t totalC = 0;
    for( i = 0 ; i < xnWv->N ; i++ ){
        totalC += pc->counts[i];
    }
    double meanI = (double)totalC / (double)xnWv->N;

    // hyperparameter for p( pi(k) )
    params->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->uPiArr[i] = 1.0;
        params->sumUPi += params->uPiArr[i];
    }

    // hyperparameter for p( A(i,j) )
    for( i = 0 ; i < sNo ; i++ ){
        params->sumUAArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                params->uAMat[i][j] = 100.0;
            } else {
                params->uAMat[i][j] = 1.0;
            }
            params->sumUAArr[i] += params->uAMat[i][j];
        }
    }

    // hyperparameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        params->aIArr[i] = 1.0;
        params->bIArr[i] = 1.0 / meanI;
    }
    
    double sumPar = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] = 1.0/(double)sNo + enoise(0.1/(double)sNo);
        sumPar += params->avgPi[i];
    }
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] /= sumPar;
        params->avgLnPi[i] = log( params->avgPi[i] );
    }

    for( i = 0 ; i < sNo ; i++ ){
        sumPar = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                params->avgA[i][j] = 100.0 + enoise(1.0);
            } else {
                params->avgA[i][j] = 1.0 + enoise(0.01);
            }
            sumPar += params->avgA[i][j];
        }
        for( j = 0 ; j < sNo ; j++ ){
            params->avgA[i][j] /= sumPar;
            params->avgLnA[i][j] = log( params->avgA[i][j] );
        }
    }

    for( i = 0 ; i < sNo ; i++ ){
        params->avgI[i] = meanI + enoise(meanI/100.0);
        params->avgLnI[i] = log( params->avgI[i] );
    }

    return params;
}

pcParameters *blankParameters_pc( sNo )
int sNo;
{
    int i;
    pcParameters *params = (pcParameters*)malloc( sizeof(pcParameters) );

    params->sNo = 0;
    params->binSize = 0.0;
    params->uPiArr = (double*)malloc( sNo * sizeof(double) );
    params->sumUPi = 0.0;
    params->uAMat = (double**)malloc( sNo * sizeof(double*) );
    params->sumUAArr = (double*)malloc( sNo * sizeof(double) );
    for( i = 0 ; i < sNo ; i++ ){
        params->uAMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    params->aIArr = (double*)malloc( sNo * sizeof(double) );
    params->bIArr = (double*)malloc( sNo * sizeof(double) );
    params->avgPi = (double *)malloc( sNo * sizeof(double) );
    params->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    params->avgA = (double **)malloc( sNo * sizeof(double*) );
    params->avgLnA = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        params->avgA[i] = (double *)malloc( sNo * sizeof(double) );
        params->avgLnA[i] = (double *)malloc( sNo * sizeof(double) );
    }
    params->avgI = (double *)malloc( sNo * sizeof(double) );
    params->avgLnI = (double *)malloc( sNo * sizeof(double) );
    params->Ni = (double *)malloc( sNo * sizeof(double) );
    params->Ci = (double *)malloc( sNo * sizeof(double) );
    params->Mi = (double *)malloc( sNo * sizeof(double) );
    params->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    
    return params;
}


void freeParameters_pc( params )
void *params;
{
    pcParameters *p = (pcParameters*)params;
    int i;

    free( p->uPiArr );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->uAMat[i] );
    }
    free( p->uAMat );
    free( p->sumUAArr );
    free( p->aIArr );
    free( p->bIArr );
    
    free( p->avgPi );
    free( p->avgLnPi );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->avgA[i] );
        free( p->avgLnA[i] );
    }
    free( p->avgA );
    free( p->avgLnA );
    free( p->avgI );
    free( p->avgLnI );

    free( p->Ni );
    free( p->Ci );
    free( p->Mi );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Nij[i] );   }
    free( p->Nij );
    
    free( p );
    p = NULL;
}

void freePcDataSet( pcTraj )
xnDataSet *pcTraj;
{
    pcData *pc = pcTraj->data;
    free( pc->counts );
    free( pcTraj->data );
    free( pcTraj );
    pcTraj = NULL;
}


double pTilde_z1_pc( i, params )
int i;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_pc( i, j, params )
int i, j;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_pc( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    pcData *xn = (pcData*)xnWv->data;
    return exp( p->avgLnI[i] * (double)xn->counts[n] - p->avgI[i] ) / gsl_sf_fact( xn->counts[n] );
}


void calcStatsVars_pc( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    pcData *xn = (pcData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Ci = p->Ci;
    double *Mi = p->Mi, **Nij = p->Nij;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i] = 1e-10;
        Ci[i] = 1e-10;
        Mi[i] = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }
    }
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            Ni[i]  += gmMat[n][i];
            Ci[i]  += gmMat[n][i] * (double)xn->counts[n];
            for( j = 0 ; j < sNo ; j++ ){
                Mi[i]    += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
    }
}

void maximization_pc( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *Ni = p->Ni, *Ci = p->Ci;
    double *Mi = p->Mi, **Nij = p->Nij;
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
        pcData *pc = xnWv->data;
        for( n = 0 ; n < xnWv->N ; n++ ){
            totalC += pc->counts[n];
        }
        double meanI = (double)totalC / (double)xnWv->N;
        avgI[i] = MIN( avgI[i], maxIntensityRatio * meanI );
        avgLnI[i] = MIN( avgLnI[i], log(maxIntensityRatio * meanI) );
#endif
    }
}

double varLowerBound_pc( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *Ni = p->Ni, *Ci = p->Ci;
    double *Mi = p->Mi, **Nij = p->Nij;
    size_t n;
    int i, j;
    
    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpI = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqI = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpI += aIArr[i] * log(bIArr[i]) - gsl_sf_lngamma(aIArr[i]);
        lnpI += (aIArr[i] - 1.0) * avgLnI[i] - bIArr[i] * avgI[i];

        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnqI += (Ci[i] + aIArr[i]) * log(Ni[i] + bIArr[i]) - gsl_sf_lngamma(Ci[i] + aIArr[i]);
        lnqI += (Ci[i] + aIArr[i] - 1.0) * avgLnI[i] - (Ni[i] + bIArr[i]) * avgI[i];

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Mi[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0) * avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi( sumUAArr[i]+Mi[i]) );
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }
        
    double val;
    val  = lnpPi + lnpA + lnpI;
    val -= lnqPi + lnqA + lnqI;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));

    return val;
}    


void reorderParameters_pc( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA;
    double **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = p->Ni, *Ci = p->Ci;
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

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ci[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ci[i] = store[i];   }

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


void outputPcResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
pcParameters *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    int i, j;
    fprintf(logFP, "  results: K = %d \n", s);
    fprintf(logFP, "   intensities: ( %g", params->avgI[0]);
    for( i = 1 ; i < cParams->sNo ; i++ )
    {   fprintf(logFP, ", %g", params->avgI[i]);   }
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
        fprintf(fp, "I");
        for( i = 0 ; i < cParams->sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g", params->avgI[i]);
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
