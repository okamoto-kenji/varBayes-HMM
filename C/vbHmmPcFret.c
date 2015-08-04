/*
 *  vbHmmPcFret.c
 *  Model-specific core functions for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include "vbHmmPcFret.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
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

void setFunctions_pcFret(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_pcFret;
    funcs.initialize_vbHmm = initialize_vbHmm_pcFret;
    funcs.freeParameters = freeParameters_pcFret;
    funcs.pTilde_z1 = pTilde_z1_pcFret;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_pcFret;
    funcs.pTilde_xn_zn = pTilde_xn_zn_pcFret;
    funcs.calcStatsVars = calcStatsVars_pcFret;
    funcs.maximization = maximization_pcFret;
    funcs.varLowerBound = varLowerBound_pcFret;
    funcs.reorderParameters = reorderParameters_pcFret;
    funcs.outputResults = outputResults_pcFret;
    setFunctions( funcs );
}


void **mallocParameterArray_pcFret( n )
size_t n;
{
    return (void**)malloc( n * sizeof(pcFretParameters*) );
}


void outputResults_pcFret( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputPcFretResults( cParams, (pcFretParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_pcFret( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    pcFretData *pc = xnWv->data;
    int sNo = cParams->sNo;
    pcFretParameters *params = blankParameters_pcFret( sNo );
    params->sNo = sNo;
    params->binSize = pc->binSize;

    int i, j;
    size_t totalC = 0;
    for( i = 0 ; i < xnWv->N ; i++ ){
        totalC += pc->dCounts[i] + pc->aCounts[i];
    }
    double meanI = (double)totalC / (double)xnWv->N;

    // hyper parameter for p( pi(i) )
    params->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->uPiArr[i] = 1.0;
        params->sumUPi += params->uPiArr[i];
    }

    // hyper parameter for p( A(i,j) )
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
    
    // hyper parameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        params->aIArr[i] = 1.0;
        params->bIArr[i] = 1.0 / meanI;
    }
    
    // hyper parameter for p( E(i) )
    for( i = 0 ; i < sNo ; i++ ){
        params->uEArr[i] = 1.0;
        params->vEArr[i] = 1.0;
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
    
    for( i = 0 ; i < sNo ; i++ ){
        params->avgE[i] = 0.5 + enoise(0.01);
        params->avgLnE[i][1] = log( params->avgE[i] );
        params->avgLnE[i][0] = log( 1.0 - params->avgE[i] );
    }

#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, "pi:%g, ", params->avgPi[i]);
        fprintf(logFP, "lnPi:%g, ", params->avgLnPi[i]);
        fprintf(logFP, "A(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", params->avgA[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "lnA(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", params->avgLnA[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "I:%g, ", params->avgI[i]);
        fprintf(logFP, "lnI:%g, ", params->avgLnI[i]);
        fprintf(logFP, "E:%g, ", params->avgE[i]);
        fprintf(logFP, "lnE:%g  \n", params->avgLnE[i][0]);
    }
    fprintf(logFP, "//\n");
}
#endif
    
    return params;
}

pcFretParameters *blankParameters_pcFret( sNo )
int sNo;
{
    int i;
    pcFretParameters *params = (pcFretParameters*)malloc( sizeof(pcFretParameters) );
    
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
    params->uEArr = (double*)malloc( sNo * sizeof(double) );
    params->vEArr = (double*)malloc( sNo * sizeof(double) );

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
    params->avgE = (double *)malloc( sNo * sizeof(double) );
    params->avgLnE = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->avgLnE[i] = (double *)malloc( 2 * sizeof(double) );   }

    params->Ni = (double *)malloc( sNo * sizeof(double) );
    params->Ci = (double *)malloc( sNo * sizeof(double) );
    params->Di = (double *)malloc( sNo * sizeof(double) );
    params->Ai = (double *)malloc( sNo * sizeof(double) );
    params->Mi = (double *)malloc( sNo * sizeof(double) );
    params->Nij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
    
    return params;
}


void freeParameters_pcFret( params )
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    int i;

    free( p->uPiArr );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->uAMat[i] );
    }
    free( p->uAMat );
    free( p->sumUAArr );
    free( p->aIArr );
    free( p->bIArr );
    free( p->uEArr );
    free( p->vEArr );
    
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
    free( p->avgE );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->avgLnE[i] );
    }
    free( p->avgLnE );

    free( p->Ni );
    free( p->Ci );
    free( p->Di );
    free( p->Ai );
    free( p->Mi );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Nij[i] );   }
    free( p->Nij );
    
    free( p );
    p = NULL;
}

void freePcFretDataSet( pcFretTraj )
xnDataSet *pcFretTraj;
{
    pcFretData *pc = pcFretTraj->data;
    free( pc->dCounts );
    free( pc->aCounts );
    free( pcFretTraj->data );
    free( pcFretTraj );
    pcFretTraj = NULL;
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

double pTilde_xn_zn_pcFret( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    pcFretData *xn = (pcFretData*)xnWv->data;
    return exp( xn->dCounts[n]*p->avgLnE[i][0] + xn->aCounts[n]*p->avgLnE[i][1] + (xn->dCounts[n]+xn->aCounts[n])*p->avgLnI[i] - p->avgI[i] ) / gsl_sf_fact( xn->dCounts[n] ) / gsl_sf_fact( xn->aCounts[n] );
}


void calcStatsVars_pcFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    pcFretData *xn = (pcFretData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Ci = p->Ci, *Di = p->Di, *Ai = p->Ai;
    double *Mi = p->Mi, **Nij = p->Nij;
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
            Di[i]  += gmMat[n][i] * (double)xn->dCounts[n];
            Ai[i]  += gmMat[n][i] * (double)xn->aCounts[n];
            for( j = 0 ; j < sNo ; j++ ){
                Mi[i]     += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
        Ci[i] = Di[i] + Ai[i];
    }

#ifdef DEBUG
#pragma omp critical
{
    for( n = 0 ; n < 20 ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(logFP, "%g,", gmMat[n][i]);
        }
        fprintf(logFP, "; ");
    }
    fprintf(logFP, "\n");
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, "Ni(%d)=%g,  ", i, Ni[i]);
        fprintf(logFP, "Ti(%d)=%g,  ", i, Ti[i]);
        fprintf(logFP, "Mi(%d)=%g,  ", i, Mi[i]);
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i )
                fprintf(logFP, "Nij(%d,%d)=%g, ", i, j, Nij[i][j]);
        }
        fprintf(logFP, "\n");
    }
}
#endif
}


void maximization_pcFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = p->Ni, *Ci = p->Ci, *Di = p->Di, *Ai = p->Ai;
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


double varLowerBound_pcFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, **avgLnE = p->avgLnE;
    double *Ni = p->Ni, *Ci = p->Ci, *Di = p->Di, *Ai = p->Ai;
    double *Mi = p->Mi, **Nij = p->Nij;
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

#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
    if( val > 100000 ){
        fprintf(logFP, "  > %g; %g; %g; %g;", lnpPi, lnpA, lnpI, lnpE);
        fprintf(logFP, " %g; %g; %g; %g; %g\n", lnqPi, lnqA, lnqI, lnqE, lnpX);
    }
}
#endif

    return val;
}    


void reorderParameters_pcFret( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    pcFretParameters *p = (pcFretParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA;
    double **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *Ni = p->Ni, *Ci = p->Ci, *Di = p->Di, *Ai = p->Ai;
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


void outputPcFretResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
pcFretParameters *params;
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

    fprintf(logFP, "   FRET efficiencies: ( %g", params->avgE[0]);
    for( i = 1 ; i < cParams->sNo ; i++ ){
        fprintf(logFP, ", %g", params->avgE[i]);
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
        fprintf(fp, "I, E");
        for( i = 0 ; i < cParams->sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g, %g", params->avgI[i], params->avgE[i]);
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
