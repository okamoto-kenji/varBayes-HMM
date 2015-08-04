/*
 *  vbHmmTs.c
 *  Model-specific core functions for VB-HMM-TS.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include "vbHmmTs.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

//#define  DEBUG

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

void setFunctions_ts(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_ts;
    funcs.initialize_vbHmm = initialize_vbHmm_ts;
    funcs.freeParameters = freeParameters_ts;
    funcs.pTilde_z1 = pTilde_z1_ts;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_ts;
    funcs.pTilde_xn_zn = pTilde_xn_zn_ts;
    funcs.calcStatsVars = calcStatsVars_ts;
    funcs.maximization = maximization_ts;
    funcs.varLowerBound = varLowerBound_ts;
    funcs.reorderParameters = reorderParameters_ts;
    funcs.outputResults = outputResults_ts;
    setFunctions( funcs );
}

void **mallocParameterArray_ts( n )
size_t n;
{
    return (void**)malloc( n * sizeof(tsParameters*) );
}


void outputResults_ts( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputTsResults( cParams, (tsParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_ts( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    int sNo = cParams->sNo;
    tsParameters *params = blankParameters( sNo );
    params->sNo = sNo;

    int i, j;

    // hyperparameter for p( pi(k) )
    params->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->uPiArr[i] = 1.0;
        params->sumUPi += params->uPiArr[i];
    }

    // hyperparameter for p( k(i,j) ) (i != j)
    for( i = 0 ; i < sNo ; i++ ){
        params->sumUKArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            params->uKMat[i][j] = 1.0;
            if( j != i ){
                params->sumUKArr[i] += params->uKMat[i][j];
            }
        }
    }

    double meanI = (double)xnWv->N / xnWv->T;
    double sumPar = 0.0;

    // hyperparameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        params->aIArr[i] = 1.0;
        params->bIArr[i] = 1.0 / meanI;
    }
    

    sumPar = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] = 1.0/(double)sNo + enoise(0.1/(double)sNo);
        sumPar += params->avgPi[i];
    }
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] /= sumPar;
        params->avgLnPi[i] = log( params->avgPi[i] );
    }

    for( i = 0 ; i < sNo ; i++ ){
        params->avgI[i] = meanI + 0.1 * enoise(meanI);
        params->avgLnI[i] = log( params->avgI[i] );
    }

    for( i = 0 ; i < sNo ; i++ ){
        params->avgK[i][i] = meanI/100.0 + enoise(meanI/1000.0);
        params->avgLnK[i][i] = log( params->avgK[i][i] );
        params->avgLnKI[i] = log( params->avgK[i][i] + params->avgI[i] );
        sumPar = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                params->avgK[i][j] = 1.0 + enoise(0.1);
                sumPar += params->avgK[i][j];
            }
        }
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                params->avgK[i][j] /= sumPar;
                params->avgLnK[i][j] = log( params->avgK[i][j] );
            }
        }
    }

#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, "pi:%g, ", params->avgPi[i]);
        fprintf(logFP, "lnPi:%g, ", params->avgLnPi[i]);
        fprintf(logFP, "k(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", params->avgK[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "lnK(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", params->avgLnK[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "lnKI:%g, ", params->avgLnKI[i]);
        fprintf(logFP, "I:%g, ", params->avgI[i]);
        fprintf(logFP, "lnI:%g  \n", params->avgLnI[i]);
    }
    fprintf(logFP, "//\n");
}
#endif

    return params;
}


tsParameters *blankParameters( sNo )
int sNo;
{
    int i;
    tsParameters *params = (tsParameters*)malloc( sizeof(tsParameters) );
    params->uPiArr = (double*)malloc( sNo * sizeof(double) );
    params->sumUPi = 0.0;
    // hyperparameter for p( k(i,j) ) (i != j)
    params->uKMat = (double**)malloc( sNo * sizeof(double*) );
    params->sumUKArr = (double*)malloc( sNo * sizeof(double) );
    for( i = 0 ; i < sNo ; i++ ){
        params->uKMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    // hyperparameter for p( I(k) )
    params->aIArr = (double*)malloc( sNo * sizeof(double) );
    params->bIArr = (double*)malloc( sNo * sizeof(double) );

    params->avgPi = (double *)malloc( sNo * sizeof(double) );
    params->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    params->avgI = (double *)malloc( sNo * sizeof(double) );
    params->avgLnI = (double *)malloc( sNo * sizeof(double) );
    params->avgK = (double **)malloc( sNo * sizeof(double*) );
    params->avgLnK = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        params->avgK[i] = (double *)malloc( sNo * sizeof(double) );
        params->avgLnK[i] = (double *)malloc( sNo * sizeof(double) );
    }
    params->avgLnKI = (double *)malloc( sNo * sizeof(double) );
    
    params->Ni = (double *)malloc( sNo * sizeof(double) );
    params->Ti = (double *)malloc( sNo * sizeof(double) );
    params->Nii = (double *)malloc( sNo * sizeof(double) );
    params->Nij = (double *)malloc( sNo * sizeof(double) );
    params->Mij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->Mij[i] = (double *)malloc( sNo * sizeof(double) );   }

    return params;
}

void freeParameters_ts( params )
//vbHmmCommonParameters *cParams;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    int i;

    free( p->uPiArr );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->uKMat[i] );
    }
    free( p->uKMat );
    free( p->sumUKArr );
    free( p->aIArr );
    free( p->bIArr );

    free( p->avgPi );
    free( p->avgLnPi );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->avgK[i] );
        free( p->avgLnK[i] );
    }
    free( p->avgK );
    free( p->avgLnK );
    free( p->avgLnKI );
    free( p->avgI );
    free( p->avgLnI );
    free( p->Ni );
    free( p->Ti );
    free( p->Nii );
    free( p->Nij );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Mij[i] );   }
    free( p->Mij );
    
    free( p );
    p = NULL;
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
    tsData *xn = (tsData*)xnWv->data;
    return exp( p->avgLnI[i] - p->avgI[i] * xn[n].dt );
}


void calcStatsVars_ts( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    tsData *xn = (tsData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Ti = p->Ti;
    double *Nii = p->Nii, *Nij = p->Nij, **Mij = p->Mij;
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
            Ti[i]  += gmMat[n][i] * (xn[n]).dt;
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

#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
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
        fprintf(logFP, "Nii(%d)=%g,  ", i, Nii[i]);
        fprintf(logFP, "Nij(%d)=%g,  ", i, Nij[i]);
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i )
                fprintf(logFP, "N'ij(%d,%d)=%g, ", i, j, Mij[i][j]);
        }
        fprintf(logFP, "\n");
    }
}
#endif
}


void maximization_ts( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *Ni = p->Ni, *Ti = p->Ti;
    double *Nii = p->Nii, *Nij = p->Nij, **Mij = p->Mij;
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

#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, "pi:%g, ", avgPi[i]);
        fprintf(logFP, "lnPi:%g, ", avgLnPi[i]);
        fprintf(logFP, "k(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", avgK[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "lnK(");
        for( j = 0 ; j < sNo ; j++ ){
            fprintf(logFP, "%g,", avgLnK[i][j]);
        }
        fprintf(logFP, "), ");
        fprintf(logFP, "lnKI:%g, ", avgLnKI[i]);
        fprintf(logFP, "I:%g, ", avgI[i]);
        fprintf(logFP, "lnI:%g  \n", avgLnI[i]);
    }
    fprintf(logFP, "//\n");
}
#endif
}


double varLowerBound_ts( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *avgLnPi = p->avgLnPi, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *Ni = p->Ni, *Ti = p->Ti;
    double *Nii = p->Nii, *Nij = p->Nij, **Mij = p->Mij;
    size_t n;
    int i, j;
    
    double lnpPi;
    lnpPi = gsl_sf_lngamma(sumUPi);
    double Ck = 1.0, lnpKii = sNo * log(Ck);
    double lnpKij = 0.0;
    double lnpI = 0.0;
    double lnqPi;
    lnqPi = gsl_sf_lngamma(sumUPi + 1);
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

#ifdef DEBUG
#pragma omp critical
{        
    FILE *logFP = stderr;
    if( val > 100000 ){
        fprintf(logFP, "  > %g; %g; %g; %g;", lnpPi, lnpKii, lnpKij, lnpI);
        fprintf(logFP, " %g; %g; %g; %g\n", lnqPi, lnqKiiI, lnqKij, lnpX);
    }
}
#endif
            
    return val;
}    


void reorderParameters_ts( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    tsParameters *p = (tsParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK;
    double **avgLnK = p->avgLnK, *avgLnKI = p->avgLnKI;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = p->Ni, *Ti = p->Ti;
    size_t n;
    int i, j;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( sNo * sizeof(double) );   }

    // index indicates order of avgI values (0=biggest avgI -- sNo=smallest avgI).
#ifdef DEBUG
#pragma omp critical
{
    FILE *logFP = stderr;
    fprintf(logFP, "I (");
    for( i = 0 ; i < sNo ; i++ )
    {   fprintf(logFP, "%g,", avgI[i]);   }
    fprintf(logFP, ") -> (");
}
#endif
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
#ifdef DEBUG
#pragma omp critical
{
    for( i = 0 ; i < sNo ; i++ )
    {   fprintf(logFP, "%d,", index[i]);   }
    fprintf(logFP, ")\n");
}
#endif

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


void outputTsResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
tsParameters *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{
    int i, j;
    if( logFP != NULL ){
        fprintf(logFP, "  results: K = %d \n", s);
        fprintf(logFP, "   intensities: ( %g", params->avgI[0]);
        for( i = 1 ; i < cParams->sNo ; i++ )
        {   fprintf(logFP, ", %g", params->avgI[i]);   }
        fprintf(logFP, " ) \n");
        fprintf(logFP, "   k_matrix: [");
        for( i = 0 ; i < cParams->sNo ; i++ ){
                fprintf(logFP, " ( %g", params->avgK[i][0]);
                for( j = 1 ; j < cParams->sNo ; j++ )
                {   fprintf(logFP, ", %g", params->avgK[i][j]);   }
                fprintf(logFP, ")");
        }
        fprintf(logFP, " ] \n\n");
    }

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.param%03d", out_name, s );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "I");
        for( i = 0 ; i < cParams->sNo ; i++ )
        {   fprintf(fp, ", K%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g", params->avgI[i]);
            for( j = 0 ; j < cParams->sNo ; j++ )
            {   fprintf(fp, ", %g", params->avgK[j][i]);   }
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
