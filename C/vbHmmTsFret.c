/*
 *  vbHmmTsFret.c
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

#include "vbHmmTsFret.h"
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


void setFunctions_tsFret(){
    commonFunctions funcs;
    funcs.mallocParameterArray = mallocParameterArray_tsFret;
    funcs.initialize_vbHmm = initialize_vbHmm_tsFret;
    funcs.freeParameters = freeParameters_tsFret;
    funcs.pTilde_z1 = pTilde_z1_tsFret;
    funcs.pTilde_zn_zn1 = pTilde_zn_zn1_tsFret;
    funcs.pTilde_xn_zn = pTilde_xn_zn_tsFret;
    funcs.calcStatsVars = calcStatsVars_tsFret;
    funcs.maximization = maximization_tsFret;
    funcs.varLowerBound = varLowerBound_tsFret;
    funcs.reorderParameters = reorderParameters_tsFret;
    funcs.outputResults = outputResults_tsFret;
    setFunctions( funcs );
}


void **mallocParameterArray_tsFret( n )
size_t n;
{
    return (void**)malloc( n * sizeof(tsFretParameters*) );
}


void outputResults_tsFret( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int s;
char *out_name;
FILE *logFP;
{    
    outputTsFretResults( cParams, (tsFretParameters*)params, results, s, out_name, logFP );
}


void *initialize_vbHmm_tsFret( xnWv, cParams )
xnDataSet *xnWv;
vbHmmCommonParameters* cParams;
{
    int sNo = cParams->sNo;
    tsFretParameters *params = blankParameters_tsFret( sNo );
    params->sNo = cParams->sNo;
    int i, j;

    // hyperparameter for p( pi(i) )
    params->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        params->uPiArr[i] = 1.0;
        params->sumUPi += params->uPiArr[i];
    }

    // hyperparameter for p( k(i,j) ) (i != j)
    for( i = 0 ; i < sNo ; i++ ){
        params->uKMat[i] = (double*)malloc( sNo * sizeof(double) );
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

    // hyperparameter for p( E(i) )
    for( i = 0 ; i < sNo ; i++ ){
        params->uEArr[i] = 1.0;
        params->vEArr[i] = 1.0;
    }


    // initial conditions
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] = 1.0/(double)sNo + enoise(0.01/(double)sNo);
        sumPar += params->avgPi[i];
    }
    for( i = 0 ; i < sNo ; i++ ){
        params->avgPi[i] /= sumPar;
        params->avgLnPi[i] = log( params->avgPi[i] );
    }

    for( i = 0 ; i < sNo ; i++ ){
        params->avgI[i] = meanI + enoise(meanI/100.0);
        params->avgLnI[i] = log( params->avgI[i] );
    }

    for( i = 0 ; i < sNo ; i++ ){
        params->avgK[i][i] = meanI/100.0 + enoise(meanI/10000.0);
        params->avgLnK[i][i] = log( params->avgK[i][i] );
        params->avgLnKI[i] = log( params->avgK[i][i] + params->avgI[i] );
        sumPar = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                params->avgK[i][j] = 1.0 + enoise(0.01);
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

    for( i = 0 ; i < sNo ; i++ ){
        params->avgE[i] = 0.5 + enoise(0.01);
        params->avgLnE[i][1] = log( params->avgE[i] );
        params->avgLnE[i][0] = log( 1.0 - params->avgE[i] );
    }

    return params;
}


tsFretParameters *blankParameters_tsFret( sNo )
int sNo;
{
    int i;
    tsFretParameters *params = (tsFretParameters*)malloc( sizeof(tsFretParameters) );
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
    // hyperparameter for p( E(i) )
    params->uEArr = (double*)malloc( sNo * sizeof(double) );
    params->vEArr = (double*)malloc( sNo * sizeof(double) );

    // parameters
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
    params->avgE = (double *)malloc( sNo * sizeof(double) );
    params->avgLnE = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->avgLnE[i] = (double *)malloc( 2 * sizeof(double) );   }

    params->Ni = (double *)malloc( sNo * sizeof(double) );
    params->Ti = (double *)malloc( sNo * sizeof(double) );
    params->eps = (double *)malloc( sNo * sizeof(double) );
    params->Nii = (double *)malloc( sNo * sizeof(double) );
    params->Nij = (double *)malloc( sNo * sizeof(double) );
    params->Mij = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   params->Mij[i] = (double *)malloc( sNo * sizeof(double) );   }

    return params;
}


void freeParameters_tsFret( params )
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    int i;

    free( p->uPiArr );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->uKMat[i] );
    }
    free( p->uKMat );
    free( p->sumUKArr );
    free( p->aIArr );
    free( p->bIArr );
    free( p->uEArr );
    free( p->vEArr );
    
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
    free( p->avgE );
    for( i = 0 ; i < p->sNo ; i++ ){
        free( p->avgLnE[i] );
    }
    free( p->avgLnE );
    free( p->Ni );
    free( p->Ti );
    free( p->eps );
    free( p->Nii );
    free( p->Nij );
    for( i = 0 ; i < p->sNo ; i++ )
    {   free( p->Mij[i] );   }
    free( p->Mij );
    
    free( p );
    p = NULL;
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

double pTilde_xn_zn_tsFret( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    tsFretData *xn = (tsFretData*)xnWv->data;
    return exp( p->avgLnI[i] - (p->avgI[i] * xn[n].dt) + p->avgLnE[i][xn[n].ch] );
}



void calcStatsVars_tsFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    tsFretData *xn = (tsFretData*)xnWv->data;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *Ni = p->Ni, *Ti = p->Ti, *eps = p->eps;
    double *Nii = p->Nii, *Nij = p->Nij, **Mij = p->Mij;
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
            Ti[i]  += gmMat[n][i] * xn[n].dt;
            eps[i] += gmMat[n][i] * xn[n].ch;
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
//        fprintf(logFP, "eps(%d)=%g,  ", i, eps[i]);
//        fprintf(logFP, "Nii(%d)=%g,  ", i, Nii[i]);
//        fprintf(logFP, "Nij(%d)=%g,  ", i, Nij[i]);
//        for( j = 0 ; j < sNo ; j++ ){
//            if( j != i )
//                fprintf(logFP, "N'ij(%d,%d)=%g, ", i, j, Mij[i][j]);
//        }
//        fprintf(logFP, "\n");
//    }
//}
//#endif
}


void maximization_tsFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK, **avgLnK = p->avgLnK;
    double *avgLnKI = p->avgLnKI, *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = p->Ni, *Ti = p->Ti, *eps = p->eps;
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

        avgE[i] = ( uEArr[i] + eps[i] ) / ( uEArr[i] + vEArr[i] + Ni[i] );
        avgLnE[i][0]  = gsl_sf_psi( vEArr[i] + Ni[i] - eps[i] );    // ln(1-E) for donor
        avgLnE[i][0] -= gsl_sf_psi( uEArr[i] + vEArr[i] + Ni[i] );
        avgLnE[i][1]  = gsl_sf_psi( uEArr[i] + eps[i] );            // ln(E) for acceptor
        avgLnE[i][1] -= gsl_sf_psi( uEArr[i] + vEArr[i] + Ni[i] );
    }

}


double varLowerBound_tsFret( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, *cn = cParams->cn;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uKMat = p->uKMat, *sumUKArr = p->sumUKArr;
    double *uEArr = p->uEArr, *vEArr = p->vEArr;
    double *avgLnPi = p->avgLnPi, **avgLnK = p->avgLnK, *avgLnKI = p->avgLnKI;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, **avgLnE = p->avgLnE;
    double *Ni = p->Ni, *Ti = p->Ti, *eps = p->eps;
    double *Nii = p->Nii, *Nij = p->Nij, **Mij = p->Mij;
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

//#ifdef DEBUG
//#pragma omp critical
//{
//    FILE *logFP = stderr;
//    if( val > 100000 ){
//        fprintf(logFP, "  > %g; %g; %g; %g; %g;", lnpPi, lnpKii, lnpKij, lnpI, lnpE);
//        fprintf(logFP, " %g; %g; %g; %g; %g\n", lnqPi, lnqKiiI, lnqKij, lnqE, lnpX);
//    }
//}
//#endif
    
    return val;
}    


void reorderParameters_tsFret( cParams, params )
vbHmmCommonParameters *cParams;
void *params;
{
    tsFretParameters *p = (tsFretParameters*)params;
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat, ***xiMat = cParams->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgK = p->avgK;
    double **avgLnK = p->avgLnK, *avgLnKI = p->avgLnKI;
    double *avgI = p->avgI, *avgLnI = p->avgLnI, *avgE = p->avgE, **avgLnE = p->avgLnE;
    double *Ni = p->Ni, *Ti = p->Ti, *eps = p->eps;
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


void outputTsFretResults( cParams, params, results, s, out_name, logFP )
vbHmmCommonParameters *cParams;
tsFretParameters *params;
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

    fprintf(logFP, "   k_matrix: [");
    for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(logFP, " ( %g", params->avgK[i][0]);
            for( j = 1 ; j < cParams->sNo ; j++ )
            {   fprintf(logFP, ", %g", params->avgK[i][j]);   }
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
        {   fprintf(fp, ", K%dx", i);   }
        fprintf(fp, "\n");

        for( i = 0 ; i < cParams->sNo ; i++ ){
            fprintf(fp, "%g, %g", params->avgI[i], params->avgE[i]);
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
