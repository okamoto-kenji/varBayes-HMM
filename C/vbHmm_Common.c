/*
 *  vbHmm_Common.c
 *  Common VB-HMM engine.
 *  Reference: 
 *    Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include "vbHmm_Common.h"
#include <string.h>
#include <float.h>


// Pointers of functions to call model-specific functions.
static mallocParameterArray_func  mallocParameterArray  = NULL;
static initialize_vbHmm_func      initialize_vbHmm      = NULL;
static freeParameters_func        freeParameters        = NULL;
static pTilde_z1_func             pTilde_z1             = NULL;
static pTilde_zn_zn1_func         pTilde_zn_zn1         = NULL;
static pTilde_xn_zn_func          pTilde_xn_zn          = NULL;
static calcStatsVars_func         calcStatsVars         = NULL;
static maximization_func          maximization          = NULL;
static varLowerBound_func         varLowerBound         = NULL;
static reorderParameters_func     reorderParameters     = NULL;
static outputResults_func         outputResults         = NULL;


// This function must be called to connect with the model before executing analysis.
void setFunctions( funcs )
commonFunctions funcs;
{
    mallocParameterArray = funcs.mallocParameterArray;
    initialize_vbHmm = funcs.initialize_vbHmm;
    freeParameters = funcs.freeParameters;
    pTilde_z1 = funcs.pTilde_z1;
    pTilde_zn_zn1 = funcs.pTilde_zn_zn1;
    pTilde_xn_zn = funcs.pTilde_xn_zn;
    calcStatsVars = funcs.calcStatsVars;
    maximization = funcs.maximization;
    varLowerBound = funcs.varLowerBound;
    reorderParameters = funcs.reorderParameters;
    outputResults = funcs.outputResults;
}


//////////////////////////////////////////////////////////////////  VB-HMM Execution Functions
int modelComparison( xnWv, sFrom, sTo, trials, maxIteration, threshold, out_name, logFP )
xnDataSet *xnWv;
int sFrom, sTo, trials;
int maxIteration;
double threshold;
char *out_name;
FILE *logFP;
{
    int s, t;
    if( logFP != NULL ){
        fprintf( logFP, "  No. of states from %d to %d,  trials = %d, ", sFrom, sTo, trials);
        fprintf( logFP, "  analyze:  maxIteration = %d,  threshold = %g \n\n", maxIteration, threshold);
    }

    double *LqVsK = malloc( trials * (sTo - sFrom + 1) * sizeof(double) );

    int maxS = 0;
    double maxLq = -DBL_MAX;
    vbHmmResults **resArray = (vbHmmResults**)malloc( trials * sizeof(vbHmmResults*) );
    vbHmmCommonParameters **cParArray = (vbHmmCommonParameters**)malloc( trials * sizeof(vbHmmCommonParameters*) );
    void **parArray = (*mallocParameterArray)( trials );
    for( s = sFrom ; s <= sTo ; s++ ){
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
        for( t = 0 ; t < trials ; t++ ){
            int st = (s - sFrom) * trials + t;
            resArray[t] = newResults();
            cParArray[t] = initialize_Common( xnWv, s );
            parArray[t] = (*initialize_vbHmm)( xnWv, cParArray[t] );
            
            LqVsK[st] = vbHmm_Main( xnWv, cParArray[t], parArray[t], resArray[t], maxIteration, threshold, logFP );
            if( LqVsK[st] > maxLq ){
                maxLq = LqVsK[st];
                maxS = s;
            }
        }
        
        double maxLqForS = 0.0;
        int maxT = 0;
        for( t = 0 ; t < trials ; t++ ){
            int st = (s - sFrom) * trials + t;
            if( LqVsK[st] > maxLqForS ){
                maxLqForS = LqVsK[st];
                maxT = t;
            }
        }
        (*outputResults)( cParArray[maxT], parArray[maxT], resArray[maxT], s, out_name, logFP );
        
        for( t = 0 ; t < trials ; t++ ){
            (*freeParameters)( parArray[t] );
            freeCommonParameters( cParArray[t] );
            freeResults( resArray[t] );
        }
        
        if( s >= (maxS+3) ){
            s++;
            break;
        }
        
    }
    sTo = s - 1;
    free( resArray );
    free( cParArray );
    free( parArray );
    
    char fn[256];
    FILE *fp = NULL;
    strncpy( fn, out_name, sizeof(fn) );
    strncat( fn, ".LqVsK", sizeof(fn) - strlen(fn) - 1 );
    if( (fp = fopen( fn, "w" )) != NULL ){
        for( s = 0 ; s < trials * (sTo - sFrom + 1) ; s++ ){
            fprintf( fp, "%2d, %.20g\n", (s/trials) + sFrom, LqVsK[s] );
        }
        fclose(fp);
    }
    free( LqVsK );
    return maxS;
}


//////////////////////////////////////////////////////////////////  VB-HMM Common Engine
double vbHmm_Main( xnWv, cParams, params, results, maxIteration, threshold, logFP )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
vbHmmResults *results;
int maxIteration;
double threshold;
FILE *logFP;
{
    double **LqArr = &results->LqArr;
    *LqArr = realloc( *LqArr, maxIteration * sizeof(double) );

    int i;
    for( i = 0 ; i < maxIteration ; i++ ){
        // E-step
        forwardBackward( xnWv, cParams, params );

        (*calcStatsVars)( xnWv, cParams, params );
        (*LqArr)[i] = (*varLowerBound)( xnWv, cParams, params );

        // End loop if derivative of variational lower bound reaches threshold.
        if( (i>0) && ( fabs( ((*LqArr)[i] - (*LqArr)[i-1]) / (*LqArr)[i] ) < threshold ) ){
            break;
        }

        // M-step
        (*maximization)( xnWv, cParams, params );
    }
    if( i == maxIteration ){
        if( logFP != NULL ){
            fprintf(logFP, "MAX iteration (%d) reached.\n", maxIteration);
        }
        i--;
    }
    (*reorderParameters)( cParams, params );
    maxSum( xnWv, cParams, params, &(results->maxSumTraj) );

    results->iteration = i+1;
    *LqArr = realloc( *LqArr, (i+1) * sizeof(double) );
    results->maxLq = (*LqArr)[i];
    if( logFP != NULL ){
        fprintf( logFP, "  iteration: %d    evidence p(x|K=%d) = %.20g \n", i+1, cParams->sNo, results->maxLq );
    }
    return results->maxLq;
}


// Initializes parameters commonly used in VB-HMM analysis
vbHmmCommonParameters *initialize_Common( xnWv, sNo )
xnDataSet *xnWv;
int sNo;
{
    vbHmmCommonParameters *cParams = (vbHmmCommonParameters*)malloc( sizeof(vbHmmCommonParameters) );
    cParams->dLen = xnWv->N;
    cParams->sNo = sNo;
    size_t dLen = cParams->dLen;
    int i, n;

    // gamma
    cParams->gmMat = (double**)malloc( dLen * sizeof(double*) );
    double **gmMat = cParams->gmMat;
    for( n = 0 ; n < dLen ; n++ )
    {
        gmMat[n] = (double*)malloc( sNo * sizeof(double) );
        memset( gmMat[n], 0, sNo * sizeof(double) );
    }
    // xi
    cParams->xiMat = (double***)malloc( dLen * sizeof(double**) );
    double ***xiMat = cParams->xiMat;
    for( n = 0 ; n < dLen ; n++ ){
        xiMat[n] = (double**)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ ){
            xiMat[n][i] = (double*)malloc( sNo * sizeof(double) );
            memset( xiMat[n][i], 0, sNo * sizeof(double) );
        }
    }

    // alpha for E-step
    cParams->aMat = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   cParams->aMat[n] = (double*)malloc( sNo * sizeof(double) );  }
    // beta for E-step
    cParams->bMat = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   cParams->bMat[n] = (double*)malloc( sNo * sizeof(double) );  }
    // scaling factor for E-step
    cParams->cn = (double*)malloc( dLen * sizeof(double) );    

    // temporary storage of calculation resutls to save time
    cParams->valpZnZn1 = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   cParams->valpZnZn1[i] = (double*)malloc( sNo * sizeof(double) );  }
    cParams->valpXnZn = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   cParams->valpXnZn[n] = (double*)malloc( sNo * sizeof(double) );  }
    
    return cParams;
}

// Frees memory for common parameters
void freeCommonParameters( cParams )
vbHmmCommonParameters *cParams;
{
    size_t dLen = cParams->dLen, sNo = cParams->sNo;
    int i, n;

    // gamma
    for( n = 0 ; n < dLen ; n++ )
    {   free( cParams->gmMat[n] );  }
    free( cParams->gmMat );
    // xi
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            free( cParams->xiMat[n][i] );
        }
        free( cParams->xiMat[n] );
    }
    free( cParams->xiMat );

    // alpha
    for( n = 0 ; n < dLen ; n++ )
    {   free( cParams->aMat[n] );  }
    free( cParams->aMat );
    // beta
    for( n = 0 ; n < dLen ; n++ )
    {   free( cParams->bMat[n] );  }
    free( cParams->bMat );
    // scaling factor
    free( cParams->cn );

    // temporary storage
    for( i = 0 ; i < sNo ; i++ )
    {   free( cParams->valpZnZn1[i] );  }
    free( cParams->valpZnZn1 );
    for( n = 0 ; n < dLen ; n++ )
    {   free( cParams->valpXnZn[n] );  }
    free( cParams->valpXnZn );

    free( cParams );
    cParams = NULL;
}

// Provide blank struct of analysis results
vbHmmResults *newResults()
{
    vbHmmResults *results = (vbHmmResults*)malloc( sizeof(vbHmmResults) );
    results->iteration = 0;
    results->maxLq = 0.0;
    results->LqArr = NULL;
    results->maxSumTraj = NULL;
    return results;
}

// Frees memory for results struct
void freeResults( results )
vbHmmResults *results;
{
    free( results->LqArr );
    free( results->maxSumTraj );
    free( results );
    results = NULL;
}


// Baum-Welch algorithm for E-step calculation
void forwardBackward( xnWv, cParams, params )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
void *params;
{
    size_t dLen = cParams->dLen;	// number of time stamp data points
    int sNo = cParams->sNo;
    double **gmMat = cParams->gmMat;
    double ***xiMat = cParams->xiMat;
    double **aMat = cParams->aMat, **bMat = cParams->bMat;
    double *cn = cParams->cn;
    double **valpZnZn1 = cParams->valpZnZn1, **valpXnZn = cParams->valpXnZn;

    size_t  n, i, j;

    // forward
    cn[0] = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        valpXnZn[0][i] = (*pTilde_xn_zn)( xnWv, 0, (int)i, params );
        aMat[0][i] = (*pTilde_z1)( (int)i, params ) * valpXnZn[0][i];

        cn[0] += aMat[0][i];

        for( j = 0 ; j < sNo ; j++ ){
            valpZnZn1[i][j] = (*pTilde_zn_zn1)( (int)i, (int)j, params );
        }
    }
    for( i = 0 ; i < sNo ; i++ ){
        aMat[0][i] /= cn[0];
    }
    for( n = 1 ; n < dLen ; n++ ){
        cn[n] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            aMat[n][j] = 0;
            for( i = 0 ; i < sNo ; i++ ){
                aMat[n][j] += aMat[n-1][i] * valpZnZn1[i][j];
            }
            valpXnZn[n][j] = (*pTilde_xn_zn)( xnWv, n, (int)j, params );
            aMat[n][j] *= valpXnZn[n][j];
            cn[n] += aMat[n][j];
        }
        for( j = 0 ; j < sNo ; j++ ){
            aMat[n][j] /= cn[n];
        }
    }

    // backward
    for( i = 0 ; i < sNo ; i++ ){
        bMat[dLen-1][i] = 1;
    }
    double betaTerm;
    for( n = dLen-1 ; n > 0 ; ){
        n--;
        for( i = 0 ; i < sNo ; i++ ){
            bMat[n][i] = 0;
            for( j = 0 ; j < sNo ; j++ ){
                betaTerm  = bMat[n+1][j];
                betaTerm *= valpZnZn1[i][j];
                betaTerm *= valpXnZn[n+1][j];
                bMat[n][i] += betaTerm;
            }
            bMat[n][i] /= cn[n+1];
        }
    }

    // update gamma
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] = aMat[n][i] * bMat[n][i];
        }
    }

    // update xi
    double xiTerm;
    for( i = 0 ; i < sNo ; i++ ){
        for( j = 0 ; j < sNo ; j++ ){
            xiMat[0][i][j] = 0;
    }   }
    for( n = 1 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            for( j = 0 ; j < sNo ; j++ ){
                xiTerm  = aMat[n-1][i];
                xiTerm *= valpXnZn[n][j];
                xiTerm *= valpZnZn1[i][j];
                xiTerm *= bMat[n][j];
                xiMat[n][i][j] = xiTerm / cn[n];
            }
        }
    }
}


// Viterbi algorithm to construct most likely trajectory
int *maxSum( xnWv, cParams, params, maxSumTraj )
xnDataSet *xnWv;
vbHmmCommonParameters *cParams;
int **maxSumTraj;
void *params;
{
    size_t dLen = cParams->dLen;
    int sNo = cParams->sNo;
    size_t n;
    int i, j;

    *maxSumTraj = (int*)realloc( *maxSumTraj, dLen * sizeof(int) );
    double **wnMat = (double **)malloc( dLen * sizeof(double*) );
    double **phiMat = (double **)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ ){
        wnMat[n] = (double*)malloc( sNo * sizeof(double) );
        phiMat[n] = (double*)malloc( sNo * sizeof(double) );
    }
    int maxI;
    double wnTest, maxWn;

    // forward
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            wnMat[n][i] = 0.0;
            phiMat[n][i] = 0.0;
        }
    }
    for( i = 0 ; i < sNo ; i++ ){
        wnMat[0][i] = log((*pTilde_z1)(i, params)) + log((*pTilde_xn_zn)(xnWv, 0, i, params));
    }
    for( n = 1 ; n < dLen ; n++ ){
        for( j = 0 ; j < sNo ; j++ ){
            maxWn = log( (*pTilde_zn_zn1)(0, j, params) ) + wnMat[n-1][0];
            maxI = 0;
            for( i = 1 ; i < sNo ; i++ ){
                wnTest = log( (*pTilde_zn_zn1)(i, j, params) ) + wnMat[n-1][i];
                if( wnTest > maxWn ){
                    maxWn = wnTest;
                    maxI = i;
                }
            }
            phiMat[n][j] = maxI;
            wnMat[n][j] = log((*pTilde_xn_zn)(xnWv, n, j, params)) + maxWn;
        }
    }

    // backward
    maxWn = wnMat[dLen-1][0];
    maxI = 0;
    for( i = 1 ; i < sNo ; i++ ){
        if( wnMat[dLen-1][i] > maxWn ){
            maxWn = wnMat[dLen-1][i];
            maxI = i;
        }
    }
    (*maxSumTraj)[dLen-1] = maxI;
    for( n = dLen-1 ; n > 0 ; n-- ){
        (*maxSumTraj)[n-1] = phiMat[n][(*maxSumTraj)[n]];
    }

    for( n = 0 ; n < dLen ; n++ ){
        free( wnMat[n] );
        free( phiMat[n] );
    }
    free( wnMat );
    free( phiMat );
    
    return *maxSumTraj;
}


//
