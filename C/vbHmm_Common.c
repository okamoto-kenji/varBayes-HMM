/*
 *  vbHmm_Common.c
 *  Common VB-HMM engine.
 *  Reference: 
 *    Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmm_Common.h"
#include <string.h>
#include <float.h>


// Pointers of functions to call model-specific functions.
new_model_parameters_func  newModelParameters    = NULL;
free_model_parameters_func freeModelParameters   = NULL;
new_model_stats_func       newModelStats         = NULL;
free_model_stats_func      freeModelStats        = NULL;
initialize_vbHmm_func      initializeVbHmm       = NULL;
pTilde_z1_func             pTilde_z1             = NULL;
pTilde_zn_zn1_func         pTilde_zn_zn1         = NULL;
pTilde_xn_zn_func          pTilde_xn_zn          = NULL;
calcStatsVars_func         calcStatsVars         = NULL;
maximization_func          maximization          = NULL;
varLowerBound_func         varLowerBound         = NULL;
reorderParameters_func     reorderParameters     = NULL;
outputResults_func         outputResults         = NULL;


// This function must be called to connect with the model before executing analysis.
void setFunctions( funcs )
commonFunctions funcs;
{
    newModelParameters      = funcs.newModelParameters;
    freeModelParameters     = funcs.freeModelParameters;
    newModelStats           = funcs.newModelStats;
    freeModelStats          = funcs.freeModelStats;
    initializeVbHmm         = funcs.initializeVbHmm;
    pTilde_z1               = funcs.pTilde_z1;
    pTilde_zn_zn1           = funcs.pTilde_zn_zn1;
    pTilde_xn_zn            = funcs.pTilde_xn_zn;
    calcStatsVars           = funcs.calcStatsVars;
    maximization            = funcs.maximization;
    varLowerBound           = funcs.varLowerBound;
    reorderParameters       = funcs.reorderParameters;
    outputResults           = funcs.outputResults;
}


//////////////////////////////////////////////////////////////////  VB-HMM Execution Functions
int modelComparison( xn, sFrom, sTo, trials, maxIteration, threshold, logFP )
xnDataSet *xn;
int sFrom, sTo, trials;
int maxIteration;
double threshold;
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
    globalVars **gvArray = (globalVars**)malloc( trials * sizeof(globalVars*) );
    indVars **ivArray = (indVars**)malloc( trials * sizeof(indVars*) );
    for( s = sFrom ; s <= sTo ; s++ ){
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
        for( t = 0 ; t < trials ; t++ ){
            int st = (s - sFrom) * trials + t;
            gvArray[t] = newGlobalVars( xn, s );
            ivArray[t] = newIndVars( xn, gvArray[t] );
            
            LqVsK[st] = vbHmm_Main( xn, gvArray[t], ivArray[t], maxIteration, threshold, logFP );
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
        (*outputResults)( xn, gvArray[maxT], ivArray[maxT], logFP );
        
        for( t = 0 ; t < trials ; t++ ){
            freeIndVars( xn, gvArray[t], &ivArray[t] );
            freeGlobalVars( xn, &gvArray[t] );
        }
        
        if( s >= (maxS+3) ){
            s++;
            break;
        }
        
    }
    sTo = s - 1;
    free( gvArray );
    free( ivArray );
    
    char fn[256];
    FILE *fp = NULL;
    strncpy( fn, xn->name, sizeof(fn) );
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


globalVars *newGlobalVars( xn, sNo )
xnDataSet *xn;
int sNo;
{
    globalVars *gv = (globalVars*)malloc( sizeof(globalVars) );
    gv->sNo = sNo;
    
    gv->iteration = 0;
    gv->maxLq = 0.0;
    gv->LqArr = NULL;
    
    gv->params = (*newModelParameters)( xn, sNo );
    
    return gv;
}

void freeGlobalVars( xn, gv )
xnDataSet *xn;
globalVars **gv;
{
    free( (*gv)->LqArr );
    
    (*freeModelParameters)( &(*gv)->params, xn, (*gv)->sNo );
    
    free( *gv );
    *gv = NULL;
}

indVars *newIndVars( xn, gv )
xnDataSet *xn;
globalVars *gv;
{
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    int i, n;
    indVars *iv = (indVars*)malloc( sizeof(indVars) );
    
    // gamma
    iv->gmMat = (double**)malloc( dLen * sizeof(double*) );
    double **gmMat = iv->gmMat;
    for( n = 0 ; n < dLen ; n++ ){
        gmMat[n] = (double*)malloc( sNo * sizeof(double) );
        memset( gmMat[n], 0, sNo * sizeof(double) );
    }
    // xi
    iv->xiMat = (double***)malloc( dLen * sizeof(double**) );
    double ***xiMat = iv->xiMat;
    for( n = 0 ; n < dLen ; n++ ){
        xiMat[n] = (double**)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ ){
            xiMat[n][i] = (double*)malloc( sNo * sizeof(double) );
            memset( xiMat[n][i], 0, sNo * sizeof(double) );
        }
    }
    
    // alpha for E-step
    iv->aMat = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   iv->aMat[n] = (double*)malloc( sNo * sizeof(double) );  }
    // beta for E-step
    iv->bMat = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   iv->bMat[n] = (double*)malloc( sNo * sizeof(double) );  }
    // scaling factor for E-step
    iv->cn = (double*)malloc( dLen * sizeof(double) );
    
    // temporary storage of calculation resutls to save time
    iv->valpZnZn1 = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   iv->valpZnZn1[i] = (double*)malloc( sNo * sizeof(double) );  }
    iv->valpXnZn = (double**)malloc( dLen * sizeof(double*) );
    for( n = 0 ; n < dLen ; n++ )
    {   iv->valpXnZn[n] = (double*)malloc( sNo * sizeof(double) );  }
    
    iv->stats = (*newModelStats)( xn, gv, iv );
    
    iv->stateTraj = NULL;
    
    return iv;
}

void freeIndVars( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars **iv;
{
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    int i, n;
    
    // gamma
    for( n = 0 ; n < dLen ; n++ )
    {   free( (*iv)->gmMat[n] );  }
    free( (*iv)->gmMat );
    // xi
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            free( (*iv)->xiMat[n][i] );
        }
        free( (*iv)->xiMat[n] );
    }
    free( (*iv)->xiMat );
    
    // alpha
    for( n = 0 ; n < dLen ; n++ )
    {   free( (*iv)->aMat[n] );  }
    free( (*iv)->aMat );
    // beta
    for( n = 0 ; n < dLen ; n++ )
    {   free( (*iv)->bMat[n] );  }
    free( (*iv)->bMat );
    // scaling factor
    free( (*iv)->cn );
    
    // temporary storage
    for( i = 0 ; i < sNo ; i++ )
    {   free( (*iv)->valpZnZn1[i] );  }
    free( (*iv)->valpZnZn1 );
    for( n = 0 ; n < dLen ; n++ )
    {   free( (*iv)->valpXnZn[n] );  }
    free( (*iv)->valpXnZn );
    
    free( (*iv)->stateTraj );
    
    (*freeModelStats)( &(*iv)->stats, xn, gv, (*iv) );
    
    free( *iv );
    *iv = NULL;
}


//////////////////////////////////////////////////////////////////  VB-HMM Common Engine
double vbHmm_Main( xn, gv, iv ,maxIteration, threshold, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
int maxIteration;
double threshold;
FILE *logFP;
{
    double **LqArr = &gv->LqArr;
    *LqArr = realloc( *LqArr, maxIteration * sizeof(double) );

    (*initializeVbHmm)( xn, gv, iv );

    int i;
    for( i = 0 ; i < maxIteration ; i++ ){
        // E-step
        forwardBackward( xn, gv, iv );

        (*calcStatsVars)( xn, gv, iv );
        (*LqArr)[i] = (*varLowerBound)( xn, gv, iv );

        // End loop if derivative of variational lower bound reaches threshold.
        if( (i>0) && ( fabs( ((*LqArr)[i] - (*LqArr)[i-1]) / (*LqArr)[i] ) < threshold ) ){
            break;
        }

        // M-step
        (*maximization)( xn, gv, iv );
    }
    if( i == maxIteration ){
        if( logFP != NULL ){
            fprintf(logFP, "MAX iteration (%d) reached.\n", maxIteration);
        }
        i--;
    }
    (*reorderParameters)( xn, gv, iv );
    maxSum( xn, gv, iv );

    gv->iteration = i+1;
    *LqArr = realloc( *LqArr, (i+1) * sizeof(double) );
    gv->maxLq = (*LqArr)[i];
    if( logFP != NULL ){
        fprintf( logFP, "  iteration: %d    evidence p(x|K=%d) = %.20g \n", i+1, gv->sNo, gv->maxLq );
    }
    return gv->maxLq;
}


// Baum-Welch algorithm for E-step calculation
void forwardBackward( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    size_t dLen = xn->N;	// number of time stamp data points
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double **aMat = iv->aMat, **bMat = iv->bMat;
    double *cn = iv->cn;
    double **valpZnZn1 = iv->valpZnZn1, **valpXnZn = iv->valpXnZn;

    size_t  n, i, j;

    // forward
    cn[0] = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        valpXnZn[0][i] = (*pTilde_xn_zn)( xn, 0, (int)i, gv->params );
        aMat[0][i] = (*pTilde_z1)( (int)i, gv->params ) * valpXnZn[0][i];

        cn[0] += aMat[0][i];

        for( j = 0 ; j < sNo ; j++ ){
            valpZnZn1[i][j] = (*pTilde_zn_zn1)( (int)i, (int)j, gv->params );
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
            valpXnZn[n][j] = (*pTilde_xn_zn)( xn, n, (int)j, gv->params );
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
int *maxSum( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    size_t n;
    int i, j;

    iv->stateTraj = (int*)realloc( iv->stateTraj, dLen * sizeof(int) );
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
        wnMat[0][i] = log((*pTilde_z1)(i, gv->params)) + log((*pTilde_xn_zn)(xn, 0, i, gv->params));
    }
    for( n = 1 ; n < dLen ; n++ ){
        for( j = 0 ; j < sNo ; j++ ){
            maxWn = log( (*pTilde_zn_zn1)(0, j, gv->params) ) + wnMat[n-1][0];
            maxI = 0;
            for( i = 1 ; i < sNo ; i++ ){
                wnTest = log( (*pTilde_zn_zn1)(i, j, gv->params) ) + wnMat[n-1][i];
                if( wnTest > maxWn ){
                    maxWn = wnTest;
                    maxI = i;
                }
            }
            phiMat[n][j] = maxI;
            wnMat[n][j] = log((*pTilde_xn_zn)(xn, n, j, gv->params)) + maxWn;
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
    iv->stateTraj[dLen-1] = maxI;
    for( n = dLen-1 ; n > 0 ; n-- ){
        iv->stateTraj[n-1] = phiMat[n][iv->stateTraj[n]];
    }

    for( n = 0 ; n < dLen ; n++ ){
        free( wnMat[n] );
        free( phiMat[n] );
    }
    free( wnMat );
    free( phiMat );
    
    return iv->stateTraj;
}


//
