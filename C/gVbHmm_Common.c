/*
 *  gVbHmm_Common.c
 *  Common VB-HMM engine for global analysis.
 *  Reference: 
 *    Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2016.11.04
 */

#include "gVbHmm_Common.h"
#include <string.h>
#include <float.h>


// Pointers of functions to call model-specific functions.
extern new_model_parameters_func   newModelParameters;
extern free_model_parameters_func  freeModelParameters;
extern new_model_stats_func        newModelStats;
extern free_model_stats_func       freeModelStats;
new_model_statsG_func              newModelStatsG          = NULL;
free_model_statsG_func             freeModelStatsG         = NULL;
initialize_vbHmmG_func             initializeVbHmmG        = NULL;
extern pTilde_z1_func              pTilde_z1;
extern pTilde_zn_zn1_func          pTilde_zn_zn1;
extern pTilde_xn_zn_func           pTilde_xn_zn;
calcStatsVarsG_func                calcStatsVarsG          = NULL;
maximizationG_func                 maximizationG           = NULL;
varLowerBoundG_func                varLowerBoundG          = NULL;
reorderParametersG_func            reorderParametersG      = NULL;
outputResultsG_func                outputResultsG          = NULL;


// This function must be called to connect with the model before executing analysis.
void setGFunctions( funcs )
gCommonFunctions funcs;
{
    newModelParameters      = funcs.newModelParameters;
    freeModelParameters     = funcs.freeModelParameters;
    newModelStats           = funcs.newModelStats;
    freeModelStats          = funcs.freeModelStats;
    newModelStatsG          = funcs.newModelStatsG;
    freeModelStatsG         = funcs.freeModelStatsG;
    initializeVbHmmG        = funcs.initializeVbHmmG;
    pTilde_z1               = funcs.pTilde_z1;
    pTilde_zn_zn1           = funcs.pTilde_zn_zn1;
    pTilde_xn_zn            = funcs.pTilde_xn_zn;
    calcStatsVarsG          = funcs.calcStatsVarsG;
    maximizationG           = funcs.maximizationG;
    varLowerBoundG          = funcs.varLowerBoundG;
    reorderParametersG      = funcs.reorderParametersG;
    outputResultsG          = funcs.outputResultsG;
}


//////////////////////////////////////////////////////////////////  VB-HMM Execution Functions
int gModelComparison( xns, sFrom, sTo, trials, maxIteration, threshold, logFP )
xnDataBundle *xns;
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
    
    int maxS = 0, rNo = xns->R;
    double maxLq = -DBL_MAX;
    globalVars **gvArray = (globalVars**)malloc( trials * sizeof(globalVars*) );
    indVarBundle **ivsArray = (indVarBundle**)malloc( trials * sizeof(indVarBundle*) );
    for( s = sFrom ; s <= sTo ; s++ ){
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
        for( t = 0 ; t < trials ; t++ ){
            int r, st = (s - sFrom) * trials + t;
            gvArray[t] = newGlobalVarsG( xns, s );

            ivsArray[t] = (indVarBundle*)malloc( sizeof(indVarBundle) );
            ivsArray[t]->indVars = (indVars**)malloc( rNo * sizeof(indVars*) );
            for( r = 0 ; r < rNo ; r++ ){
                ivsArray[t]->indVars[r] = newIndVars( xns->xn[r], gvArray[t] );
            }
            ivsArray[t]->stats = (*newModelStatsG)( xns, gvArray[t], ivsArray[t] );

            LqVsK[st] = gVbHmm_Main( xns, gvArray[t], ivsArray[t], maxIteration, threshold, logFP );
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
        (*outputResultsG)( xns, gvArray[maxT], ivsArray[maxT], logFP );
        
        for( t = 0 ; t < trials ; t++ ){
            int r;
            for( r = 0 ; r < rNo ; r++ ){
                freeIndVars( xns->xn[r], gvArray[t], &ivsArray[t]->indVars[r] );
            }
            free( ivsArray[t]->indVars );
            free( ivsArray[t]->stats );
            free( ivsArray[t] );
            ivsArray[t] = NULL;
            freeGlobalVarsG( xns, &gvArray[t] );
        }
        
        if( s >= (maxS+3) ){
            s++;
            break;
        }
        
    }
    sTo = s - 1;
    free( gvArray );
    free( ivsArray );

    char fn[256];
    FILE *fp = NULL;
    strncpy( fn, xns->xn[0]->name, sizeof(fn) );
    strncat( fn, ".g.LqVsK", sizeof(fn) - strlen(fn) - 1 );
    if( (fp = fopen( fn, "w" )) != NULL ){
        for( s = 0 ; s < trials * (sTo - sFrom + 1) ; s++ ){
            fprintf( fp, "%2d, %.20g\n", (s/trials) + sFrom, LqVsK[s] );
        }
        fclose(fp);
    }
    free( LqVsK );
    return maxS;
}


// Initializes parameters commonly used in VB-HMM analysis
globalVars *newGlobalVarsG( xns, sNo )
xnDataBundle *xns;
int sNo;
{
    globalVars *gv = (globalVars*)malloc( sizeof(globalVars) );
    gv->sNo = sNo;
    gv->iteration = 0;
    gv->maxLq = 0.0;
    gv->LqArr = NULL;
    gv->params = (*newModelParameters)( NULL, sNo );
    return gv;
}

// Frees memory for common parameters
void freeGlobalVarsG( xns, gv )
xnDataBundle *xns;
globalVars **gv;
{
    free( (*gv)->LqArr );
    (*freeModelParameters)( &(*gv)->params, NULL, (*gv)->sNo );
    free( *gv );
    *gv = NULL;
}


//////////////////////////////////////////////////////////////////  VB-HMM Common Engine
double gVbHmm_Main( xns, gv, ivs, maxIteration, threshold, logFP )
xnDataBundle *xns;
globalVars *gv;
indVarBundle *ivs;
int maxIteration;
double threshold;
FILE *logFP;
{
    double **LqArr = &gv->LqArr;
    *LqArr = realloc( *LqArr, maxIteration * sizeof(double) );

    (*initializeVbHmmG)( xns, gv, ivs );

    int i, r, rNo = xns->R;
    for( i = 0 ; i < maxIteration ; i++ ){

        // E-step
        for( r = 0 ; r < rNo ; r++ ){
            forwardBackward( xns->xn[r], gv, ivs->indVars[r] );
        }
        (*calcStatsVarsG)( xns, gv, ivs );
        (*LqArr)[i] = (*varLowerBoundG)( xns, gv, ivs );
        
        // End loop if derivative of variational lower bound reaches threshold.
        if( (i>0) && ( fabs( ((*LqArr)[i] - (*LqArr)[i-1]) / (*LqArr)[i] ) < threshold ) ){
            break;
        }
        
        // M-step
        (*maximizationG)( xns, gv, ivs );
    }
    if( i == maxIteration ){
        if( logFP != NULL ){
            fprintf(logFP, "MAX iteration (%d) reached.\n", maxIteration);
        }
        i--;
    }
    (*reorderParametersG)( xns, gv, ivs );
    for( r = 0 ; r < rNo ; r++ ){
#ifdef OUTPUT_MAX_GAMMA
        maxGamma( xns->xn[r], gv, ivs->indVars[r] );
#endif
        maxSum( xns->xn[r], gv, ivs->indVars[r] );
    }

    gv->iteration = i+1;
    *LqArr = realloc( *LqArr, (i+1) * sizeof(double) );
    gv->maxLq = (*LqArr)[i];
    if( logFP != NULL ){
        fprintf( logFP, "  iteration: %d    evidence p(x|K=%d) = %.20g \n", i+1, gv->sNo, gv->maxLq );
    }
    return gv->maxLq;
}


//
