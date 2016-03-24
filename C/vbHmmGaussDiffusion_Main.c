/*
 *  vbHmmGaussDiffusion_Main.c
 *  Model-specific main function for VB-HMM-GAUSS-DIFFUSION.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2016
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2016.02.08
 */

#include "vbHmmGaussDiffusion_Main.h"
#include "vbHmmGaussDiffusion.h"
#include "vbHmmGaussDiffusionDataHandler.h"
#include <string.h>
#include <time.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

int main( int argc, char *argv[] ){

    FILE  *logFP = stderr;
    char  logFilename[256];
    time_t  startTime = time((time_t *)NULL);
    initRan();

    int sFrom = 0, sTo = 0, trials = 0, maxIteration = 0;
    double threshold = 0.0;

    // Get command line arguments.
    if( argc < 7 ){
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model Gauss-Diffusion Analysis\n" );
        fprintf( stderr, "    built on 2015.09.xx\n" );
        fprintf( stderr, "Syntax : vbHmmGaussDiff sFrom sMax trials maxIteration threshold data_filename ...\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition\n" );
        fprintf( stderr, "                        (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         data_filename   : name of data file\n" );
        fprintf( stderr, "Example: vbHmmGaussDiff 2 20 5 1000 1e-10 data-file\n\n" );
        freeRan();
        exit(1);
        
    } else if( argc >= 7 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );

        strncpy( logFilename, argv[6], sizeof(logFilename) );
        strncat( logFilename, ".log", sizeof(logFilename) - strlen(logFilename) - 1 );
        logFP = fopen( logFilename, "w" );
    }

    // Prepare data, typically by loading from file(s).
    xnDataBundle *xns = (xnDataBundle*)malloc( sizeof(xnDataBundle) );
    xns->R = 0;
    xns->xn = NULL;
    int i;
    for( i = 6 ; i < argc ; i++ ){
        xnDataSet *xn = readGaussDiffText( argv[i], logFP );
        if( xn != NULL ){
            xns->xn = (xnDataSet**)realloc( xns->xn, (xns->R + 1) * sizeof(xnDataSet*) );
            xns->xn[xns->R] = xn;
            xns->R++;
        }
    }

    // If data can be obtained, execute analysis.
    if( xns->R > 0 ){
        int optK;
        if( xns->R == 1 ){

            setFunctions_gaussDiff();
            xnDataSet *xn = xns->xn[0];
            free( xns->xn );
            free( xns );
            optK = modelComparison( xn, sFrom, sTo, trials, maxIteration, threshold, logFP );
            freeXnDataSet_gaussDiff( &xn );
            
        } else {
            
            setGFunctions_gaussDiff();
            optK = gModelComparison( xns, sFrom, sTo, trials, maxIteration, threshold, logFP );
            for( i = 0 ; i < xns->R ; i++ ){
                freeXnDataSet_gaussDiff( &xns->xn[i] );
            }
            free( xns->xn );
            free( xns );
            
        }
        fprintf( logFP, " No. of state = %d was chosen.\n", optK);
    }
    
    fprintf( logFP, "FINISH: %d sec. spent.\n", (int)(time((time_t *)NULL) - startTime) );
    if( logFP != stderr )  fclose( logFP );

    freeRan();
    return 1;
}


//
