/*
 *  vbHmmPc_Main.c
 *  Model-specific main function for VB-HMM-PC.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2015.09.17
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "vbHmmPc_Main.h"
#include "vbHmmPcDataHandler.h"
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
    if( argc != 7 ){
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model Photon-Counting Analysis\n" );
        fprintf( stderr, "    built on 2010.09.24\n" );
        fprintf( stderr, " Syntax : vbHmmPc sFrom sMax trials maxIteration threshold filename\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition\n" );
        fprintf( stderr, "                        (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         filename     : name of data file\n" );
        fprintf( stderr, " Example: vbHmmPc 2 20 5 1000 1e-10 datafile\n\n" );
        freeRan();
        exit(1);
        
    } else if( argc == 7 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );

        // Set the log filename
        strncpy( logFilename, argv[6], sizeof(logFilename) );
        strncat( logFilename, ".log", sizeof(logFilename) - strlen(logFilename) - 1 );
        logFP = fopen( logFilename, "w" );
    }

    // Prepare data, typically by loading from file(s).
    xnDataSet *xn = readPcBinary( argv[6], logFP );

    // If data can be obtained, execute analysis.
    if( xn != NULL ){
        setFunctions_pc();
        int optK = modelComparison( xn, sFrom, sTo, trials, maxIteration, threshold, logFP );
        fprintf( logFP, " No. of state = %d was chosen.\n", optK);

        freeXnDataSet_pc( &xn );
    }
    
    fprintf( logFP, "FINISH: %d sec. spent.\n", (int)(time((time_t *)NULL) - startTime) );
    if( logFP != stderr )  fclose( logFP );

    freeRan();
    return 1;
}

//
