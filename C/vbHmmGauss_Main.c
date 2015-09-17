/*
 *  vbHmmGauss_Main.c
 *  Model-specific main function for VB-HMM-GAUSS.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmGauss_Main.h"
#include "vbHmmGauss.h"
#include "vbHmmGaussDataHandler.h"
#include <string.h>
#include <time.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

int main( int argc, char *argv[] ){

    char  in_name[255], out_name[255];
    FILE  *logFP = stderr;
    char  logFilename[256];
    time_t  startTime = time((time_t *)NULL);
    initRan();

    int sFrom = 0, sTo = 0, trials = 0, maxIteration = 0;
    double threshold = 0.0;

    // Get command line arguments.
    if( argc != 7 ){
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model Gauss Analysis\n" );
        fprintf( stderr, "    built on 2015.09.xx\n" );
        fprintf( stderr, "Syntax : vbHmmGauss sFrom sMax trials maxIteration threshold data_filename\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition\n" );
        fprintf( stderr, "                        (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         data_filename   : name of data file\n" );
        fprintf( stderr, "Example: vbHmmGauss 2 20 5 1000 1e-10 data-file\n\n" );
        freeRan();
        exit(1);
        
    } else if( argc == 7 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );

        strncpy( in_name, argv[6], sizeof(in_name) );
        strncpy( out_name, argv[6], sizeof(out_name) );

        strncpy( logFilename, out_name, sizeof(logFilename) );
        strncat( logFilename, ".log", sizeof(logFilename) - strlen(logFilename) - 1 );
        logFP = fopen( logFilename, "w" );
    }

    // Prepare data, typically by loading from file(s).
    xnDataSet *xnWv = readGaussText( in_name, logFP );

    // If data can be obtained, execute analysis.
    if( xnWv != NULL ){
        setFunctions_gauss();
        int optK = modelComparison( xnWv, sFrom, sTo, trials, maxIteration, threshold, out_name, logFP );
        fprintf( logFP, " No. of state = %d was chosen.\n", optK);

        free( xnWv->data );
        free( xnWv );
    }
    
    fprintf( logFP, "FINISH: %d sec. spent.\n", (int)(time((time_t *)NULL) - startTime) );
    if( logFP != stderr )  fclose( logFP );

    freeRan();
    return 1;
}


//
