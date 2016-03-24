/*
 *  vbHmm_template_Main.c
 *  Template of model-specific main function.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include "vbHmm_template_Main.h"
#include "vbHmm_Template_DataHandler.h"
#include <string.h>
#include <time.h>


// For random number generation.
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

int main( int argc, char *argv[] ){
    
    FILE  *logFP = stderr;
    char  logFilename[256];
    time_t  startTime = time((time_t *)NULL);
    initRan();
    
    int sFrom, sTo, trials, maxIteration;
    double threshold;

    // Get command line arguments.
    if (argc < 7) {
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model\n" );
        fprintf( stderr, "Syntax : (vbHmm) sFrom sMax trials maxIteration threshold filename ...\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         filename     : name of data file\n" );
        fprintf( stderr, "Example: (vbHmm) 2 20 5 1000 1e-10 data\n\n" );
        freeRan();
        exit(1);

    } else if( argc >= 7 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );
        
        // Set the output filename
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
        xnDataSet *xn = read_template_data( argv[i], logFP );
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
            
            setFunctions__();
            xnDataSet *xn = xns->xn[0];
            free( xns->xn );
            free( xns );
            optK = modelComparison( xn, sFrom, sTo, trials, maxIteration, threshold, logFP );
            freeXnDataSet__( &xn );
            
        } else {
            
            setGFunctions__();
            optK = gModelComparison( xns, sFrom, sTo, trials, maxIteration, threshold, logFP );
            for( i = 0 ; i < xns->R ; i++ ){
                freeXnDataSet__( &xns->xn[i] );
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
