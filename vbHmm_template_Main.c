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
    
    char  in_name[255], out_name[255];
    FILE  *logFP = stderr;
    char  logFilename[256];
    time_t  startTime = time((time_t *)NULL);
    initRan();
    
    int sFrom, sTo, trials, maxIteration;
    double threshold;

    // Get command line arguments.
    if (argc != 7) {
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model\n" );
//        fprintf( stderr, "    built on 2011.03.07\n" );
        fprintf( stderr, "Syntax : (vbHmm) sFrom sMax trials maxIteration threshold filename\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         filename     : name of data file\n" );
        fprintf( stderr, "Example: (vbHmm) 2 20 5 1000 1e-10 data\n\n" );
        freeRan();
        exit(1);

    } else if( argc == 7 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );
        
        // Set the output filename
        strncpy( in_name, argv[6], sizeof(in_name) );
        strncpy( out_name, argv[6], sizeof(out_name) );

        strncpy( logFilename, out_name, sizeof(logFilename) );
        strncat( logFilename, ".log", sizeof(logFilename) );
        logFP = fopen( logFilename, "w" );
    }
    
    // Prepare data, typically by loading from file(s).
    xnDataSet *xnWv = read_template_data( in_name, logFP );

    // If data can be obtained, execute analysis.
    if( xnWv != NULL ){
        setFunctions__();
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
