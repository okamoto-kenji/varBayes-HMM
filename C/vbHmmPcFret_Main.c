/*
 *  vbHmmPcFret_Main.c
 *  Model-specific main function for VB-HMM-PC-FRET.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2011.04.19
 */

#include "vbHmmPcFret_Main.h"
#include "vbHmmPcFret.h"
#include "vbHmmPcFretDataHandler.h"
#include <string.h>
#include <time.h>
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

int main( int argc, char *argv[] ){

    char  in_name1[255], in_name2[255], out_name[255];
    FILE  *logFP = stderr;
    char  logFilename[256];
    time_t  startTime = time((time_t *)NULL);
    initRan();

    int sFrom = 0, sTo = 0, trials = 0, maxIteration = 0;
    double threshold = 0.0;

    // Get command line arguments.
    if( argc != 8 ){
        // Invalid number of arguments
        fprintf( stderr, "\nVariational-Bayesian Hidden-Markov-Model Photon-Counting FRET Analysis\n" );
        fprintf( stderr, "    built on 2010.09.24\n" );
        fprintf( stderr, "Syntax : vbHmmPcFret sFrom sMax trials maxIteration threshold donor_filename acceptor_filename\n" );
        fprintf( stderr, "         sFrom        : minimum number of state to analyze\n" );
        fprintf( stderr, "         sMax         : maximum number of state to analyze\n" );
        fprintf( stderr, "         trials       : number of repetition\n" );
        fprintf( stderr, "                        (giving chance to change initial conditions)\n" );
        fprintf( stderr, "         maxIteration : maximum iteration if inference does not reach threshold\n" );
        fprintf( stderr, "         threshold    : stop iteration if inference reaches this value\n" );
        fprintf( stderr, "         d_filename   : name of data file for donor signal\n" );
        fprintf( stderr, "         a_filename   : name of data file for acceptor signal\n" );
        fprintf( stderr, "Example: vbHmmPcFret 2 20 5 1000 1e-10 D-data A-data\n\n" );
        freeRan();
        exit(1);
        
    } else if( argc == 8 ){
        // Set parameters
        sFrom = atoi( argv[1] );
        sTo = atoi( argv[2] );
        trials = atoi( argv[3] );
        maxIteration = atoi( argv[4] );
        threshold = atof( argv[5] );

        // Set the output filename
        // assumes two files, for donor and acceptor signals, respectively, begin with
        // the same strings and different suffxes.
        strncpy( in_name1, argv[6], sizeof(in_name1) );
        strncpy( in_name2, argv[7], sizeof(in_name2) );
        int i;
        for( i = 0 ; (in_name1[i]!='\0')&&(in_name2[i]!='\0') ; i++ ){
            if( in_name1[i]==in_name2[i] ){
                out_name[i] = in_name1[i];
            } else {
                break;
            }
        }
        if( (i > 0) && (out_name[i-1] == '_') )    i--;
        if( i == 0 ){
            strncpy( out_name, "output", sizeof(out_name) );
        } else {
            out_name[i] = '\0';
        }
        
        strncpy( logFilename, out_name, sizeof(logFilename) );
        strncat( logFilename, ".log", sizeof(logFilename) - strlen(logFilename) - 1 );
        logFP = fopen( logFilename, "w" );
    }

    // Prepare data, typically by loading from file(s).
    xnDataSet *xnWv = readPcFretBinary( in_name1, in_name2, logFP );

    // If data can be obtained, execute analysis.
    if( xnWv != NULL ){
        setFunctions_pcFret();
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
