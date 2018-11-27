/*
 *  vbHmmGaussDiffusion_Main.cpp
 *  Model-specific main function for VB-HMM-GAUSS-DIFFUSION.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2018
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2018.11.26
 */

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "vbHmmGaussDiffusion_Main.h"
#include "rand.h"


int main( int argc, char *argv[] ){

    time_t  startTime = time((time_t *)NULL);

    // Get command line arguments.
    if( argc < 7 ){
        // Invalid number of arguments
        cerr << endl << "Variational-Bayesian Hidden-Markov-Model Gauss-Diffusion Analysis" << endl;
        cerr << "    built on 2015.09.xx" << endl;
        cerr << "Syntax : vbHmmGaussDiff sFrom sMax trials maxIteration threshold data_filename ..." << endl;
        cerr << "         sFrom        : minimum number of state to analyze" << endl;
        cerr << "         sMax         : maximum number of state to analyze" << endl;
        cerr << "         trials       : number of repetition" << endl;
        cerr << "                        (giving chance to change initial conditions)" << endl;
        cerr << "         maxIteration : maximum iteration if inference does not reach threshold" << endl;
        cerr << "         threshold    : stop iteration if inference reaches this value" << endl;
        cerr << "         data_filename   : name of data file" << endl;
        cerr << "Example: vbHmmGaussDiff 2 20 5 1000 1e-10 data-file" << endl << endl;
        freeRan();
        exit(1);
        
    }  // else if( argc >= 7 ){

    initRan();
    disableGslErrorHandler();

    // Set parameters
    vbHmmCond_GaussDiff *c = new vbHmmCond_GaussDiff();
    c->sFrom = atoi( argv[1] );
    c->sTo = atoi( argv[2] );
    c->trials = atoi( argv[3] );
    c->maxIteration = atoi( argv[4] );
    c->lbTh = atof( argv[5] );
    c->initType = initType_random;

    fstream  logFS;
    string  logFilename = argv[6];
    logFilename += ".log";
    logFS.open( logFilename, ios::out );

    // Prepare data, typically by loading from file(s).
    vbHmmData *d = new vbHmmData();
    int i;
    for( i = 6 ; i < argc ; i++ ){
        vbHmmTrace_GaussDiff *t = new vbHmmTrace_GaussDiff();
        t->readGaussDiffText( argv[i], &logFS );
        if( t->N > 0 ){
            d->appendTrace( t );
        } else {
            delete t;
        }
    }

    // If data can be obtained, execute analysis.
    if( d->R > 0 ){
        int optK;
        vbHmmModel_GaussDiff *m = new vbHmmModel_GaussDiff();
        optK = modelComparison( c, d, m, &logFS );
        delete m;
    
        logFS <<" No. of state = " << optK << " was chosen." << endl;
    }
    delete d;
    delete c;

    char str[1024];
    sprintf( str, "FINISH: %d sec. spent.", (int)(time((time_t *)NULL) - startTime) );
    logFS << str << endl;
    logFS.close();

    freeRan();
    return 1;
}


//
