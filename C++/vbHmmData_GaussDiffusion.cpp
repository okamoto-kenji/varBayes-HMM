/*
 *  vbHmmData_GaussDiffusion.c
 *  Model-specific functions for VB-HMM-GAUSS-DIFFUSION.
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
#include <cstdlib>
#include <fstream>
#include <string>
#include "vbHmm_GaussDiffusion.h"


//////////  trace
vbHmmTrace_GaussDiff::vbHmmTrace_GaussDiff() : vbHmmTrace() {
    xn = NULL;
}

vbHmmTrace_GaussDiff::vbHmmTrace_GaussDiff(const vbHmmTrace_GaussDiff &other) : vbHmmTrace( other ) {  // do not copy
    xn = NULL;
}

vbHmmTrace_GaussDiff::~vbHmmTrace_GaussDiff(){
    delete xn;
}

double vbHmmTrace_GaussDiff::kMeans_x( int n ){
    return xn->p[n];
}


void vbHmmTrace_GaussDiff::readGaussDiffText( const char *filename, fstream *logFS ){
    int dlen;
    
    Vec1D *doubleData = readDoubleArrayFromFile( filename, &dlen, logFS );
    if( doubleData->d1 > 0 ){
        name = filename;
        N = doubleData->d1;
        xn = doubleData;

        *logFS << "  Total of " << N << " points read in." << endl << endl;
    } else {
        delete doubleData;
        *logFS << "  No data was read from file '" << filename << "'." << endl << endl;
    }
}


Vec1D *vbHmmTrace_GaussDiff::readDoubleArrayFromFile( const char *filename, int *dlen, fstream *logFS ){
    fstream fs;
    fs.open( filename, ios::in );
    if( !fs.is_open() ){
        // file does not exist, exit with error
        *logFS << "File [" << filename << "] does not exist. Exiting ..." << endl;
        exit(1);
    }
    
    // file exists. start reading in data and correlate
    *logFS << "  Reading data from '" << filename << "' ... ";
    
    Vec1D *doubleData = new Vec1D( 0 );
    double data;
    string str;
    while(1){
        fs >> str;
        if( fs.eof() )  break;
        data = atof( str.c_str() );
        doubleData->append( data );
    }
    fs.close();

    *logFS << "done." << endl;

    return doubleData;
}


//
