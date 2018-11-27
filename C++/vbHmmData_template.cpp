/*
 *  vbHmmData_template.c
 *  Model-specific functions for VB-HMM-TEMPLATE.
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
#include "vbHmm_template.h"


//////////  trace
vbHmmTrace_template::vbHmmTrace_template() : vbHmmTrace() {
}

vbHmmTrace_template::vbHmmTrace_template(const vbHmmTrace_template &other) : vbHmmTrace( other ) {  // do not copy
}

vbHmmTrace_template::~vbHmmTrace_template(){
}

double vbHmmTrace_template::kMeans_x( int n ){
    return 0.0;
}

void vbHmmTrace_template::importDataFromFile(const char *filename, fstream *logFS){
    // implement function to load data from a file
}

//
