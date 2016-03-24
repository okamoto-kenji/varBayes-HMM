/*
 *  vbHmmGaussDataHandler.c
 *  File loader for VB-HMM-GAUSS data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmGaussDataHandler.h"
#include <string.h>

#define DATA_READ_NUMBER 1000

xnDataSet *readGaussText( filename, logFP )
char *filename;
FILE *logFP;
{
    size_t dlen;
    xnDataSet *xn = NULL;

    double *doubleData = readDoubleArrayFromFile( filename, &dlen, logFP );
    if( dlen > 0 ){
        xn = newXnDataSet_gauss( filename );
        xn->N = dlen;
        gaussData *dat = xn->data;
        dat->v = (double*)malloc( dlen * sizeof(double) );
        size_t i;
        for( i = 0 ; i < dlen ; i++ ){
            dat->v[i] = doubleData[i];
        }
        fprintf( logFP, "  Total of %u points read in.\n\n", (unsigned int)xn->N);

        free( doubleData );
    }
    return xn;
}


double *readDoubleArrayFromFile( filename, dlen, logFP )
char *filename;
size_t *dlen;
FILE *logFP;
{
    FILE *fp;
    if( (fp = fopen( filename, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( logFP, "File [%s] does not exist. Exiting ...\n", filename );
        exit(1);
    }

    // file exists. start reading in data and correlate
    size_t totalRead = 0;
    fprintf( logFP, "  Reading data from '%s' ... ", filename );

    double *doubleData = NULL;
    double data;
    size_t slen;
    char *str = (char*)malloc(DATA_READ_NUMBER * sizeof(char));
    while(1){
        str = fgets( str, DATA_READ_NUMBER, fp);
        if( feof(fp) )  break;
        if( (slen = strlen(str)) <= 0 )  break;
        if( str[slen-1] == '\n' )  str[slen-1] = '\0';
        if( (strlen(str) <= 0) )  break;
        data = atof(str);
        totalRead++;
        doubleData =(double*)realloc( doubleData, totalRead * sizeof(double) );
        doubleData[totalRead-1] = data;
    }
    fclose(fp);
    free(str);
    
    *dlen = totalRead;

    fprintf( logFP, "done.\n" );
    return doubleData;
}    

//

