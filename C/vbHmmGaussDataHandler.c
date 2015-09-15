/*
 *  vbHmmGaussDataHandler.c
 *  File loader for VB-HMM-PC-FRET data.
 *
 *  Created by OKAMOTO Kenji and SAKO Yasushi
 *  Copyright 2011
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.xx
 */

#include "vbHmmGaussDataHandler.h"

#define DATA_READ_NUMBER 1000

xnDataSet *readGaussText( filename, logFP )
char *filename;
FILE *logFP;
{
    size_t dlen;

    unsigned short *shortData = readGaussTextLine( filename, &dlen, logFP );

    xnDataSet *traj = (xnDataSet*)malloc( sizeof(xnDataSet) );
    traj->N = dlen;
    traj->T = 0.0;
    traj->data = (gaussData*)malloc( sizeof( gaussData ) );
    gaussData *dat = traj->data;
    dat->v = (double*)malloc( dlen * sizeof(double) );
    size_t i;
    for( i = 0 ; i < dlen ; i++ ){
        dat->v[i] = shortData[i];
    }
    fprintf( logFP, "  Total of %u points read in.\n\n", (unsigned int)traj->N);

    free( shortData );
    return traj;
}


unsigned short *readGaussTextLine( filename, dlen, logFP )
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
//    unsigned char *buffer = NULL;  //, *data= NULL;
//    size_t i;
    fprintf( logFP, "  Reading data from '%s' ... ", filename );

    unsigned short *shortData = NULL;  //(unsigned short*)malloc( *dlen * sizeof(unsigned short) );
    double data;
    unsigned char *buffer = (unsigned char*)malloc(DATA_READ_NUMBER * sizeof(unsigned char));
    do{
//        buffer = fTextRead( buffer, sizeof(unsigned char), DATA_READ_NUMBER, fp);
//        if( eof() )  break;
        data = atof(buffer);
        totalRead++;
        shortData =(unsigned short*)realloc( shortData, totalRead * sizeof(unsigned short) );
//        for( i = 0 ; i < dataRead ; i++ ){
//            data[totalRead+i] = buffer[i];
//        }
//        totalRead += dataRead;
        shortData[totalRead-1] = data;
    }while( eof() );  // dataRead >= DATA_READ_NUMBER );
    fclose(fp);
    free(buffer);
    
    *dlen = totalRead;
//    for( i = 0 ; i < *dlen ; i++ ){
//        shortData[i] = 256*(unsigned short)data[2*i] + (unsigned short)data[2*i+1];
//    }
//    free(data);     data = NULL;
    
    fprintf( logFP, "done.\n" );
    return shortData;
}    

//

