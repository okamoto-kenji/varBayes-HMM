/*
 *  vbHmmTsDataHandler.c
 *  File loader for VB-HMM-TS data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmTsDataHandler.h"

#define DATA_READ_NUMBER 1000

xnDataSet *readTimeStampBinary( filename, logFP )
char *filename;
FILE *logFP;
{
    FILE *fp;
    if( (fp = fopen( filename, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( stderr, "File [%s] does not exist. Exiting ...\n", filename );
        exit(1);
        
    } else {
        // file exists. start reading in data and correlate
        fprintf( logFP, "  Reading data from '%s'.\n", filename );
        size_t dataRead, totalRead = 0;
        unsigned char *buffer, *data;
        size_t i;
        buffer = (unsigned char*)malloc(DATA_READ_NUMBER * sizeof(unsigned char));
        data = NULL;
        do{
            dataRead = fread( buffer, sizeof(unsigned char), DATA_READ_NUMBER, fp);
            data = (unsigned char*)realloc( data, sizeof(unsigned char)*(totalRead+dataRead) );
            for( i = 0 ; i < dataRead ; i++ ){
                data[totalRead+i] = buffer[i];
            }
            totalRead += dataRead;
        }while( dataRead >= DATA_READ_NUMBER );
        fclose(fp);
        free(buffer);
        
        size_t dlen = (totalRead/4);
        unsigned long *longData = (unsigned long*)malloc( dlen * sizeof(unsigned long) );
        for( i = 0 ; i < dlen ; i++ ){
            longData[i] = 16777216*(long)data[i*4] + 65536*(long)data[i*4+1] + 256*(long)data[i*4+2] + (long)data[i*4+3];
        }
        free(data);
        
        xnDataSet *traj = tsDataSetFromLongArray( longData, dlen, logFP );

        free( longData );        
        return traj;
    }
}


xnDataSet *tsDataSetFromLongArray( longData, dlen, logFP )
unsigned long *longData;
size_t dlen;
FILE *logFP;
{
    if( longData == NULL ){
        return NULL;
    }
    
    double freq = (double)longData[0];
    fprintf( logFP, "  sample freq:%g\n", freq);
    xnDataSet *x = (xnDataSet*)malloc( sizeof(xnDataSet) );
    x->data = (tsData*)malloc( (dlen-1) * sizeof( tsData ) );
    tsData *data = x->data;

    double dt, T = 0.0;
    size_t count, n;
    for( n = 1 ; n < dlen ; n++ ){
        count = longData[n];
        
        dt = (double)count/freq;
        T += dt;
        data[n-1].dt = dt;
        data[n-1].time = T;
    }
    x->N = n-1;
    x->T = T;
    fprintf( logFP, "  Total of %u photons, %.3lf seconds read in.\n\n", (unsigned int)x->N, x->T);
    
    return x;
}

//
