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
        
        xnDataSet *xn = NULL;
        if( longData != NULL ){
            xn = newXnDataSet_ts( filename );
            tsDataSetFromLongArray( xn, longData, dlen, logFP );
        }
        free( longData );
        return xn;
    }
}


void tsDataSetFromLongArray( xn, longData, dlen, logFP )
xnDataSet *xn;
unsigned long *longData;
size_t dlen;
FILE *logFP;
{
    double freq = (double)longData[0];
    fprintf( logFP, "  sample freq:%g\n", freq);
    tsData *d = xn->data;
    d->dt = (double*)malloc( (dlen-1) * sizeof(double) );
//    d->time = (double*)malloc( (dlen-1) * sizeof(double) );
    d->time = NULL;

    double dt, T = 0.0;
    size_t count, n;
    for( n = 1 ; n < dlen ; n++ ){
        count = longData[n];
        
        dt = (double)count/freq;
        T += dt;
        d->dt[n-1] = dt;
//        d->time[n-1] = T;
    }
    xn->N = n-1;
    d->T = T;
    fprintf( logFP, "  Total of %u photons, %.3lf seconds read in.\n\n", (unsigned int)xn->N, d->T);
    
    return;
}

//
