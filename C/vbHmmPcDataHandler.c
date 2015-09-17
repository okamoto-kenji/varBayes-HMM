/*
 *  vbHmmPcDataHandler.c
 *  File loader for VB-HMM-PC data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmPcDataHandler.h"

#define DATA_READ_NUMBER 1000

xnDataSet *readPcBinary( filename, logFP )
char *filename;
FILE *logFP;
{
    FILE *fp;
    if( (fp = fopen( filename, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( stderr, "File [%s] does not exist. Exiting ...\n", filename );
        return NULL;

    }

    // file exists. start reading in data and correlate
    fprintf( logFP, "  Reading data from '%s'.\n", filename );
    size_t dataRead, totalRead = 0;
    unsigned char *buffer, *data = NULL;
    size_t i;
    buffer = (unsigned char*)malloc(DATA_READ_NUMBER * sizeof(unsigned char));
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

    size_t dlen = (totalRead/2);
    unsigned short *shortData = (unsigned short*)malloc( dlen * sizeof(unsigned short) );
    for( i = 0 ; i < dlen ; i++ ){
        shortData[i] = 256*(unsigned short)data[2*i] + (unsigned short)data[2*i+1];
    }
    free(data);

    double acqTime= (shortData[0]*65536 + shortData[1]) * 1.0e-9;

    xnDataSet *traj = (xnDataSet*)malloc( sizeof(xnDataSet) );
    traj->N = dlen - 2;
    traj->T = acqTime * traj->N;
    traj->data = (pcData*)malloc( sizeof( pcData ) );
    pcData *pc = traj->data;
    pc->binSize = acqTime;
    pc->counts = (unsigned int*)malloc( (dlen-2) * sizeof(unsigned int) );
#ifdef SIMU_STATS
    pc->states = NULL;
    pc->pn = NULL;
#endif

    for( i = 2 ; i < dlen ; i++ ){
        pc->counts[i-2] = shortData[i];
    }
    fprintf( logFP, "  Total of %u bins, %.3lf seconds read in.\n\n", (unsigned int)traj->N, traj->T);

    free( shortData );        
    return traj;
}


//
