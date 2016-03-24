/*
 *  vbHmmTsFretDataHandler.c
 *  File loader for VB-HMM-PC-FRET data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmPcFretDataHandler.h"

#define DATA_READ_NUMBER 1000

xnDataSet *readPcFretBinary( filename1, filename2, out_name, logFP )
char *filename1;
char *filename2;
char *out_name;
FILE *logFP;
{
    size_t dlen1, dlen2;
    double acqTime1, acqTime2;

    unsigned short *shortData1 = readPcBinary( filename1, &dlen1, &acqTime1, logFP );
    unsigned short *shortData2 = readPcBinary( filename2, &dlen2, &acqTime2, logFP );

    if( (dlen1 != dlen2) || (acqTime1 != acqTime2) ){
        fprintf( logFP, "  Data length and/or bin size doesn't match between donor and acceptor signals.\n" );
        fprintf( logFP, "      data length : D) %d,  A)%d \n", (int)dlen1, (int)dlen2 );
        fprintf( logFP, "      bin size : D) %g,  A)%g \n", acqTime1, acqTime2 );
        exit(1);
    }

    xnDataSet *xn = newXnDataSet_pcFret( out_name );
    xn->N = dlen1 - 2;
    xn->data = (pcFretData*)malloc( sizeof( pcFretData ) );
    pcFretData *d = xn->data;
    d->binSize = acqTime1;
    d->dCounts = (unsigned int*)malloc( (dlen1-2) * sizeof(unsigned int) );
    d->aCounts = (unsigned int*)malloc( (dlen1-2) * sizeof(unsigned int) );
    size_t i;
    for( i = 2 ; i < dlen1 ; i++ ){
        d->dCounts[i-2] = shortData1[i];
        d->aCounts[i-2] = shortData2[i];
    }
#ifdef SIMU_STATS
    pc->states = NULL;
    pc->pn = NULL;
#endif

    fprintf( logFP, "  Total of %u bins, with bin size %g s, %.3lf seconds read in.\n\n", (unsigned int)xn->N, acqTime1, acqTime1 * xn->N);

    free( shortData1 );
    free( shortData2 );
    return xn;
}


unsigned short *readPcBinary( filename, dlen, acqTime, logFP )
char *filename;
size_t *dlen;
double *acqTime;
FILE *logFP;
{
    FILE *fp;
    if( (fp = fopen( filename, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( logFP, "File [%s] does not exist. Exiting ...\n", filename );
        exit(1);
        
    }

    // file exists. start reading in data and correlate
    size_t dataRead, totalRead = 0;
    unsigned char *buffer = NULL, *data= NULL;
    size_t i;
    fprintf( logFP, "  Reading data from '%s' ... ", filename );

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
    
    *dlen = (totalRead/2);
    unsigned short *shortData = (unsigned short*)malloc( *dlen * sizeof(unsigned short) );
    for( i = 0 ; i < *dlen ; i++ ){
        shortData[i] = 256*(unsigned short)data[2*i] + (unsigned short)data[2*i+1];
    }
    free(data);     data = NULL;
    
    *acqTime= (shortData[0]*65536 + shortData[1]) * 1.0e-9;
    fprintf( logFP, "done.\n" );
    return shortData;
}    

//

