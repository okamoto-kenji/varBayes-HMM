/*
 *  vbHmmTsFretDataHandler.c
 *  File loader for VB-HMM-TS-FRET data.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.0.0
 *  Last modified on 2015.09.17
 */

#include "vbHmmTsFretDataHandler.h"

#define DATA_READ_NUMBER 1000

xnDataSet *readTsFretBinary( filename1, filename2, out_name, logFP )
char *filename1;
char *filename2;
char *out_name;
FILE *logFP;
{
    FILE *fp;
    size_t dataRead, totalRead;
    unsigned char *buffer, *data;
    size_t i, dlen1, dlen2;

    if( (fp = fopen( filename1, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( logFP, "File [%s] does not exist. Exiting ...\n", filename1 );
        exit(1);
        
    }
    // file exists. start reading in data and correlate
    fprintf( logFP, "  Reading donor data from '%s' ... ", filename1 );
    totalRead = 0;
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

    dlen1 = (totalRead/4);  // 4 is long variable size of LabVIEW
    unsigned long *longData1 = (unsigned long*)malloc( dlen1 * sizeof(unsigned long) );
    for( i = 0 ; i < dlen1 ; i++ ){
        longData1[i] = 16777216*(long)data[i*4] + 65536*(long)data[i*4+1] + 256*(long)data[i*4+2] + (long)data[i*4+3];
    }
    fprintf( logFP, " %d data points read.\n", (int)dlen1-1 );
    free(data);     data = NULL;
    
    if( (fp = fopen( filename2, "r" )) == NULL ){
        // file does not exist, exit with error
        fprintf( logFP, "File [%s] does not exist. Exiting ...\n", filename2 );
        free( longData1 );
        exit(1);
        
    }
    // file exists. start reading in data and correlate
    fprintf( logFP, "  Reading acceptor data from '%s' ... ", filename2 );
    totalRead = 0;
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

    dlen2 = (totalRead/4);
    unsigned long *longData2 = (unsigned long*)malloc( dlen2 * sizeof(unsigned long) );
    for( i = 0 ; i < dlen2 ; i++ ){
        longData2[i] = 16777216*(long)data[i*4] + 65536*(long)data[i*4+1] + 256*(long)data[i*4+2] + (long)data[i*4+3];
    }
    fprintf( logFP, " %d data points read.\n", (int)dlen2-1 );
    free(data);

    xnDataSet *traj = tsFretDataSetFromLongArray( longData1, longData2, dlen1, dlen2, out_name, logFP );

    free( longData1 );
    free( longData2 );
    return traj;
}


xnDataSet *tsFretDataSetFromLongArray( longData1, longData2, dlen1, dlen2, out_name, logFP )
unsigned long *longData1;
unsigned long *longData2;
size_t dlen1;
size_t dlen2;
char *out_name;
FILE *logFP;
{
    if( (longData1 == NULL)||(longData2==NULL) ){
        return NULL;
    }
    
    double freq = (double)longData1[0], freq2 = (double)longData2[0];
    if( freq == freq2 ){
        fprintf( logFP, "  sample freq:%g\n", freq);
    } else {
        fprintf( logFP, "  frequencies are different:%g/%g\n", freq, freq2);
        return NULL;
    }
    xnDataSet *xn = newXnDataSet_tsFret( out_name );
    tsFretData *d = xn->data;
    d->dt = (double*)malloc( (dlen1 + dlen2-2) * sizeof(double) );
//    d->time = (double*)malloc( (dlen1 + dlen2-2) * sizeof(double) );
    d->ch = (int*)malloc( (dlen1 + dlen2-2) * sizeof(int) );

    size_t c1=0, c2=0, c=0, c0=0;   // total counts
    size_t n1=1, n2=1, n=0;         // index for array
    int ch;
    c1 = longData1[n1];     c2 = longData2[n2];
    for( ; (n1 < dlen1) || (n2 < dlen2) ; n++ ){

        if( n1 >= dlen1 ){
            ch = 1;

        } else if( n2 >= dlen2 ){
            ch = 0;

        } else{
            if( c1 <= c2 ){
                ch = 0;
            } else {
                ch = 1;
            }
        }

        if( ch == 0 ){
            c = c1;
            n1++;
            c1 += longData1[n1];
        } else {
            c = c2;
            n2++;
            c2 += longData2[n2];
        }
        
        d->dt[n]   = (double)(c - c0)/freq;
//        d->time[n] = (double)c/freq;
        d->ch[n]   = ch;
        c0 = c;
    }
    xn->N = n;
//    d->T = d->time[n-1];
    d->T = (double)c/freq;
    fprintf( logFP, "  Total of %u photons, %.3lf seconds read in.\n\n", (unsigned int)xn->N, d->T);

    return xn;
}
