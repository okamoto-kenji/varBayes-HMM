/*
 *  vbHmmPc.c
 *  Model-specific core functions for VB-HMM-PC.
 *
 *  Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
 *  Copyright 2011-2015
 *  Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
 *  All rights reserved.
 *
 *  Ver. 1.1.0
 *  Last modified on 2015.09.17
 */

#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_pow_int.h>
#include <string.h>
#include "vbHmmPc.h"
#include "rand.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define  MAX(a,b)  ((a)>(b)?(a):(b))
#define  MIN(a,b)  ((a)<(b)?(a):(b))

//static int isGlobalAnalysis = 0;

void setFunctions_pc(){
    commonFunctions funcs;
    funcs.newModelParameters    = newModelParameters_pc;
    funcs.freeModelParameters   = freeModelParameters_pc;
    funcs.newModelStats         = newModelStats_pc;
    funcs.freeModelStats        = freeModelStats_pc;
    funcs.initializeVbHmm       = initializeVbHmm_pc;
    funcs.pTilde_z1             = pTilde_z1_pc;
    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_pc;
    funcs.pTilde_xn_zn          = pTilde_xn_zn_pc;
    funcs.calcStatsVars         = calcStatsVars_pc;
    funcs.maximization          = maximization_pc;
    funcs.varLowerBound         = varLowerBound_pc;
    funcs.reorderParameters     = reorderParameters_pc;
    funcs.outputResults         = outputResults_pc;
    setFunctions( funcs );
}

//void setGFunctions_pc(){
//    gCommonFunctions funcs;
//    funcs.newModelParameters    = newModelParameters_pc;
//    funcs.freeModelParameters   = freeModelParameters_pc;
//    funcs.newModelStats         = newModelStats_pc;
//    funcs.freeModelStats        = freeModelStats_pc;
//    funcs.newModelStatsG        = newModelStatsG_pc;
//    funcs.freeModelStatsG       = freeModelStatsG_pc;
//    funcs.initializeVbHmmG      = initializeVbHmmG_pc;
//    funcs.pTilde_z1             = pTilde_z1_pc;
//    funcs.pTilde_zn_zn1         = pTilde_zn_zn1_pc;
//    funcs.pTilde_xn_zn          = pTilde_xn_zn_pc;
//    funcs.calcStatsVarsG        = calcStatsVarsG_pc;
//    funcs.maximizationG         = maximizationG_pc;
//    funcs.varLowerBoundG        = varLowerBoundG_pc;
//    funcs.reorderParametersG    = reorderParametersG_pc;
//    funcs.outputResultsG        = outputResultsG_pc;
//    setGFunctions( funcs );
//    isGlobalAnalysis = 1;
//}


void outputResults_pc( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    outputPcResults( xn, gv, iv, logFP );
}

//void outputResultsG_pc( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//    outputPcResultsG( xns, gv, ivs, logFP );
//}


void *newModelParameters_pc( xn, sNo )
xnDataSet *xn;
int sNo;
{
    int i;
    pcParameters *p = (void*)malloc( sizeof(pcParameters) );
    
    p->uPiArr = (double*)malloc( sNo * sizeof(double) );
    p->sumUPi = 0.0;
    p->uAMat = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->uAMat[i] = (double*)malloc( sNo * sizeof(double) );
    }
    p->sumUAArr = (double*)malloc( sNo * sizeof(double) );

    p->avgPi = (double *)malloc( sNo * sizeof(double) );
    p->avgLnPi = (double *)malloc( sNo * sizeof(double) );
    p->avgA = (double **)malloc( sNo * sizeof(double*) );
    p->avgLnA = (double **)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ ){
        p->avgA[i] = (double *)malloc( sNo * sizeof(double) );
        p->avgLnA[i] = (double *)malloc( sNo * sizeof(double) );
    }

    p->avgI = (double *)malloc( sNo * sizeof(double) );
    p->avgLnI = (double *)malloc( sNo * sizeof(double) );

    p->aIArr = (double*)malloc( sNo * sizeof(double) );
    p->bIArr = (double*)malloc( sNo * sizeof(double) );

    return p;
}

void freeModelParameters_pc( p, xn, sNo )
void **p;
xnDataSet *xn;
int sNo;
{
    pcParameters *gp = *p;
    int i;
    
    free( gp->uPiArr );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->uAMat[i] );
    }
    free( gp->uAMat );
    free( gp->sumUAArr );
    
    free( gp->avgPi );
    free( gp->avgLnPi );
    for( i = 0 ; i < sNo ; i++ ){
        free( gp->avgA[i] );
        free( gp->avgLnA[i] );
    }
    free( gp->avgA );
    free( gp->avgLnA );

    free( gp->avgI );
    free( gp->avgLnI );
    
    free( gp->aIArr );
    free( gp->bIArr );
    
    free( *p );
    *p = NULL;
}


void *newModelStats_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        pcStats *s = (pcStats*)malloc( sizeof(pcStats) );
        
        int i;
        s->Ni = (double *)malloc( sNo * sizeof(double) );
        s->Nij = (double **)malloc( sNo * sizeof(double*) );
        for( i = 0 ; i < sNo ; i++ )
        {   s->Nij[i] = (double *)malloc( sNo * sizeof(double) );   }
        s->Ci = (double *)malloc( sNo * sizeof(double) );
        s->Mi = (double *)malloc( sNo * sizeof(double) );
        
        return s;
        
//    } else {
//        
//        return NULL;
//        
//    }
}

void freeModelStats_pc( s, xn, gv, iv )
void **s;
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
//    if( isGlobalAnalysis == 0 ){
        int sNo = gv->sNo;
        pcStats *gs = *s;
        int i;
        free( gs->Ni );
        for( i = 0 ; i < sNo ; i++ )
        {   free( gs->Nij[i] );   }
        free( gs->Nij );
        free( gs->Ci );
        free( gs->Mi );
    
        free( gs );
        *s = NULL;
//    }
}

//void *newModelStatsG_pc( xns, gv, ivs)
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    pcGlobalStats *gs = (pcGlobalStats*)malloc( sizeof(pcGlobalStats) );
//    
//    return gs;
//}

//void freeModelStatsG_pc( gs, xns, gv, ivs )
//void **gs;
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//    int sNo = gv->sNo;
//    pcGlobalStats *ggs = *gs;
//    
//    free( *gs );
//    *gs = NULL;
//}


void initializeVbHmm_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcData *d = xn->data;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    pcParameters *p = gv->params;

    int i, j;
    size_t totalC = 0;
    for( i = 0 ; i < dLen ; i++ ){
        totalC += d->counts[i];
    }
    double meanI = (double)totalC / (double)dLen;

    // hyperparameter for p( pi(k) )
    p->sumUPi = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        p->uPiArr[i] = 1.0;
        p->sumUPi += p->uPiArr[i];
    }

    // hyperparameter for p( A(i,j) )
    for( i = 0 ; i < sNo ; i++ ){
        p->sumUAArr[i] = 0.0;
        for( j = 0 ; j < sNo ; j++ ){
            if( j == i ){
                p->uAMat[i][j] = 100.0;
            } else {
                p->uAMat[i][j] = 1.0;
            }
            p->sumUAArr[i] += p->uAMat[i][j];
        }
    }

    // hyperparameter for p( I(k) )
    for( i = 0 ; i < sNo ; i++ ){
        p->aIArr[i] = 1.0;
        p->bIArr[i] = 1.0 / meanI;
    }
    
    initialize_indVars_pc( xn, gv, iv );
    
    calcStatsVars_pc( xn, gv, iv );
    maximization_pc( xn, gv, iv );
}

//void initializeVbHmmG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void initialize_indVars_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    
    int i;
    size_t n;
    double sumPar;
    for( n = 0 ; n < dLen ; n++ ){
        sumPar = 0.0;
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] = enoise(1.0) + 1.0;
            sumPar += gmMat[n][i];
        }
        for( i = 0 ; i < sNo ; i++ ){
            gmMat[n][i] /= sumPar;
        }
    }
}


xnDataSet *newXnDataSet_pc( filename )
const char *filename;
{
    xnDataSet *xn = (xnDataSet*)malloc( sizeof(xnDataSet) );
    xn->name = (char*)malloc( strlen(filename) + 2 );
    strncpy( xn->name, filename, strlen(filename)+1 );
    xn->data = (pcData*)malloc( sizeof(pcData) );
    pcData *d = (pcData*)xn->data;
    d->binSize = 0.0;
    d->counts = NULL;
    return xn;
}

void freeXnDataSet_pc( xn )
xnDataSet **xn;
{
    pcData *d = (pcData*)(*xn)->data;
    free( d->counts );
    free( (*xn)->data );
    free( (*xn)->name );
    free( *xn );
    *xn = NULL;
}


double pTilde_z1_pc( i, params )
int i;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    return exp( p->avgLnPi[i] );
}

double pTilde_zn_zn1_pc( i, j, params )
int i, j;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    return exp( p->avgLnA[i][j] );
}

double pTilde_xn_zn_pc( xnWv, n, i, params )
xnDataSet *xnWv;
size_t n;
int i;
void *params;
{
    pcParameters *p = (pcParameters*)params;
    pcData *xn = (pcData*)xnWv->data;
    return exp( p->avgLnI[i] * (double)xn->counts[n] - p->avgI[i] ) / gsl_sf_fact( xn->counts[n] );
}


void calcStatsVars_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcData *d = (pcData*)xn->data;
    pcStats *s = (pcStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *Ni = s->Ni, **Nij = s->Nij;
    double *Mi = s->Mi, *Ci = s->Ci;
    size_t n;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        Ni[i] = 1e-10;
        Ci[i] = 1e-10;
        Mi[i] = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            Nij[i][j] = 1e-10;
        }
    }
    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){
            Ni[i]  += gmMat[n][i];
            Ci[i]  += gmMat[n][i] * (double)d->counts[n];
            for( j = 0 ; j < sNo ; j++ ){
                Mi[i]     += xiMat[n][i][j];
                Nij[i][j] += xiMat[n][i][j];
            }
        }
    }
}

//void calcStatsVarsG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void maximization_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcParameters *p = (pcParameters*)gv->params;
    pcStats *s = (pcStats*)iv->stats;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *Ni = s->Ni, **Nij = s->Nij;
    double *Mi = s->Mi, *Ci = s->Ci;
    int i, j;

    for( i = 0 ; i < sNo ; i++ ){
        avgPi[i] = ( uPiArr[i] + gmMat[0][i] ) / ( sumUPi + 1.0 );
        avgLnPi[i] = gsl_sf_psi( uPiArr[i] + gmMat[0][i] ) - gsl_sf_psi( sumUPi + 1.0 );
        
        for( j = 0 ; j < sNo ; j++ ){
            avgA[i][j] = ( uAMat[i][j] + Nij[i][j] ) / ( sumUAArr[i] + Mi[i] );
            avgLnA[i][j] = gsl_sf_psi( uAMat[i][j] + Nij[i][j] ) - gsl_sf_psi( sumUAArr[i] + Mi[i] );
        }
        
        avgI[i] = (Ci[i] + aIArr[i]) / (Ni[i] + bIArr[i]);
        avgLnI[i] = gsl_sf_psi( Ci[i] + aIArr[i] ) - log( Ni[i] + bIArr[i] );
#ifdef  INTENSITY_CAP
        size_t n, totalC = 0;
        pcData *pc = xnWv->data;
        for( n = 0 ; n < xnWv->N ; n++ ){
            totalC += pc->counts[n];
        }
        double meanI = (double)totalC / (double)xnWv->N;
        avgI[i] = MIN( avgI[i], maxIntensityRatio * meanI );
        avgLnI[i] = MIN( avgLnI[i], log(maxIntensityRatio * meanI) );
#endif
    }
}

//void maximizationG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


double varLowerBound_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcParameters *p = (pcParameters*)gv->params;
    pcStats *s = (pcStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, *cn = iv->cn;
    double *avgLnPi = p->avgLnPi, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *uPiArr = p->uPiArr, sumUPi = p->sumUPi, *aIArr = p->aIArr, *bIArr = p->bIArr;
    double **uAMat = p->uAMat, *sumUAArr = p->sumUAArr;
    double *Ni = s->Ni, **Nij = s->Nij;
    double *Mi = s->Mi, *Ci = s->Ci;
    size_t n;
    int i, j;
    
    double lnpPi = gsl_sf_lngamma(sumUPi);
    double lnpA = 0.0;
    double lnpI = 0.0;
    double lnqPi = gsl_sf_lngamma(sumUPi + 1.0);
    double lnqA = 0.0;
    double lnqI = 0.0;
    for( i = 0 ; i < sNo ; i++ ){
        lnpPi += (uPiArr[i]-1.0) * avgLnPi[i] - gsl_sf_lngamma(uPiArr[i]);

        lnpI += aIArr[i] * log(bIArr[i]) - gsl_sf_lngamma(aIArr[i]);
        lnpI += (aIArr[i] - 1.0) * avgLnI[i] - bIArr[i] * avgI[i];

        lnqPi += (uPiArr[i]+gmMat[0][i]-1.0) * (gsl_sf_psi(uPiArr[i]+gmMat[0][i]) - gsl_sf_psi(sumUPi+1.0));
        lnqPi -= gsl_sf_lngamma(uPiArr[i] + gmMat[0][i]);

        lnqI += (Ci[i] + aIArr[i]) * log(Ni[i] + bIArr[i]) - gsl_sf_lngamma(Ci[i] + aIArr[i]);
        lnqI += (Ci[i] + aIArr[i] - 1.0) * avgLnI[i] - (Ni[i] + bIArr[i]) * avgI[i];

        lnpA += gsl_sf_lngamma(sumUAArr[i]);
        lnqA += gsl_sf_lngamma(sumUAArr[i] + Mi[i]);
        for( j = 0 ; j < sNo ; j++ ){
            lnpA += (uAMat[i][j]-1.0) * avgLnA[i][j] - gsl_sf_lngamma(uAMat[i][j]);

            lnqA += (uAMat[i][j] + Nij[i][j] - 1.0) * (gsl_sf_psi(uAMat[i][j]+Nij[i][j]) - gsl_sf_psi( sumUAArr[i]+Mi[i]) );
            lnqA -= gsl_sf_lngamma( uAMat[i][j] + Nij[i][j] );
        }
    }

    double lnpX = 0.0;
    for( n = 0 ; n < dLen ; n++ ){
        lnpX += log( cn[n] );
    }
        
    double val;
    val  = lnpPi + lnpA + lnpI;
    val -= lnqPi + lnqA + lnqI;
    val += lnpX;
    val += log(gsl_sf_fact(sNo));

    return val;
}    

//double varLowerBoundG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void reorderParameters_pc( xn, gv, iv )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
{
    pcParameters *p = (pcParameters*)gv->params;
    pcStats *s = (pcStats*)iv->stats;
    size_t dLen = xn->N;
    int sNo = gv->sNo;
    double **gmMat = iv->gmMat, ***xiMat = iv->xiMat;
    double *avgPi = p->avgPi, *avgLnPi = p->avgLnPi, **avgA = p->avgA, **avgLnA = p->avgLnA;
    double *avgI = p->avgI, *avgLnI = p->avgLnI;
    double *Ni = s->Ni, *Ci = s->Ci;
    size_t n;
    int i, j;

    int *index = (int*)malloc( sNo * sizeof(int) );
    double *store = (double*)malloc( sNo * sizeof(double) );
    double **s2D = (double**)malloc( sNo * sizeof(double*) );
    for( i = 0 ; i < sNo ; i++ )
    {   s2D[i] = (double*)malloc( sNo * sizeof(double) );   }

    // index indicates order of avgI values (0=biggest avgI -- sNo=smallest avgI).
    for( i = 0 ; i < sNo ; i++ ){
        index[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgI[i] < avgI[j] ){
                    index[i]--;
                } else if( avgI[i] == avgI[j] ){
                    if( i < j )
                    {   index[i]--;   }
                }
            }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnPi[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnPi[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgI[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgI[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = avgLnI[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnI[i] = store[i];   }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgA[i][j] = s2D[i][j];   }
    }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = avgLnA[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgLnA[i][j] = s2D[i][j];   }
    }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ni[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ni[i] = store[i];   }

    for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = Ci[i];   }
    for( i = 0 ; i < sNo ; i++ ){   Ci[i] = store[i];   }

    for( n = 0 ; n < dLen ; n++ ){
        for( i = 0 ; i < sNo ; i++ ){   store[index[i]] = gmMat[n][i];   }
        for( i = 0 ; i < sNo ; i++ ){   gmMat[n][i] = store[i];   }
    }

    for( n = 0 ; n < dLen ; n++ ){
        for( j = 0 ; j < sNo ; j++ ){
            for( i = 0 ; i < sNo ; i++ ){   s2D[index[i]][index[j]] = xiMat[n][i][j];   }
        }
        for( j = 0 ; j < sNo ; j++ ){
            for( i = 0 ; i < sNo ; i++ ){   xiMat[n][i][j] = s2D[i][j];   }
        }
    }

    for( i = 0 ; i < sNo ; i++ ){   free( s2D[i] );   }
    free( s2D );
    free( store );
    free( index );
}

//void reorderParametersG_pc( xns, gv, ivs )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//{
//}


void outputPcResults( xn, gv, iv, logFP )
xnDataSet *xn;
globalVars *gv;
indVars *iv;
FILE *logFP;
{
    pcParameters *p = (pcParameters*)gv->params;
    int sNo = gv->sNo;

    int i, j;
    fprintf(logFP, "  results: K = %d \n", sNo);
    
    fprintf(logFP, "   intensity: ( %g", p->avgI[0]);
    for( i = 1 ; i < sNo ; i++ )
    {   fprintf(logFP, ", %g", p->avgI[i]);   }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   pi: ( %g", p->avgPi[0]);
    for( i = 1 ; i < sNo ; i++ ){
        fprintf(logFP, ", %g", p->avgPi[i]);
    }
    fprintf(logFP, " ) \n");
    
    fprintf(logFP, "   A_matrix: [");
    for( i = 0 ; i < sNo ; i++ ){
        fprintf(logFP, " ( %g", p->avgA[i][0]);
        for( j = 1 ; j < sNo ; j++ )
        {   fprintf(logFP, ", %g", p->avgA[i][j]);   }
        fprintf(logFP, ")");
    }
    fprintf(logFP, " ] \n\n");

    char fn[256];
    FILE *fp;
    size_t n;

    sprintf( fn, "%s.param%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        fprintf(fp, "I, pi");
        for( i = 0 ; i < sNo ; i++ )
        {   fprintf(fp, ", A%dx", i);   }
        fprintf(fp, "\n");
        
        for( i = 0 ; i < sNo ; i++ ){
            fprintf(fp, "%g, %g", p->avgI[i], p->avgPi[i]);
            for( j = 0 ; j < sNo ; j++ )
            {   fprintf(fp, ", %g", p->avgA[j][i]);   }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    
    sprintf( fn, "%s.Lq%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < gv->iteration ; n++ ){
            fprintf( fp, "%24.20e\n", gv->LqArr[n] );
        }
        fclose(fp);
    }
    
    sprintf( fn, "%s.maxS%03d", xn->name, sNo );
    if( (fp = fopen( fn, "w")) != NULL ){
        for( n = 0 ; n < xn->N ; n++ ){
            fprintf( fp, "%d\n", iv->stateTraj[n] );
        }
        fclose(fp);
    }

}

//void outputPcResultsG( xns, gv, ivs, logFP )
//xnDataBundle *xns;
//globalVars *gv;
//indVarBundle *ivs;
//FILE *logFP;
//{
//}

//
