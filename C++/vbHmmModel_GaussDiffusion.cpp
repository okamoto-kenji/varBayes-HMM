/*
 *  vbHmmModel_GaussDiffusion.c
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


//////////  conditions
vbHmmCond_GaussDiff::vbHmmCond_GaussDiff() : vbHmmCond() {
}

vbHmmCond_GaussDiff::vbHmmCond_GaussDiff(const vbHmmCond_GaussDiff &other) : vbHmmCond( other ) {  // do not copy
}

vbHmmCond_GaussDiff::~vbHmmCond_GaussDiff() {
}


//////////  model
vbHmmModel_GaussDiff::vbHmmModel_GaussDiff() : vbHmmModel() {
    avgDlt = NULL;
    avgLnDlt = NULL;
    uAArr = NULL;
    uBArr = NULL;
    aDlt = NULL;
    bDlt = NULL;
    
    NiR = NULL;
    RiR = NULL;
}

vbHmmModel_GaussDiff::vbHmmModel_GaussDiff(const vbHmmModel_GaussDiff &other) : vbHmmModel(other) {  // do not copy
    avgDlt = NULL;
    avgLnDlt = NULL;
    uAArr = NULL;
    uBArr = NULL;
    aDlt = NULL;
    bDlt = NULL;
    
    NiR = NULL;
    RiR = NULL;
}

vbHmmModel_GaussDiff::~vbHmmModel_GaussDiff(){
    delete avgDlt;
    delete avgLnDlt;
    delete uAArr;
    delete uBArr;
    delete aDlt;
    delete bDlt;

    delete NiR;
    delete RiR;
}

vbHmmModel *vbHmmModel_GaussDiff::newInstance(){
    return new vbHmmModel_GaussDiff();
}

void vbHmmModel_GaussDiff::setSNo( int _sNo ){
    vbHmmModel::setSNo( _sNo );

    avgDlt = new Vec1D( sNo );
    avgLnDlt = new Vec1D( sNo );
    uAArr = new Vec1D( sNo );
    uBArr = new Vec1D( sNo );
    aDlt = new Vec1D( sNo );
    bDlt = new Vec1D( sNo );

    NiR = new Vec1D( sNo );
    RiR = new Vec1D( sNo );
}


void vbHmmModel_GaussDiff::initialize_modelParams( vbHmmCond *c, vbHmmData *d ){
    // hyper parameter for p( delta(k) )
    int i;
    for( i = 0 ; i < sNo ; i++ ){
        uAArr->p[i] = 1.0;
        uBArr->p[i] = 0.0005;
    }
}

double vbHmmModel_GaussDiff::pTilde_xn_zn( vbHmmTrace *_t, int n, int i ){
    double val;
    vbHmmTrace_GaussDiff *t = (vbHmmTrace_GaussDiff*)_t;
    val  = avgLnDlt->p[i] - log(2.0);
    val -= log( t->xn->p[n] ) + avgDlt->p[i] * pow( t->xn->p[n], 2.0) / 4.0;
    return exp(val);
}

void vbHmmModel_GaussDiff::calcStatsVars( vbHmmData *d ){
    int i, j, r, n;
    double *NiRp = NiR->p, *RiRp = RiR->p, *NiiRp = NiiR->p, **NijRp = NijR->p, *z1iRp = z1iR->p;
    double **gmMat, ***xiMat;
    vbHmmTrace_GaussDiff *t;
    for( i = 0 ; i < sNo ; i++ ){
        NiRp[i]   = 1e-10;
        RiRp[i]   = 1e-10;
        NiiRp[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijRp[i][j] = 1e-10;
        }  // j
        z1iRp[i] = 1e-10;
    }  // i
    for( r = 0 ; r < R ; r++ ){
        t = (vbHmmTrace_GaussDiff*)d->traces[r];
        gmMat = ivs[r]->gmMat->p;
        xiMat = ivs[r]->xiMat->p;

        for( i = 0 ; i < sNo ; i++ ){
            z1iRp[i] += gmMat[0][i];
            for( n = 0 ; n < t->N ; n++ ){
                NiRp[i] += gmMat[n][i];
                RiRp[i] += gmMat[n][i] * pow( t->xn->p[n], 2.0 );
                for( j = 0 ; j < sNo ; j++ ){
                    NiiRp[i]    += xiMat[n][i][j];
                    NijRp[i][j] += xiMat[n][i][j];
                }  // j
            }  // n
        }  // i
    }  // r
}

int vbHmmModel_GaussDiff::maximization_modelVars( int s ){
    aDlt->p[s] = uAArr->p[s] + NiR->p[s];
    bDlt->p[s] = uBArr->p[s] + RiR->p[s] / 4.0;

    avgDlt->p[s]   = aDlt->p[s] / bDlt->p[s];
    try{
        avgLnDlt->p[s] = psi( aDlt->p[s] ) - log( bDlt->p[s] );
    }catch (int status){
        return status;
   }
    return 0;
}


double vbHmmModel_GaussDiff::varLowerBound_modelTerm( vbHmmData *d ){
    int s;
    double lnp = 0.0;
    double lnq = - sNo / 2.0;
    for( s = 0 ; s < sNo ; s++ ){
        try{
            lnp += - lnGamma(uAArr->p[s]) + uAArr->p[s] * log(uBArr->p[s]);
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }
        lnp += (uAArr->p[s] - 1.0) * avgLnDlt->p[s] - uBArr->p[s] * avgDlt->p[s];

        try{
            lnq  = - lnGamma(aDlt->p[s]) + aDlt->p[s] * log(bDlt->p[s]);
        }catch (int status){
            return numeric_limits<float>::quiet_NaN();
        }
        lnq += (aDlt->p[s] - 1.0) * avgLnDlt->p[s] - aDlt->p[s];
    }
    return lnp - lnq;
}


void vbHmmModel_GaussDiff::reorderParameters( vbHmmData *d ){
    int r, n;
    int i, j;
    
    Vec1I index( sNo );
    Vec1D store( sNo );
    Mat2D s2D( sNo, sNo );
    
    // index indicates order of avgDlt values (0=biggest avgDlt -- sNo=smallest avgDlt).
    for( i = 0 ; i < sNo ; i++ ){
        index.p[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( avgDlt->p[i] < avgDlt->p[j] ){
                    index.p[i]--;
                } else if( avgDlt->p[i] == avgDlt->p[j] ){
                    if( j > i )
                    {   index.p[i]--;   }
                }
            }
        }
    }
    
    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = avgPi->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgPi->p[i] = store.p[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = avgLnPi->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnPi->p[i] = store.p[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = avgDlt->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgDlt->p[i] = store.p[i];   }
    
    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = avgLnDlt->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   avgLnDlt->p[i] = store.p[i];   }

    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D.p[index.p[i]][index.p[j]] = avgA->p[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgA->p[i][j] = s2D.p[i][j];   }
    }
    
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   s2D.p[index.p[i]][index.p[j]] = avgLnA->p[i][j];   }
    }
    for( j = 0 ; j < sNo ; j++ ){
        for( i = 0 ; i < sNo ; i++ ){   avgLnA->p[i][j] = s2D.p[i][j];   }
    }

    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = NiR->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   NiR->p[i] = store.p[i];   }

    for( r = 0 ; r < R ; r++ ){
        for( n = 0 ; n < ivs[r]->N ; n++ ){
            for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = ivs[r]->gmMat->p[n][i];   }
            for( i = 0 ; i < sNo ; i++ ){   ivs[r]->gmMat->p[n][i] = store.p[i];   }
        }
    }
    
    for( r = 0 ; r < R ; r++ ){
        for( n = 0 ; n < ivs[r]->N ; n++ ){
            for( j = 0 ; j < sNo ; j++ ){
                for( i = 0 ; i < sNo ; i++ ){   s2D.p[index.p[i]][index.p[j]] = ivs[r]->xiMat->p[n][i][j];   }
            }
            for( j = 0 ; j < sNo ; j++ ){
                for( i = 0 ; i < sNo ; i++ ){   ivs[r]->xiMat->p[n][i][j] = s2D.p[i][j];   }
            }
        }
    }
    
}


void vbHmmModel_GaussDiff::outputResults( vbHmmData *d, fstream *logFS ){
    int i, j;
    char str[1024];

    *logFS << "  results: K = " << sNo << endl;
    
    *logFS << "   delta: ( " << avgDlt->p[0];
    for( i = 1 ; i < sNo ; i++ ){
        sprintf(str, ", %g", avgDlt->p[i]);
        *logFS << str;
    }
    *logFS << " )" << endl;
    
    *logFS << "   pi: ( " << avgPi->p[0];
    for( i = 1 ; i < sNo ; i++ ){
        sprintf(str, ", %g", avgPi->p[i]);
        *logFS << str;
    }
    *logFS << " )" << endl;
    
    *logFS << "   A_matrix: [";
    for( i = 0 ; i < sNo ; i++ ){
        sprintf(str, " ( %g", avgA->p[i][0]);
        *logFS << str;
        for( j = 1 ; j < sNo ; j++ ){
            sprintf(str, ", %g", avgA->p[i][j]);
            *logFS << str;
        }
        *logFS << ")";
    }
    *logFS << " ]" << endl << endl;
    
    char fn[1024];
    fstream fs;
    int n;
    
    sprintf( fn, "%s.param%03d", d->name.c_str(), sNo );
    fs.open( fn, ios::out );
    if( fs.is_open() ){
        fs << "delta, pi";
        for( i = 0 ; i < sNo ; i++ ){
           sprintf(str, ", A%dx", i);
            fs << str;
        }
        fs << endl;
        
        for( i = 0 ; i < sNo ; i++ ){
            sprintf(str, "%g, %g", avgDlt->p[i], avgPi->p[i]);
            fs << str;
            for( j = 0 ; j < sNo ; j++ ){
               sprintf(str, ", %g", avgA->p[j][i]);
                fs << str;
            }
            fs << endl;
        }
        fs.close();
    }
    
    sprintf( fn, "%s.Lq%03d", d->name.c_str(), sNo );
    fs.open( fn, ios::out );
    if( fs.is_open() ){
        for( n = 0 ; n < iteration ; n++ ){
            sprintf( str, "%24.20e", LqArr->p[n] );
            fs << str << endl;
        }
        fs.close();
    }
    
    int r;
    sprintf( fn, "%s.maxS%03d", d->name.c_str(), sNo );
    int flag = 0;
    fs.open( fn, ios::out );
    if( fs.is_open() ){
        n = 0;
        do{
            flag = 1;
            for( r = 0 ; r < R ; r++ ){
                if( r > 0 ){
                    fs << ",";
                }
                if( n < d->traces[r]->N ){
                    sprintf( str, "%d", ivs[r]->stateTraj->p[n] );
                    fs << str;
                }
                flag &= (n >= (d->traces[r]->N - 1));
            }
            fs << endl;
            n++;
        }while( !flag );
        fs.close();
    }
}


//
