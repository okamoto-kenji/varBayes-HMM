/*
 *  vbHmmModel_template.c
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


//////////  conditions
vbHmmCond_template::vbHmmCond_template() : vbHmmCond() {
}

vbHmmCond_template::vbHmmCond_template(const vbHmmCond_template &other) : vbHmmCond( other ) {  // do not copy
}

vbHmmCond_template::~vbHmmCond_template() {
}


//////////  model
vbHmmModel_template::vbHmmModel_template() : vbHmmModel() {
}

vbHmmModel_template::vbHmmModel_template(const vbHmmModel_template &other) : vbHmmModel(other) {  // do not copy
}

vbHmmModel_template::~vbHmmModel_template(){
}

vbHmmModel *vbHmmModel_template::newInstance(){
    return new vbHmmModel_template();
}

void vbHmmModel_template::setSNo( int _sNo ){
    vbHmmModel::setSNo( _sNo );

    // model-specific process here
}


void vbHmmModel_template::initialize_modelParams( vbHmmCond *c, vbHmmData *d ){
}

double vbHmmModel_template::pTilde_xn_zn( vbHmmTrace *_t, int n, int i ){
    return 1.0;
}

void vbHmmModel_template::calcStatsVars( vbHmmData *d ){
    // rewriting the whole function may be advantagous in preformance to reduce loops

    int i, j, n, r, N;
    
    double *NiiRp = NiiR->p, **NijRp = NijR->p, *z1iRp = z1iR->p;
    double **gmMat, ***xiMat;
    for( i = 0 ; i < sNo ; i++ ){
        NiiRp[i]  = 1e-10;
        for( j = 0 ; j < sNo ; j++ ){
            NijRp[i][j] = 1e-10;
        }
        z1iRp[i] = 1e-10;
    }  // i
    for( r = 0 ; r < R ; r++ ){
        for( i = 0 ; i < sNo ; i++ ){
            N = ivs[r]->N;
            gmMat = ivs[r]->gmMat->p;
            xiMat = ivs[r]->xiMat->p;
            
            z1iRp[i] += gmMat[0][i];
            for( n = 0 ; n < N ; n++ ){
                for( j = 0 ; j < sNo ; j++ ){
                    NiiRp[i]    += xiMat[n][i][j];
                    NijRp[i][j] += xiMat[n][i][j];
                }  // j
            }  // n
        }  // i
    }  // r
}

int vbHmmModel_template::maximization_modelVars( int s ){
    return 0;
}


double vbHmmModel_template::varLowerBound_modelTerm( vbHmmData *d ){
    double lnp = 0.0;
    double lnq = 0.0;
    return lnp - lnq;
}


void vbHmmModel_template::reorderParameters( vbHmmData *d ){
    int r, n;
    int i, j;
    
    Vec1I index( sNo );
    Vec1D store( sNo );
    Mat2D s2D( sNo, sNo );
    
    // index indicates order of NiiR values (0=biggest NiiR -- sNo=smallest NiiR).
    for( i = 0 ; i < sNo ; i++ ){
        index.p[i] = sNo - 1;
        for( j = 0 ; j < sNo ; j++ ){
            if( j != i ){
                if( NiiR->p[i] < NiiR->p[j] ){
                    index.p[i]--;
                } else if( NiiR->p[i] == NiiR->p[j] ){
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

    for( i = 0 ; i < sNo ; i++ ){   store.p[index.p[i]] = NiiR->p[i];   }
    for( i = 0 ; i < sNo ; i++ ){   NiiR->p[i] = store.p[i];   }

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


void vbHmmModel_template::outputResults( vbHmmData *d, fstream *logFS ){
    int i, j;
    char str[1024];

    *logFS << "  results: K = " << sNo << endl;
    
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
        fs << "pi";
        for( i = 0 ; i < sNo ; i++ ){
           sprintf(str, ", A%dx", i);
            fs << str;
        }
        fs << endl;
        
        for( i = 0 ; i < sNo ; i++ ){
            sprintf(str, "%g", avgPi->p[i]);
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
