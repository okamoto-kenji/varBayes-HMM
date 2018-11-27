/*
 *  vbHmmCore.cpp
 *  Common VB-HMM engine.
 *  Reference: 
 *    Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006
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
#include <fstream>
#include <cfloat>
#include <cmath>
#include "vbHmmCore.h"


//////////////////////////////////////////////////////////////////  VB-HMM Execution Functions
int modelComparison( vbHmmCond *c, vbHmmData *d, vbHmmModel *m, fstream *logFS ){

    *logFS << "  No. of states from " << c->sFrom << " to " << c->sTo << ",  trials = " << c->trials << ", ";
    *logFS << "  analyze:  maxIteration = " << c->maxIteration << ",  threshold = " << c->lbTh << endl << endl;

    int r, s, t;

    vbHmmModel **ms = (vbHmmModel **)malloc( c->trials * sizeof(vbHmmModel*) );
    for( t = 0 ; t < c->trials ; t++ )  ms[t] = NULL;

    Vec1D LqVsK( c->trials * (c->sTo - c->sFrom + 1) );
    int maxS = 0;
    double maxLq = -DBL_MAX;

    for( s = c->sFrom ; s <= c->sTo ; s++ ){
        for( t = 0 ; t < c->trials ; t++ ){
            if( t < (c->trials / 2) ){
                c->initType = initType_kMeans;
            } else {
                c->initType = initType_random;
            }

            int st = (s - c->sFrom) * c->trials + t;
            ms[t] = m->newInstance();
            ms[t]->setSNo( s );
            for( r = 0 ; r < d->R ; r++ ){
                ms[t]->appendIndVars( d->traces[r]->N );
            }

            int status;
                status = vbHmm_Main( c, d, ms[t], logFS );
            if( status == 0 ){
                LqVsK.p[st] = ms[t]->maxLq;
            } else {
                LqVsK.p[st] = numeric_limits<float>::quiet_NaN();
            }
            if( LqVsK.p[st] > maxLq ){
                maxLq = LqVsK.p[st];
                maxS = s;
            }
        }
        
        double maxLqForS = -DBL_MAX;
        vbHmmModel *maxM = NULL;
        for( t = 0 ; t < c->trials ; t++ ){
            int st = (s - c->sFrom) * c->trials + t;
            if( LqVsK.p[st] > maxLqForS ){
                maxLqForS = LqVsK.p[st];
                maxM = ms[t];
            }
        }
        if( maxM != NULL )  maxM->outputResults( d, logFS );
        maxM = NULL;

        for( t = 0 ; t < c->trials ; t++ )  delete ms[t];
        
        if( s >= (maxS+3) ){
            s++;
            break;
        }
        
    }
    c->sTo = s - 1;
    free( ms );

    char str[1024];
    fstream fs;
    fs.open( d->name + ".LqVsK", ios::out );
    if( fs.is_open() ){
        for( s = 0 ; s < (c->trials * (c->sTo - c->sFrom + 1)) ; s++ ){
            sprintf( str, "%2d, %.20g", (s/c->trials) + c->sFrom, LqVsK.p[s] );
            fs << str << endl;
        }
        fs.close();
    }
    return maxS;
}


//////////////////////////////////////////////////////////////////  VB-HMM Common Engine
int vbHmm_Main( vbHmmCond *c, vbHmmData *d, vbHmmModel *m, fstream *logFS ){
    m->LqArr->resize( c->maxIteration );
    double *LqArr = m->LqArr->p;

    m->initialize( c, d );

    int i, r;
    for( i = 0 ; i < c->maxIteration ; i++ ){
        int status;

        // E-step
        forwardBackward( c, d, m );

        m->calcStatsVars( d );
        LqArr[i] = m->varLowerBound( d );
        if( !isfinite( LqArr[i] ) ){
            cerr << "Error caught in varLowerBound()." << endl;
            *logFS << "Error caught in varLowerBound()." << endl;
            return 1;
        }

        // End loop if derivative of variational lower bound reaches threshold.
        if( (i > 0) && ( fabs( (LqArr[i] - LqArr[i-1]) / LqArr[i] ) < c->lbTh ) ){
            break;
        }

        // M-step
        status = m->maximization();
        if( status != 0 ){
            cerr << "Error " << status << " caught in maximization()." << endl;
            *logFS << "Error " << status << " caught in maximization()." << endl;
            return status;
        }
    }
    if( i == c->maxIteration ){
        *logFS << "MAX iteration (" << c->maxIteration << ") reached." << endl;
        i--;
    }
    m->reorderParameters( d );

#ifdef OUTPUT_MAX_GAMMA
    for( r = 0 ; r < m->R ; r++ ){
        maxGamma( c, d, m );
    }
#endif
    for( r = 0 ; r < m->R ; r++ ){
        maxSum( c, d, m );
    }

    m->iteration = i+1;
    m->LqArr->resize( i+1 );
    m->maxLq = LqArr[i];

    char str[1024];
    sprintf( str, "  iteration: %d    evidence p(x|K=%d) = %.20g ", i+1, m->sNo, m->maxLq );
    *logFS << str << endl;

    return 0;
}


void forwardBackward( vbHmmCond *c, vbHmmData *d, vbHmmModel *m ){
    int  r, n, i, j, N;
    int sNo = m->sNo;
    double **gmMat, ***xiMat, **aMat, **bMat, *cn, **valpZnZn1, **valpXnZn;
    for( r = 0 ; r < m->R ; r++ ){
        N = m->ivs[r]->N;
        gmMat = m->ivs[r]->gmMat->p;
        xiMat = m->ivs[r]->xiMat->p;
        aMat = m->ivs[r]->aMat->p;
        bMat = m->ivs[r]->bMat->p;
        cn = m->ivs[r]->cn->p;
        valpZnZn1 = m->ivs[r]->valpZnZn1->p;
        valpXnZn = m->ivs[r]->valpXnZn->p;

        // forward
        cn[0] = 0.0;
        for( i = 0 ; i < sNo ; i++ ){
            valpXnZn[0][i] = m->pTilde_xn_zn( d->traces[r], 0, i );
            aMat[0][i] = m->pTilde_z1( i ) * valpXnZn[0][i];

            cn[0] += aMat[0][i];

            for( j = 0 ; j < sNo ; j++ ){
                valpZnZn1[i][j] = m->pTilde_zn_zn1( i, j );
            }  // j
        }  // i
        for( i = 0 ; i < sNo ; i++ ){
            aMat[0][i] /= cn[0];
        }  // i
        for( n = 1 ; n < N ; n++ ){
            cn[n] = 0.0;
            for( j = 0 ; j < sNo ; j++ ){
                aMat[n][j] = 0;
                for( i = 0 ; i < sNo ; i++ ){
                    aMat[n][j] += aMat[n-1][i] * valpZnZn1[i][j];
                }  // i
                valpXnZn[n][j] = m->pTilde_xn_zn( d->traces[r], n, j );
                aMat[n][j] *= valpXnZn[n][j];
                cn[n] += aMat[n][j];
            }  // j
            for( j = 0 ; j < sNo ; j++ ){
                aMat[n][j] /= cn[n];
            }  // j
        }  // n

        // backward
        for( i = 0 ; i < sNo ; i++ ){
            bMat[N-1][i] = 1.0;
        }  // i
        double betaTerm;
        for( n = N-1 ; n > 0 ; ){
            n--;
            for( i = 0 ; i < sNo ; i++ ){
                bMat[n][i] = 0;
                for( j = 0 ; j < sNo ; j++ ){
                    betaTerm  = bMat[n+1][j];
                    betaTerm *= valpZnZn1[i][j];
                    betaTerm *= valpXnZn[n+1][j];
                    bMat[n][i] += betaTerm;
                }  // j
                bMat[n][i] /= cn[n+1];
            }  // i
        }  // n

        // update gamma
        for( n = 0 ; n < N ; n++ ){
            for( i = 0 ; i < sNo ; i++ ){
                gmMat[n][i] = aMat[n][i] * bMat[n][i];
            }  // i
        }  // n

        // update xi
        double xiTerm;
        for( i = 0 ; i < sNo ; i++ ){
            for( j = 0 ; j < sNo ; j++ ){
                xiMat[0][i][j] = 0;
        }   }  // j, i
        for( n = 1 ; n < N ; n++ ){
            for( i = 0 ; i < sNo ; i++ ){
                for( j = 0 ; j < sNo ; j++ ){
                    xiTerm  = aMat[n-1][i];
                    xiTerm *= valpXnZn[n][j];
                    xiTerm *= valpZnZn1[i][j];
                    xiTerm *= bMat[n][j];
                    xiMat[n][i][j] = xiTerm / cn[n];
                }  // j
            }  // i
        }  // n
    }  // r

}

#ifdef OUTPUT_MAX_GAMMA
//// construct most likely trajectory to trace max Gamma
void maxGamma( vbHmmCond *c, vbHmmData *d, vbHmmModel *m ){
    int  r, n, i, N;
    int sNo = m->sNo;
    double **gmMat;
    Vec1I *gammaTraj;
    for( r = 0 ; r < m->R ; r++ ){
        N = m->ivs[r]->N;
        gmMat = m->ivs[r]->gmMat->p;
    
        gammaTraj = m->ivs[r]->gammaTraj;
        gammaTraj->resize( N );

        int maxI = -1;
        double maxG;

        for( n = 0 ; n < N ; n++ ){
            maxG = - 1.0;
            for( i = 0 ; i < sNo ; i++ ){
                if( gmMat[n][i] > maxG ){
                    maxG = gmMat[n][i];
                    maxI = i;
                }
            }  // i
            gammaTraj->p[n] = maxI;
        }  // n
    }  // r
}
#endif


// Viterbi algorithm to construct most likely trajectory
void maxSum( vbHmmCond *c, vbHmmData *d, vbHmmModel *m ){
    int  r, n, i, j, N;
    int sNo = m->sNo;
    Vec1I *stateTraj;
    int maxI;
    double wnTest, maxWn;
    for( r = 0 ; r < m->R ; r++ ){
        N = d->traces[r]->N;

        stateTraj = m->ivs[r]->stateTraj;
        stateTraj->resize( N );
        Mat2D wnMat( N, sNo ), phiMat( N, sNo );

        // forward
        for( i = 0 ; i < sNo ; i++ ){
            wnMat.p[0][i] = log( m->pTilde_z1( i ) ) + log( m->pTilde_xn_zn( d->traces[r], 0, i ) );
        }  // i
        for( n = 1 ; n < N ; n++ ){
            for( j = 0 ; j < sNo ; j++ ){
                maxWn = log( m->pTilde_zn_zn1(0, j) ) + wnMat.p[n-1][0];
                maxI = 0;
                for( i = 1 ; i < sNo ; i++ ){
                    wnTest = log( m->pTilde_zn_zn1(i, j) ) + wnMat.p[n-1][i];
                    if( wnTest > maxWn ){
                        maxWn = wnTest;
                        maxI = i;
                    }
                }  // i
                phiMat.p[n][j] = maxI;
                wnMat.p[n][j] = log( m->pTilde_xn_zn(d->traces[r], n, j) ) + maxWn;
            }  // j
        }  // n

        // backward
        maxWn = wnMat.p[N-1][0];
        maxI = 0;
        for( i = 1 ; i < sNo ; i++ ){
            if( wnMat.p[N-1][i] > maxWn ){
                maxWn = wnMat.p[N-1][i];
                maxI = i;
            }
        }  // i
        stateTraj->p[N-1] = maxI;
        for( n = N-1 ; n > 0 ; n-- ){
            stateTraj->p[n-1] = phiMat.p[n][stateTraj->p[n]];
        }  // n
    }  // r
}


//
