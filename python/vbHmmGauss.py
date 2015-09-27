#! /usr/local/bin/python3.5
# -*- coding: utf-8 -*-
"""
vbHmmGauss.py
Model-specific core functions for VB-HMM-GAUSS.

Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
Copyright 2011-2015
Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
All rights reserved.

@author: kenji
"""

import math
import numpy as np
import scipy as sp
import scipy.special as spsp
import vbHmm_Common

class vbHmm_Gauss_Data(vbHmm_Common.vbHmm_Data):
    def __init__(self, dLen):
        super().__init__(dLen)
        self.v = np.zeros(self.N, dtype=np.float_)

class vbHmm_Gauss_Params(vbHmm_Common.vbHmm_Params):
    def __init__(self, sNo, dLen):
        super().__init__(sNo, dLen)

        self.uPiArr = np.zeros(sNo, dtype=np.float_)
        self.sumUPi = 0.0;
        self.uAMat = np.zeros((sNo, sNo), dtype=np.float_)
        self.sumUAArr = np.zeros(sNo, dtype=np.float_)

        self.avgPi = np.zeros(sNo, dtype=np.float_)
        self.avgLnPi = np.zeros(sNo, dtype=np.float_)
        self.avgA = np.zeros((sNo, sNo), dtype=np.float_)
        self.avgLnA = np.zeros((sNo, sNo), dtype=np.float_)
        self.avgMu = np.zeros(sNo, dtype=np.float_)
        self.avgLm = np.zeros(sNo, dtype=np.float_)
        self.avgLnLm = np.zeros(sNo, dtype=np.float_)

        self.uBtArr = np.zeros(sNo, dtype=np.float_)
        self.uMuArr = np.zeros(sNo, dtype=np.float_)
        self.uAArr = np.zeros(sNo, dtype=np.float_)
        self.uBArr = np.zeros(sNo, dtype=np.float_)
        self.btMu = np.zeros(sNo, dtype=np.float_)
        self.aLm = np.zeros(sNo, dtype=np.float_)
        self.bLm = np.zeros(sNo, dtype=np.float_)
        self.mu0 = np.zeros(sNo, dtype=np.float_)
        
        self.Ni = np.zeros(sNo, dtype=np.float_)
        self.Nij = np.zeros((sNo, sNo), dtype=np.float_)
        self.Nii = np.zeros(sNo, dtype=np.float_)
        self.barX = np.zeros(sNo, dtype=np.float_)
        self.NiSi = np.zeros(sNo, dtype=np.float_)

class vbHmm_Gauss:
    def __init__(self, data, sNo):
        dLen = data.N

        self.d = data
        self.p = vbHmm_Gauss_Params(sNo, dLen)      # params
        self.r = vbHmm_Common.vbHmm_Results()       # results

#        totalX = self.d.v.sum()
        meanX = self.d.v.mean()
        varX = self.d.v.var()
        precX = 1.0 / varX
#        sys.stderr.write(" mean:{0},  var:{1},  prec{2}\n".format(meanX, varX, precX))

#    // hyper parameter for p( pi(i) )
        for i in range(len(self.p.uPiArr)):
            self.p.uPiArr[i] = 1.0
        self.p.sumUPi = self.p.uPiArr.sum()

#    // hyper parameter for p( A(i,j) )
        for i in range(sNo):
            self.p.sumUAArr[i] = 0.0
            for j in range(sNo):
                if j == i:
                    self.p.uAMat[i][j] = 5.0
                else:
                    self.p.uAMat[i][j] = 1.0
                self.p.sumUAArr[i] += self.p.uAMat[i][j]
    
#    // hyper parameter for p( mu(k), lm(k) )
        for i in range(sNo):
            self.p.uBtArr[i] = precX / 250.0
            self.p.uMuArr[i] = meanX
            self.p.uAArr[i] = 2.5
            self.p.uBArr[i] = 0.01

        for n in range(dLen):
            sumPar = 0.0
            for i in range(sNo):
                self.p.gmMat[n][i] = 2.0 * np.random.random_sample()
                sumPar += self.p.gmMat[n][i]
            for i in range(sNo):
                self.p.gmMat[n][i] /= sumPar

        self.calcStatsVars()
        self.maximization()


    def pTilde_z1(self, i):
        return math.exp(self.p.avgLnPi[i])

    def pTilde_zn_zn1(self, i, j):
        return math.exp(self.p.avgLnA[i][j])

    def pTilde_xn_zn(self, n, i):
        val  = self.p.avgLnLm[i] - math.log(2.0 * np.pi)
        val -= 1.0/self.p.btMu[i] + self.p.aLm[i] / self.p.bLm[i] * (self.d.v[n] - self.p.mu0[i])**2.0
        return math.exp(val / 2.0)


    def calcStatsVars(self):
        sNo = self.p.sNo
        dLen = self.p.dLen

        for i in range(sNo):
            self.p.Ni[i]   = 1.0e-10
            self.p.Nii[i]  = 1.0e-10
            self.p.barX[i] = 1.0e-10
            self.p.NiSi[i] = 1.0e-10
            for j in range(sNo):
                self.p.Nij[i][j] = 1.0e-10

            for n in range(dLen):
                self.p.Ni[i]   += self.p.gmMat[n][i]
                self.p.barX[i] += self.p.gmMat[n][i] * self.d.v[n]
                for j in range(sNo):
                    self.p.Nii[i]    += self.p.xiMat[n][i][j]
                    self.p.Nij[i][j] += self.p.xiMat[n][i][j]
            self.p.barX[i] /= self.p.Ni[i]
            for n in range(dLen):
                self.p.NiSi[i] += self.p.gmMat[n][i] * (self.d.v[n] - self.p.barX[i])**2.0

    def maximization(self):
        sNo = self.p.sNo

        for i in range(sNo):
            self.p.avgPi[i] = (self.p.uPiArr[i] + self.p.gmMat[0][i]) / (self.p.sumUPi + 1.0)
            self.p.avgLnPi[i]  = spsp.psi(self.p.uPiArr[i] + self.p.gmMat[0][i])
            self.p.avgLnPi[i] -= spsp.psi(self.p.sumUPi + 1.0)

            for j in range(sNo):
                self.p.avgA[i][j] = (self.p.uAMat[i][j] + self.p.Nij[i][j]) / (self.p.sumUAArr[i] + self.p.Nii[i])
                self.p.avgLnA[i][j]  = spsp.psi(self.p.uAMat[i][j] + self.p.Nij[i][j]) 
                self.p.avgLnA[i][j] -= spsp.psi(self.p.sumUAArr[i] + self.p.Nii[i])

            self.p.btMu[i] = self.p.uBtArr[i] + self.p.Ni[i]
            self.p.mu0[i]  = (self.p.uBtArr[i] * self.p.uMuArr[i] + self.p.Ni[i] * self.p.barX[i]) / self.p.btMu[i]
            self.p.aLm[i]  = self.p.uAArr[i] + (self.p.Ni[i] + 1.0) / 2.0
            self.p.bLm[i]  = self.p.uBArr[i] + (self.p.NiSi[i] / 2.0)
            self.p.bLm[i] += self.p.uBtArr[i] * self.p.Ni[i] * (self.p.barX[i] - self.p.uMuArr[i])**2.0 / 2.0 / (self.p.uBtArr[i] + self.p.Ni[i])

            self.p.avgMu[i]    = self.p.mu0[i]
            self.p.avgLm[i]    = self.p.aLm[i] / self.p.bLm[i]
            self.p.avgLnLm[i]  = spsp.psi(self.p.aLm[i]) - math.log(self.p.bLm[i])


    def varLowerBound(self):
        sNo = self.p.sNo
        dLen = self.p.dLen

        lnpPi = spsp.gammaln(self.p.sumUPi)
        lnpA = 0.0
        lnpMuLm = 0.0
        lnqPi = spsp.gammaln(self.p.sumUPi + 1.0)
        lnqA = 0.0
        lnqMuLm = - sNo / 2.0

        for i in range(sNo):
            lnpPi += (self.p.uPiArr[i] - 1.0) * self.p.avgLnPi[i] - spsp.gammaln(self.p.uPiArr[i])

            lnpMuLm += math.log(self.p.uBtArr[i]) / 2.0
            lnpMuLm -= self.p.uBtArr[i] / 2.0 * (1.0 / self.p.btMu[i] + self.p.aLm[i] * (self.p.mu0[i] - self.p.uMuArr[i])**2.0 / self.p.bLm[i])
            lnpMuLm += - spsp.gammaln(self.p.uAArr[i]) + self.p.uAArr[i] * math.log(self.p.uBArr[i])
            lnpMuLm += (self.p.uAArr[i] - 0.5) * self.p.avgLnLm[i] - self.p.uBArr[i] * self.p.avgLm[i]
            
            lnqPi += (self.p.uPiArr[i] + self.p.gmMat[0][i] - 1.0) * (spsp.psi(self.p.uPiArr[i] + self.p.gmMat[0][i]) - spsp.psi(self.p.sumUPi + 1.0))
            lnqPi -= spsp.gammaln(self.p.uPiArr[i] + self.p.gmMat[0][i])

            lnpA += spsp.gammaln(self.p.sumUAArr[i])

            lnqA += spsp.gammaln(self.p.sumUAArr[i] + self.p.Nii[i])

            for j in range(sNo):
                lnpA += (self.p.uAMat[i][j] - 1.0) * self.p.avgLnA[i][j] - spsp.gammaln(self.p.uAMat[i][j])

                lnqA += (self.p.uAMat[i][j] + self.p.Nij[i][j] - 1.0) * (spsp.psi(self.p.uAMat[i][j] + self.p.Nij[i][j]) - spsp.psi(self.p.sumUAArr[i] + self.p.Nii[i]))
                lnqA -= spsp.gammaln(self.p.uAMat[i][j] + self.p.Nij[i][j])

            lnqMuLm += math.log(self.p.btMu[i]) / 2.0
            lnqMuLm += -spsp.gammaln(self.p.aLm[i]) + self.p.aLm[i] * math.log(self.p.bLm[i])
            lnqMuLm += (self.p.aLm[i] - 0.5) * self.p.avgLnLm[i] - self.p.aLm[i]

        lnpX = 0.0
        for n in range(dLen):
            lnpX += math.log(self.p.cn[n])

        val  = lnpPi + lnpA + lnpMuLm
        val -= lnqPi + lnqA + lnqMuLm
        val += lnpX
        val += spsp.gammaln(sNo)
        return val


    def reorderParameters(self):
        sNo = self.p.sNo

        index = np.zeros(sNo, dtype=np.int)
        store = np.zeros(sNo, dtype=np.float_)
        s2D = np.zeros((sNo, max(sNo,2)), dtype=np.float_)

#        // index indicates order of avgE values (0=biggest avgE -- sNo=smallest avgE).
        for i in range(sNo):
            index[i] = sNo - 1
            for j in range(sNo):
                if j != i:
                    if self.p.avgMu[i] < self.p.avgMu[j]:
                        index[i] -= 1
                    elif self.p.avgMu[i] == self.p.avgMu[j]:
                        if j > i:
                            index[i] -= 1

        for i in range(sNo):
            store[index[i]] = self.p.avgPi[i]
        for i in range(sNo):
            self.p.avgPi[i] = store[i]

        for i in range(sNo):
            store[index[i]] = self.p.avgLnPi[i]
        for i in range(sNo):
            self.p.avgLnPi[i] = store[i]

        for i in range(sNo):
            store[index[i]] = self.p.avgMu[i]
        for i in range(sNo):
            self.p.avgMu[i] = store[i]

        for i in range(sNo):
            store[index[i]] = self.p.avgLm[i]
        for i in range(sNo):
            self.p.avgLm[i] = store[i]

        for i in range(sNo):
            store[index[i]] = self.p.avgLnLm[i]
        for i in range(sNo):
            self.p.avgLnLm[i] = store[i]

        for j in range(sNo):
            for i in range(sNo):
                s2D[index[i]][index[j]] = self.p.avgA[i][j]
        for j in range(sNo):
            for i in range(sNo):
                self.p.avgA[i][j] = s2D[i][j]

        for j in range(sNo):
            for i in range(sNo):
                s2D[index[i]][index[j]] = self.p.avgLnA[i][j]
        for j in range(sNo):
            for i in range(sNo):
                self.p.avgLnA[i][j] = s2D[i][j]

        for i in range(sNo):
            store[index[i]] = self.p.Ni[i]
        for i in range(sNo):
            self.p.Ni[i] = store[i]

        for n in range(sNo):
            for i in range(sNo):
                store[index[i]] = self.p.gmMat[n][i]
            for i in range(sNo):
                self.p.gmMat[n][i] = store[i]

        for n in range(sNo):
            for j in range(sNo):
                for i in range(sNo):
                    s2D[index[i]][index[j]] = self.p.xiMat[n][i][j]
            for j in range(sNo):
                for i in range(sNo):
                    self.p.xiMat[n][i][j] = s2D[i][j]


    def outputResults(self, out_name, logFP):
        sNo = self.p.sNo
        dLen = self.p.dLen

        logFP.write("  results: K = {0} \n".format(sNo))

        logFP.write("   means: ( {0}".format(self.p.avgMu[0]))
        for i in range(1,sNo-1):
            logFP.write(", {0}".format(self.p.avgMu[i]))
        logFP.write(" ) \n")

        logFP.write("   lambda: ( {0}".format(self.p.avgLm[0]))
        for i in range(1,sNo-1):
            logFP.write(", {0}".format(self.p.avgLm[i]))
        logFP.write(" ) \n")

        logFP.write("   A_matrix: [")
        for i in range(sNo):
            logFP.write(" ( {0}".format(self.p.avgA[i][0]))
            for j in range(1,sNo-1):
                logFP.write(", {0}".format(self.p.avgA[i][j]))
            logFP.write(")")
        logFP.write(" ] \n\n")

        fn = "{0}.param{1:03d}".format(out_name, sNo)
        with open(fn, "w") as fp:
            fp.write("mu, lambda")
            for i in range(sNo):
                fp.write(", A{0}x".format(i))
            fp.write("\n")

            for i in range(sNo):
                fp.write("{0}, {1}".format(self.p.avgMu[i], self.p.avgLm[i]))
                for j in range(sNo):
                    fp.write(", {0}".format(self.p.avgA[j][i]))
                fp.write("\n")

        fn = "{0}.Lq{1:03d}".format(out_name, sNo)
        with open(fn, "w") as fp:
            for n in range(self.r.iteration):
                fp.write("{0:24.20f}\n".format(self.r.LqArr[n]))

        fn = "{0}.maxS{1:03d}".format(out_name, sNo)
        with open(fn, "w") as fp:
            for n in range(dLen):
                fp.write("{0:d}\n".format(self.r.maxSumTraj[n]))


import sys
import time

if __name__ == "__main__":

    startTime = time.time()
    
    sFrom = 0
    sTo = 0
    trials = 0
    maxIteration = 0
    threshold = 0.0
    
#    // Get command line arguments.
    if len(sys.argv) != 7:
#        // Invalid number of arguments
        sys.stderr.write("\nVariational-Bayesian Hidden-Markov-Model Gauss Analysis\n")
        sys.stderr.write("    built on 2015.09.xx\n" )
        sys.stderr.write("Syntax : vbHmmGauss sFrom sMax trials maxIteration threshold data_filename\n" )
        sys.stderr.write("         sFrom        : minimum number of state to analyze\n" )
        sys.stderr.write("         sMax         : maximum number of state to analyze\n" )
        sys.stderr.write("         trials       : number of repetition\n" )
        sys.stderr.write("                        (giving chance to change initial conditions)\n" )
        sys.stderr.write("         maxIteration : maximum iteration if inference does not reach threshold\n" )
        sys.stderr.write("         threshold    : stop iteration if inference reaches this value\n" )
        sys.stderr.write("         data_filename   : name of data file\n" )
        sys.stderr.write("Example: vbHmmGauss 2 20 5 1000 1e-10 data-file\n\n" )
        sys.exit(1)

    elif len(sys.argv) == 7:
#        // Set parameters
        sFrom = int(sys.argv[1])
        sTo = int(sys.argv[2])
        trials = int(sys.argv[3])
        maxIteration = int(sys.argv[4])
        threshold = float(sys.argv[5])
        
        in_name = sys.argv[6]
        out_name = sys.argv[6]
        
        logFilename = out_name + ".log"
        try:
            logFP = open(logFilename, "w")
        except (FileNotFoundError, TypeError) as e:
            logFP = sys.stderr


#    // Prepare data, typically by loading from file(s).
    import vbHmmGaussDataHandler as dh
    data = dh.readGaussText(in_name, logFP)

#    // If data can be obtained, execute analysis.
    if data != None:
        vbhmm = vbHmm_Common.vbHmmMain(vbHmm_Gauss)
        optK = vbhmm.modelComparison(data, sFrom, sTo, trials, maxIteration, threshold, out_name, logFP)
        logFP.write(" No. of state = {0:d} was chosen.\n".format(optK))

    logFP.write("FINISH: {0:d} sec. spent.\n".format(int(time.time() - startTime)))
    if logFP != sys.stderr:
        logFP.close()
    

if __name__ == "__main__":
    pass


#########################################################
