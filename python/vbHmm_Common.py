#! /usr/local/bin/python3.5
# -*- coding: utf-8 -*-
"""
vbHmm_Common.py
Common VB-HMM engine.

Reference: 
  Christopher M. Bishop, "Pattern Recognition and Machine Learning", Springer, 2006

Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
Copyright 2011-2015
Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
All rights reserved.

@author: kenji
"""

import math
import numpy as np
import sys

## Structure encapsulating data to analyze
class vbHmm_Data:
    def __init__(self, dLen):
        self.N = dLen       # number of data points
        self.T = None        # total length of data (may not be always necessary)

## Commonly used parameters
class vbHmm_Params:
    def __init__(self, sNo, dLen):
        self.dLen = dLen                                        # number of data points
        self.sNo = sNo                                         # number of states
        
        self.gmMat = np.zeros((dLen, sNo), dtype=np.float_)      # gamma matrix for E-step
        self.xiMat = np.zeros((dLen, sNo, sNo), dtype=np.float_)  # xi matrix for E-step
        self.aMat = np.zeros((dLen, sNo), dtype=np.float_)       # alpha matrix for Baum-Welch
        self.bMat = np.zeros((dLen, sNo), dtype=np.float_)       # beta matrix for Baum-Welch
        self.cn = np.zeros(dLen, dtype=np.float_)               # scaling factor for Baum-Welch

        # temporary storage of calculation to save time
        self.valpZnZn1 = np.zeros((sNo, sNo), dtype=np.float_)
        self.valpXnZn = np.zeros((dLen, sNo), dtype=np.float_)

## Results of analysis
class vbHmm_Results:
    def __init__(self):
        self.iteration = 0              # calculation steps
        self.maxLq = 0.0                # final lower bound
        self.LqArr = None               # time series of lower bound
        self.maxSumTraj = None          # resulting state transition trajectory



#//////////////////////////////////////////////////////////////////  VB-HMM Execution Functions
class vbHmmMain:

    def __init__(self, modelClass):
        self.modelClass = modelClass

    def modelComparison(self, data, sFrom, sTo, trials, maxIteration, threshold, out_name, logFP):
        logFP.write("  No. of states from {0:d} to {1:d},  trials = {2:d}, ".format(sFrom, sTo, trials))
        logFP.write("  analyze:  maxIteration = {0:d},  threshold = {1} \n\n".format(maxIteration, threshold))

        LqVsK = np.zeros(trials * (sTo - sFrom + 1), dtype=np.float_)

        maxS = 0
        maxLq = -sys.float_info.max
        for s in range(sFrom, sTo+1):
            modelArray = []
            for t in range(trials):
                st = (s - sFrom) * trials + t
                model = self.modelClass(data, s)
                modelArray.append(model)
                
                LqVsK[st] = self.vbHmm_Main(model, maxIteration, threshold, logFP)

                if LqVsK[st] > maxLq:
                    maxLq = LqVsK[st]
                    maxS = s

            maxLqForS = LqVsK[(s - sFrom) * trials]
            maxT = 0
            for t in range(1, trials):
                st = (s - sFrom) * trials + t
                if LqVsK[st] > maxLqForS:
                    maxLqForS = LqVsK[st]
                    maxT = t

            modelArray[maxT].outputResults(out_name, logFP)

            if s >= (maxS+3):
                sTo = s
                break

        fn = out_name + ".LqVsK"
        with open(fn, "w") as fp:
            for s in range(trials * (sTo - sFrom + 1)):
                fp.write("{0:2d}, {1:.20f}\n".format((s//trials) + sFrom, LqVsK[s]))

        return maxS


    #//////////////////////////////////////////////////////////////////  VB-HMM Common Engine
    def vbHmm_Main(self, m, maxIteration, threshold, logFP):

        m.r.LqArr = np.zeros(maxIteration, dtype=np.float_)

        for i in range(maxIteration):
#            // E-step
            self.forwardBackward(m)

            m.calcStatsVars()
            m.r.LqArr[i] = m.varLowerBound()

#            // End loop if derivative of variational lower bound reaches threshold.
            if i > 0  and (math.fabs((m.r.LqArr[i] - m.r.LqArr[i-1]) / m.r.LqArr[i]) < threshold):
                maxIteration = i + 1
                break

#            // M-step
            m.maximization()

        else:
            logFP.write("MAX iteration ({0:d}) reached.\n".format(maxIteration))

        m.reorderParameters()
        self.maxSum(m)

        m.r.iteration = maxIteration
        m.r.LqArr.resize(maxIteration)
        m.r.maxLq = m.r.LqArr[maxIteration-1]
        logFP.write("  iteration: {0:d}    evidence p(x|K={1:d}) = {2:.20f} \n".format(maxIteration, m.p.sNo, m.r.maxLq))

        return m.r.maxLq


#// Baum-Welch algorithm for E-step calculation
    def forwardBackward(self, m):
        sNo = m.p.sNo
        dLen = m.p.dLen

#        // forward
        m.p.cn[0] = 0.0
        for i in range(sNo):
            m.p.valpXnZn[0][i] = m.pTilde_xn_zn(0, i)
            m.p.aMat[0][i] = m.pTilde_z1(i) * m.p.valpXnZn[0][i]

            m.p.cn[0] += m.p.aMat[0][i]

            for j in range(sNo):
                m.p.valpZnZn1[i][j] = m.pTilde_zn_zn1(i, j)

        for i in range(sNo):
            m.p.aMat[0][i] /= m.p.cn[0]

        for n in range(1, dLen):
            m.p.cn[n] = 0.0
            for j in range(sNo):
                m.p.aMat[n][j] = 0.0
                for i in range(sNo):
                    m.p.aMat[n][j] += m.p.aMat[n-1][i] * m.p.valpZnZn1[i][j]

                m.p.valpXnZn[n][j] = m.pTilde_xn_zn(n, j)
                m.p.aMat[n][j] *= m.p.valpXnZn[n][j]
                m.p.cn[n] += m.p.aMat[n][j]

            for j in range(sNo):
                m.p.aMat[n][j] /= m.p.cn[n]

#        // backward
        for i in range(sNo):
            m.p.bMat[dLen-1][i] = 1.0

        for nn in range(dLen-1):
            n = (dLen - 2) - nn
            for i in range(sNo):
                m.p.bMat[n][i] = 0.0
                for j in range(sNo):
                    m.p.bMat[n][i] += m.p.bMat[n+1][j] *  m.p.valpZnZn1[i][j] *  m.p.valpXnZn[n+1][j]
                m.p.bMat[n][i] /= m.p.cn[n+1]

#        // update gamma
        for n in range(dLen):
            for i in range(sNo):
                m.p.gmMat[n][i] = m.p.aMat[n][i] * m.p.bMat[n][i]

#        // update xi
        for i in range(sNo):
            for j in range(sNo):
                m.p.xiMat[0][i][j] = 0.0

        for n in range(1, dLen):
            for i in range(sNo):
                for j in range(sNo):
                    xiTerm  = m.p.aMat[n-1][i]    * m.p.valpXnZn[n][j]
                    xiTerm *= m.p.valpZnZn1[i][j] * m.p.bMat[n][j]
                    m.p.xiMat[n][i][j] = xiTerm / m.p.cn[n]


#// Viterbi algorithm to construct most likely trajectory
    def maxSum(self, m):
        sNo = m.p.sNo
        dLen = m.p.dLen

        m.r.maxSumTraj = np.zeros(dLen, dtype=np.int_)
        wnMat = np.zeros((dLen, sNo), dtype=np.float_)
        phiMat = np.zeros((dLen, sNo), dtype=np.float_)

#        // forward
        for i in range(sNo):
            wnMat[0][i] = math.log(m.pTilde_z1(i)) + math.log(m.pTilde_xn_zn(0, i))

        for n in range(1, dLen):
            for j in range(sNo):
                maxWn = math.log(m.pTilde_zn_zn1(0, j)) + wnMat[n-1][0]
                maxI = 0
                for i in range(1, sNo):
                    wnTest = math.log(m.pTilde_zn_zn1(i, j)) + wnMat[n-1][i]
                    if wnTest > maxWn:
                        maxWn = wnTest
                        maxI = i

                phiMat[n][j] = maxI
                try:
                    wnMat[n][j] = math.log(m.pTilde_xn_zn(n, j)) + maxWn
                except ValueError:
#                    sys.stderr.write("  -- n:{0:d}, j:{1:d}/{2:d} - {3}\n".format(n, j, sNo, m.pTilde_xn_zn(n, j)))
                    wnMat[n][j] = -sys.float_info.max

#        // backward
        maxWn = wnMat[dLen-1][0]
        maxI = 0
        for i in range(1, sNo):
            if wnMat[dLen-1][i] > maxWn:
                maxWn = wnMat[dLen-1][i]
                maxI = i

        m.r.maxSumTraj[dLen-1] = maxI
        for nn in range(dLen-1):
            n = (dLen - 1) - nn
            m.r.maxSumTraj[n-1] = phiMat[n][m.r.maxSumTraj[n]]


if __name__ == "__main__":
    pass


#########################################################
