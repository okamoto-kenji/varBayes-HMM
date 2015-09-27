#! /usr/local/bin/python3.5
# -*- coding: utf-8 -*-
"""
vbHmmGaussDataHandler.py
File loader for VB-HMM-GAUSS data.

Created by OKAMOTO Kenji, SAKO Yasushi and RIKEN
Copyright 2011-2015
Cellular Informatics Laboratory, Advance Science Institute, RIKEN, Japan.
All rights reserved.

@author: kenji
"""

import numpy as np
import vbHmmGauss

#define DATA_READ_NUMBER 1000

def readGaussText(filename, logFP):
#    floatArray = readDoubleArrayFromFile(filename, logFP)

    with open(filename, "r") as fp:
        logFP.write("  Reading data from '{0}' ... ".format(filename))
            
#    // file exists. start reading in data and correlate
        array = []
        for line in fp:
            array.append(float(line))

    logFP.write("done.\n")

    if array != None:
        dLen = len(array)
        data = vbHmmGauss.vbHmm_Gauss_Data(dLen)
        for i in range(dLen):
            data.v[i] = array[i]
        logFP.write("  Total of {0:d} points read in.\n\n".format(data.N))

    else:
        data = None

    return data


if __name__ == "__main__":
    pass


##########################################################
