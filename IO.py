#!/usr/bin/python
"""
Simple functions for loading or saving a single array as a *.txt file.
"""
################################################################################

import numpy as np

################################################################################

def loadArray(fileName) :
    x = np.array([])
    with open(fileName) as fh :
        for line in fh :
            ell = line.strip()
            x = np.hstack((x, np.float64(ell)))
    return x

################################################################################

def saveArray(fileName, values) :
    with open(fileName, "w") as fh :
        for i in range(len(values)) :
            fh.write('{0:17.14f}\n'.format(values[i]))

