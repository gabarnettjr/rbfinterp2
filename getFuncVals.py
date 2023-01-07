#!/usr/bin/python
"""
If you already have defined nodes (x, y) and evaluation points (xe, ye), then
this script will use a smooth function to define data at these points.
The script generates f.txt to go with (x, y), and fe.txt to go with (xe, ye).
In order for this to work, you need to already have x.txt, y.txt, xe.txt, and
ye.txt saved in the folder at location $dataDir.

Greg Barnett
January 2023
"""
################################################################################

import os
from sys import path, argv
import numpy as np

path.append(".")
import IO

################################################################################

dataDir = os.path.join("randomCoords", "smoothData")
ftype = "peaks and valleys"

argv = argv[1:]
if len(argv) > 0 :
    dataDir = argv[0];  argv = argv[1:]
    if not os.path.isdir(dataDir) :
        raise ValueError("First input must be a data directory.")

if len(argv) > 0 : ftype = argv[0].lower();  argv = argv[1:]

################################################################################

x  = IO.loadArray(os.path.join(dataDir, "..", "x.txt"))
y  = IO.loadArray(os.path.join(dataDir, "..", "y.txt"))
xe = IO.loadArray(os.path.join(dataDir, "..", "xe.txt"))
ye = IO.loadArray(os.path.join(dataDir, "..", "ye.txt"))

ab = np.hstack((x, xe))
a = np.min(ab)
b = np.max(ab)
cd = np.hstack((y, ye))
c = np.min(cd)
d = np.max(cd)

################################################################################

def eff(ftype, x, y, a, b, c, d) :
    w = (b - a)
    ell = (d - c)
    s = (w + ell) / 2
    if (ftype == "peaks and valleys") or (ftype == "1") :
        z = np.cos(2*np.pi * x / (w/2)) * np.sin(2*np.pi * y / (ell/2))
    elif (ftype == "narrow stripes") or (ftype == "2") :
        z = np.cos(2*np.pi * x / (w/2)) * np.cos(2*np.pi * y / (ell/2)) \
        +   np.sin(2*np.pi * x / (w/2)) * np.sin(2*np.pi * y / (ell/2))
    elif (ftype == "wide stripes") or (ftype == "3") :
        z = np.cos(2*np.pi * x / w) * np.cos(2*np.pi * y / ell) \
        +   np.sin(2*np.pi * x / w) * np.sin(2*np.pi * y / ell)
    elif (ftype == "bells") or (ftype == "4") :
        z = np.exp(-(1/(.16*s))**2 * ((x - (a + .20*w))**2 + (y - (c + .70*ell))**2)) \
        +   np.exp(-(1/(.16*s))**2 * ((x - (a + .64*w))**2 + (y - (c + .10*ell))**2)) \
        +   np.exp(-(1/(.10*s))**2 * ((x - (a + .82*w))**2 + (y - (c + .76*ell))**2))
    else :
        s = "Invalid choice for ftype.  Try again using one of these:\n"
        s += "(1) \"peaks and valleys\"\n"
        s += "(2) \"narrow stripes\"\n"
        s += "(3) \"wide stripes\"\n"
        s += "(4) \"bells\"\n"
        raise ValueError(s)
    return z

################################################################################

IO.saveArray(os.path.join(dataDir, "f.txt"), eff(ftype, x, y, a, b, c, d))

IO.saveArray(os.path.join(dataDir, "fe.txt"), eff(ftype, xe, ye, a, b, c, d))

