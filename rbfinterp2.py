#!/usr/bin/python
"""
This is a python script that can be called to interpolate (x, y, f) to 
(xe, ye, fe_approx).  In order for this script to work, you need to have some
text files saved in a common location, so they can be loaded in by the script.

These four files should be saved in a coordinate directory:
x.txt, y.txt, xe.txt, ye.txt (nodes and evaluation points)

The following should be in a subfolder ($dataDir) of the coordinate directory:
f.txt (function values at the nodes)

In addition, if you wish to compare to some known exact values at the
evaluation points, then you should also include this file in $dataDir:
fe.txt (function values at the evaluation points)

If you have all of this, then the script will produce the file fe_approx.txt,
which will contain estimated values of the function at the evaluation points.

Greg Barnett
January 2023
"""
################################################################################

import os
from sys import argv, path
from time import time

path.append(".")
import IO
import rbf2

################################################################################

def helpString() :
    s = "\n"
    s += "\"rbfinterp2.pl\": A script for interpolating (x,y,f) to (xe,ye,fe_approx).\n\n"
    s += "To run the script, you need a coordinates directory that contains these:\n"
    s += "x.txt, y.txt, xe.txt, ye.txt\n\n"
    s += "Inside the coordinates directory should be a function value subdirectory, containing these:\n"
    s += "f.txt  (required)\n"
    s += "fe.txt (optional, but needed for error calculation)\n\n"
    s += "This script accepts up to 6 command-line inputs:\n"
    s += "(1) The path to the folder that contains your function values  (default: .\\randomCoords\\smoothData).\n"
    s += "(2) Whether or not to calculate the error, y or n              (default: n).\n"
    s += "(3) The rbf exponent, an odd integer                           (default: 3).\n"
    s += "(4) The polynomial degree, an integer from 0 up to 4           (default: 1).\n"
    s += "(5) The number of subdomains across, a positive integer        (default: auto calculate).\n"
    s += "(6) The number of subdomains going down, a positive integer    (default: auto calculate).\n\n"
    return s

################################################################################

# Process the input.

dataDir = os.path.join("randomCoords", "smoothData")
checkError = "y"
rbfPow = 3
deg = 1
nSubd = -1
mSubd = -1

argv = argv[1:]
if len(argv) > 0 :
    tmp = argv[0].lower()
    if ("help" == tmp) or ("--help" == tmp) or ("-h" == tmp) :
        print(helpString())
        exit()
    dataDir = argv[0];  argv = argv[1:]
    if not os.path.isdir(dataDir) :
        raise ValueError("First input must be a data directory.")
    elif not os.path.isfile(os.path.join(dataDir, "f.txt")) :
        raise ValueError("Data directory must contain function values.")
    checkError = "n"

if len(argv) > 0 :
    checkError = argv[0].lower();  argv = argv[1:]
    if ("y" in checkError) and ("n" not in checkError) :
        checkError = "y"
    elif ("n" in checkError) and ("y" not in checkError) :
        checkError = "n"
    else :
        ValueError("Invalid second input (checkError).  Should be \"y\" or \"n\".")
	
if len(argv) > 0 : rbfPow = int(argv[0]);  argv = argv[1:]
if len(argv) > 0 :    deg = int(argv[0]);  argv = argv[1:]
if len(argv) > 0 :  nSubd = int(argv[0]);  argv = argv[1:]
if len(argv) > 0 :  mSubd = int(argv[0]);  argv = argv[1:]

if len(argv) > 0 :
    raise ValueError("Too many inputs.  Max number of inputs is 6.")

################################################################################

if os.path.isfile(os.path.join(dataDir, "fe_approx.txt")) :
    os.remove(os.path.join(dataDir, "fe_approx.txt"))

# Load everything and interpolate in python.
x  = IO.loadArray(os.path.join(dataDir, "..", "x.txt"))
y  = IO.loadArray(os.path.join(dataDir, "..", "y.txt"))
f  = IO.loadArray(os.path.join(dataDir, "f.txt"))
xe = IO.loadArray(os.path.join(dataDir, "..", "xe.txt"))
ye = IO.loadArray(os.path.join(dataDir, "..", "ye.txt"))
computeTime = time()
fe_approx = rbf2.interp(x, y, f, xe, ye, rbfPow=rbfPow, deg=deg, nSubd=nSubd, mSubd=mSubd)
computeTime = time() - computeTime
print("computeTime = " + str(computeTime))
IO.saveArray(os.path.join(dataDir, "fe_approx.txt"), fe_approx)

# # Load everything and interpolate in perl.
# os.system("perl " + os.path.join("perl", "rbfinterp2.pl") + " " + dataDir + " " + \
# checkError + " " + str(rbfPow) + " " + str(deg) + " " + str(nSubd) + " " + str(mSubd))

# # Load everything and interpolate in julia.
# os.system("julia " + os.path.join("julia", "rbfinterp2.jl") + " " + dataDir + " " + \
#  str(rbfPow) + " " + str(deg) + " " + str(nSubd) + " " + str(mSubd))

if not os.path.isfile(os.path.join(dataDir, "fe_approx.txt")) :
    raise ValueError("Please investigate error during fe_approx.txt creation.")

################################################################################

# Use separate script to visualize results.
os.system("python plotResults.py " + dataDir + " " + checkError)

