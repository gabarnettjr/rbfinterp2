# windows

# Uses python to do the actual calculations for interpolation.

# Greg Barnett
# January 2023

################################################################################

from sys import argv, path
from time import time

path.append(".")
import IO

path.append(".\\python")
import rbf2

################################################################################

# Load everything.

dataDir = argv[1]

x  = IO.loadArray(dataDir + "\\..\\x.txt")
y  = IO.loadArray(dataDir + "\\..\\y.txt")
f  = IO.loadArray(dataDir + "\\f.txt")
xe = IO.loadArray(dataDir + "\\..\\xe.txt")
ye = IO.loadArray(dataDir + "\\..\\ye.txt")

rbfPow = int(argv[2])
deg = int(argv[3])
nSubd = int(argv[4])
mSubd = int(argv[5])

################################################################################

# Interpolate using polyharmonic spline (PHS) radial basis functions (RBFs).
computeTime = time()
fe_approx = rbf2.interp(x, y, f, xe, ye, rbfPow=rbfPow, deg=deg, nSubd=nSubd, mSubd=mSubd)
computeTime = time() - computeTime
print("computeTime = " + str(computeTime))

# Save the interpolated values to a file.
IO.saveArray(dataDir + "\\fe_approx.txt", fe_approx)

