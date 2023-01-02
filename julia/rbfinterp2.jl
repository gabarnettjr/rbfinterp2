# windows

# Uses julia to do the actual calculations for interpolation.

# Greg Barnett
# January 2023

################################################################################

include("io.jl")
include("rbf2.jl")

################################################################################

# Load everything.

dataDir = ARGS[1]

x  = io_loadArray(dataDir * "\\..\\x.txt")
y  = io_loadArray(dataDir * "\\..\\y.txt")
f  = io_loadArray(dataDir * "\\f.txt")
xe = io_loadArray(dataDir * "\\..\\xe.txt")
ye = io_loadArray(dataDir * "\\..\\ye.txt")

rbfPow = parse(Int, ARGS[2])
deg = parse(Int, ARGS[3])
nSubd = parse(Int, ARGS[4])
mSubd = parse(Int, ARGS[5])

################################################################################

# Interpolate using polyharmonic spline (PHS) radial basis functions (RBFs).
computeTime = time()
fe_approx = rbf2_interp(x, y, f, xe, ye; rbfPow=rbfPow, deg=deg, nSubd=nSubd, mSubd=mSubd)
computeTime = time() - computeTime
print("computeTime = " * string(computeTime))

# Save the interpolated values to a file.
io_saveArray(dataDir * "\\fe_approx.txt", fe_approx)

