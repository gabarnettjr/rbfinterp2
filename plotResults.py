#!/usr/bin/python
"""
Creates one 2x2 array of subplots, each displaying some useful
information about how well fe_approx approximates fe.  For comparison,
the known values at the evaluation points must be given in fe.

Greg Barnett
December 2022
"""
from sys import path, argv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri as mtri

path.append(".")
import IO

plotTriangles = True

dataDir = argv[1]
if argv[2] == "y" :
    checkError = True
elif argv[2] == "n" :
    checkError = False
else :
    exit("Please use either \"y\" or \"n\" for the second input.")

################################################################################

x = IO.loadArray(dataDir + "\\..\\x.txt")
y = IO.loadArray(dataDir + "\\..\\y.txt")
xe = IO.loadArray(dataDir + "\\..\\xe.txt")
ye = IO.loadArray(dataDir + "\\..\\ye.txt")

figNum = 1
f = IO.loadArray(dataDir + "\\f.txt")
fe_approx = IO.loadArray(dataDir + "\\fe_approx.txt")
if checkError :
    fe = IO.loadArray(dataDir + "\\fe.txt")
else :
    fe = ()

ab = np.hstack((x, xe))
a = np.min(ab)
b = np.max(ab)
cd = np.hstack((y, ye))
c = np.min(cd)
d = np.max(cd)

triang = mtri.Triangulation(x, y)                                        # nodes
TRIANG = mtri.Triangulation(xe, ye)                          # evaluation points
ms = 1                                                             # marker size
nc = 9                                                        # number of colors
box = [a, b, c, d]                                             # plotting window
lw = .1                                                             # line width

################################################################################

def getContourLevels(vals, useMeanOf = (), minDiff = 0, nColors = 64) :
    """
    Get the z-values to be used to make the contour levels in the 2D surf plots.
    """
    # if len(useMeanOf) == 0 :
        # useMeanOf = vals
    # m = np.mean(useMeanOf)
    # D = np.max([np.max(vals) - m, m - np.min(vals), minDiff])
    # clevels = np.linspace(m - D, m + D, nColors + 1)
    m = np.min(vals)
    M = np.max(vals)
    clevels = np.linspace(m, M, nColors + 1)
    return clevels

################################################################################

# Get vertices of triangles, which will be plotted over and over again.

# Stack first triangle vertex at the end, so entire triangle is plotted.
tmp = np.hstack((triang.triangles[0], triang.triangles[0][0]))

# Get x and y coordinates at the indices of each triangle in the triangulation.
xt = x[tmp]
yt = y[tmp]
for i in range(1, len(triang.triangles)) :
    tmp = np.hstack((triang.triangles[i], triang.triangles[i][0]))
    xt = np.vstack((xt, x[tmp]))
    yt = np.vstack((yt, y[tmp]))
    
################################################################################

# def plotThings(figNum, errXY, FG, errXYe, titleString) :

if checkError :
    fig = plt.figure(figNum, figsize = (13, 9.5))
    plt.subplots_adjust(top=0.963, bottom=0.041, left=0.012, right=0.988 \
    , hspace=0.145, wspace=0.0)
else :
    fig = plt.figure(figNum, figsize = (19, 9.5))
    plt.subplots_adjust(top=0.978, bottom=0.025, left=0.042, right=0.992 \
    , hspace=0.145, wspace=0.094)

# For interpolation (default)
clevels_e = getContourLevels(np.hstack((f, fe_approx, fe)), nColors = nc)
clevels_a = clevels_e

# # For looking at the first derivative.
# clevels_e = getContourLevels(np.hstack((f, fe)), nColors = nc)
# clevels_a = getContourLevels(fe_approx, nColors = nc)

# How you want to format the min/max values printed in the titles.
fmt = "6.3f"

# The known values on the nodes.
if checkError :
    ax = fig.add_subplot(221)
else :
    ax = fig.add_subplot(121)
cs = ax.tricontourf(triang, f, levels = clevels_e)
if plotTriangles :
    for i in range(len(xt)) :
        ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
ax.plot(x, y, 'ko', markersize = ms)
ax.axis('image')
ax.axis(box)
fig.colorbar(cs)
plt.title(('Input Data [{0:' + fmt + '}, {1:' + fmt + '}]').format(np.min(f), np.max(f)))

# The interpolant evaluated at the eval pts.
if checkError :
    ax = fig.add_subplot(222)
else :
    ax = fig.add_subplot(122)
cs = ax.tricontourf(TRIANG, fe_approx, levels = clevels_a)
if plotTriangles :
    for i in range(len(xt)) :
        ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
ax.plot(x, y, 'ko', markersize = ms)
# ax.plot(xe, ye, 'r.', markersize = ms/2)
ax.axis('image')
ax.axis(box)
fig.colorbar(cs)
plt.title(('RBF Interpolant [{0:' + fmt + '}, {1:' + fmt + '}]').format(np.min(fe_approx), np.max(fe_approx)))

if checkError :

    # The known values on the grid.
    ax = fig.add_subplot(223)
    cs = ax.tricontourf(TRIANG, fe, levels = clevels_e)
    if plotTriangles :
        for i in range(len(xt)) :
            ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
    ax.plot(x, y, 'ko', markersize = ms)
    # ax.plot(xe, ye, 'r.', markersize = ms/2)
    ax.axis('image')
    ax.axis(box)
    fig.colorbar(cs)
    plt.title(('Known Values [{0:' + fmt + '}, {1:' + fmt + '}]').format(np.min(fe), np.max(fe)))

    # The error relative to the known values on the grid.
    tmp = (fe_approx - fe) / np.max(np.abs(fe))
    clevels = getContourLevels(tmp, nColors = nc)
    # , useMeanOf = np.array([0]), minDiff = 0, nColors = nc)

    ax = fig.add_subplot(224)
    cs = ax.tricontourf(TRIANG, tmp, levels = clevels)
    if plotTriangles :
        for i in range(len(xt)) :
            ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
    ax.plot(x, y, 'ko', markersize = ms)
    # ax.plot(xe, ye, 'r.', markersize = ms/2)
    ax.axis('image')
    ax.axis(box)
    fig.colorbar(cs)
    plt.title("Relative Error")

plt.show()
