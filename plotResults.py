"""
Creates one 2x2 array of subplots, each displaying some useful
information about how well FG approximates errXY.  For comparison,
the known values at the evaluation points must be given in errXYe.
"""
from sys import path, argv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri as mtri

import IO

dataDir = argv[1]
if argv[2] == "1" :
    checkError = True
else :
    checkError = False

################################################################################

x = IO.loadArray(dataDir + "\\..\\x.txt")
y = IO.loadArray(dataDir + "\\..\\y.txt")
xe = IO.loadArray(dataDir + "\\..\\xe.txt")
ye = IO.loadArray(dataDir + "\\..\\ye.txt")

figNum = 1
errXY = IO.loadArray(dataDir + "\\f.txt")
FG = IO.loadArray(dataDir + "\\fe_approx.txt")
if checkError :
    errXYe = IO.loadArray(dataDir + "\\fe.txt")
else :
    errXYe = ()

a = np.min(np.hstack((x, xe)))
b = np.max(np.hstack((x, xe)))
c = np.min(np.hstack((y, ye)))
d = np.max(np.hstack((y, ye)))

triang = mtri.Triangulation(x, y)                                        # nodes
TRIANG = mtri.Triangulation(xe, ye)                          # evaluation points
ms = 1                                                             # marker size
nc = 9                                                        # number of colors
box = [a, b, c, d]                                             # plotting window
lw = .1

################################################################################

def getContourLevels(vals, useMeanOf = (), minDiff = 2, nColors = 64) :
    """
    Get the z-values to be used to make the contour levels in the 2D surf plots.
    """
    if len(useMeanOf) == 0 :
        useMeanOf = vals
    m = np.mean(useMeanOf)
    D = np.max([np.max(vals) - m, m - np.min(vals), minDiff])
    clevels = np.linspace(m - D, m + D, nColors + 1)
    return clevels

################################################################################

# Get vertices of triangles, which will be plotted over and over again.

tmp = np.hstack((triang.triangles[0], triang.triangles[0][0]))
xt = x[tmp]
yt = y[tmp]
for i in range(1, len(triang.triangles)) :
    tmp = np.hstack((triang.triangles[i], triang.triangles[i][0]))
    xt = np.vstack((xt, x[tmp]))
    yt = np.vstack((yt, y[tmp]))
    
plotTriangles = True
    
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

clevels = getContourLevels(np.hstack((errXY, FG, errXYe))  \
, useMeanOf = np.array([0]), minDiff = 0, nColors = nc)

# The known values on the nodes.
if checkError :
    ax = fig.add_subplot(221)
else :
    ax = fig.add_subplot(121)
cs = ax.tricontourf(triang, errXY, levels = clevels)
if plotTriangles :
    for i in range(len(xt)) :
        ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
ax.plot(x, y, 'ko', markersize = ms)
ax.axis('image')
ax.axis(box)
fig.colorbar(cs)
plt.title("Triangulation of Input Data at Nodes")

# The interpolant evaluated at the eval pts.
if checkError :
    ax = fig.add_subplot(222)
else :
    ax = fig.add_subplot(122)
cs = ax.tricontourf(TRIANG, FG, levels = clevels)
if plotTriangles :
    for i in range(len(xt)) :
        ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
ax.plot(x, y, 'ko', markersize = ms)
ax.plot(xe, ye, 'r.', markersize = ms/2)
ax.axis('image')
ax.axis(box)
fig.colorbar(cs)
plt.title("RBF Interpolant at Eval Pts")

if checkError :

    # The known values on the grid.
    ax = fig.add_subplot(223)
    cs = ax.tricontourf(TRIANG, errXYe, levels = clevels)
    if plotTriangles :
        for i in range(len(xt)) :
            ax.plot(xt[i], yt[i], 'k-', linewidth = lw)
    ax.plot(x, y, 'ko', markersize = ms)
    # ax.plot(xe, ye, 'r.', markersize = ms/2)
    ax.axis('image')
    ax.axis(box)
    fig.colorbar(cs)
    plt.title("Known Values at Eval Pts")

    # The error relative to the known values on the grid.
    tmp = (FG - errXYe) / np.max(np.abs(errXYe))
    clevels = getContourLevels(tmp \
    , useMeanOf = np.array([0]), minDiff = 0, nColors = nc)

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
    plt.title("Relative Error at Eval Pts")

plt.show()
