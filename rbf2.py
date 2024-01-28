#!/usr/bin/python
"""
This package contains subroutines for two-dimensional RBF interpolation using
polyharmonic spline (PHS) radial basis functions (RBFs), together with
polynomial functions up to a specified degree.  The main subroutine is called
"interp", but it uses the other subroutines to achieve a local, yet accurate
approximation.  In the typical situation, the user will have nodes (x, y), and
corresponding known function values f.  The goal is to use these known values
to predict the value of the function on some other set of points (xe, ye),
called "evaluation points".  If the user is not interested in adjusting
parameters, then they can use the interp function by passing in five inputs:
x, y, f, xe, ye.
The smallest rectangle that contains all nodes (x,y) will become the overall
computational domain.  If the user does not supply extra inputs, then the
domain will automatically be broken into a number of rectangular subdomains,
so that many small "local" interpolation problems can be solved rather than
one large "global" problem.  The main computational effort is solving for the
coefficients that determine how much of each basis function are needed to
match the function values at the nodes.

Greg Barnett
January 2023
"""
################################################################################

import numpy as np

################################################################################

def normalize(x, y, xe, ye) :
    """
    Shift and scale coordinates for a well-conditioned linear system.
    """
    # x                                                        x-coords of nodes
    # y                                                        y-coords of nodes
    # xe                                           x-coords of evaluation points
    # ye                                           y-coords of evaluation points
    
    # Shift so that (0,0) is the center of the computational domain.
    xavg = np.sum(x) / len(x)
    yavg = np.sum(y) / len(y)
    x = x - xavg
    y = y - yavg
    xe = xe - xavg
    ye = ye - yavg
    
    # Re-scale the x and y coordinates (needed for high order poly).
    alp = (np.linalg.norm(x, np.inf) + np.linalg.norm(y, np.inf)) / 2
    x = x / alp
    y = y / alp
    xe = xe / alp
    ye = ye / alp
    
    return x, y, xe, ye

################################################################################

def jostle(nx, ny, alp, a, b, c, d) :
    """
    Create "jostled" (not corners) Cartesian nodes, which are randomly moved.
    """
    # nx                         number of nodes going across (columns of nodes)
    # ny                              number of nodes going down (rows of nodes)
    # alp                                max proportion of node spacing for move
    # a                                                            left boundary
    # b                                                           right boundary
    # c                                                          bottom boundary
    # d                                                             top boundary
    
    xx = np.linspace(a, b, nx)
    yy = np.linspace(c, d, ny)
    xx, yy = np.meshgrid(xx, yy)
    xx = xx.flatten()
    yy = yy.flatten()
    
    eps = 0.0001 * ((b - a) + (d - c)) / 2
    alp = ((b - a) / (nx - 1) * alp + (d - c) / (ny - 1) * alp) / 2
    
    ne = len(xx)
    
    for i in range(ne) :
        x = xx[i]
        y = yy[i]
        if (x > (a + eps)) and (x < (b - eps)) and (y > (c + eps)) and (y < (d - eps)) :
            # Interior
            r = -alp + 2 * alp * np.random.rand()
            xx[i] = xx[i] + r
            r = -alp + 2 * alp * np.random.rand()
            yy[i] = yy[i] + r
        elif ((x < (a + eps)) or (x > (b - eps))) and (y > (c + eps)) and (y < (d - eps)) :
            # Left and right boundaries, but no corners.
            r = -alp + 2 * alp * np.random.rand()
            yy[i] = yy[i] + r
        elif ((y < (c + eps)) or (y > (d - eps))) and (x > (a + eps)) and (x < (b - eps)) :
            # Top and bottom boundaries, but no corners.
            r = -alp + 2 * alp * np.random.rand()
            xx[i] = xx[i] + r
    
    return xx, yy

################################################################################

def inrectangle(x, y, xmci, ymci, ell, w) :
    """
    Find index of all points that lie in a particular rectangular subdomain.
    """
    # x                                                        array of x-coords
    # y                                                        array of y-coords
    # xmci                     single x-coord of center of rectangular subdomain
    # ymci                     single y-coord of center of rectangular subdomain
    # ell                                        length of rectangular subdomain
    # w                                           width of rectangular subdomain

    ind = []

    for j in range(len(x)) :
        if abs(x[j] - xmci) <= w :
            if abs(y[j] - ymci) <= ell :
                ind.append(j)
    
    return ind

################################################################################

def rectangles(x, y, xe, ye, nSubd=-1, mSubd=-1, deg=-1) :
    """
    Find the center (xmc, ymc) and dimensions of each rectangular subdomain.
    """
    # x                                                        x-coords of nodes
    # y                                                        y-coords of nodes
    # xe                                                    x-coords of eval pts
    # ye                                                    y-coords of eval pts
    if (nSubd != -1) and (mSubd != -1) :
        deg = -1                                             # polynomial degree
    elif deg != -1 :
        nSubd = -1                           # number of subdomains going across
        mSubd = -1                             # number of subdomains going down
    else :
        s = "Must supply either (nSubd, mSubd) or deg, but not both."
        raise ValueError(s)

    # The rectangular computational domain, [a,b] x [c,d].
    ab = np.hstack((x, xe))
    a = np.min(ab)
    b = np.max(ab)
    cd = np.hstack((y, ye))
    c = np.min(cd)
    d = np.max(cd)

    # Initial number of subdomains across and down is small.
    if deg != -1 :
        if (np.linalg.norm(cd, np.inf) > np.linalg.norm(ab, np.inf)) :
            nSubd = 2
            mSubd = int(round(np.linalg.norm(cd, np.inf) / np.linalg.norm(ab, np.inf) * 2))
        else :
            mSubd = 2
            nSubd = int(round(np.linalg.norm(ab, np.inf) / np.linalg.norm(cd, np.inf) * 2))
        # Number of polynomial functions.
        numP = int(round((deg + 1) * (deg + 2) / 2))
    
    while True :

        # (xmc,ymc) are coordinates of the center of each rectangular subdomain.
        eps = 0.0001 * ((b - a) + (d - c)) / 2
        dx = (b - a + 2*eps) / nSubd
        dy = (d - c + 2*eps) / mSubd
        xmc = np.linspace(a - eps + dx/2, b + eps - dx/2, nSubd)
        ymc = np.linspace(c - eps + dy/2, d + eps - dy/2, mSubd)
        xmc, ymc = np.meshgrid(xmc, ymc)
        xmc = xmc.flatten()
        ymc = ymc.flatten()

        # Half-width and half-length of each rectangular subdomain.
        w = (b - a + 2*eps) / nSubd / 2
        ell = (d - c + 2*eps) / mSubd / 2
        
        # If nSubd and mSubd are function inputs, then return the results now.
        if deg == -1 :
            print('{0:1d} x {1:1d} subdomains'.format(nSubd, mSubd))
            return xmc, ymc, w, ell
        
        # Find minimum number of nodes in a subdomain or adjacent subdomains.
        minNodes = 999999
        for i in range(len(xmc)) :
            ind = inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
            if len(ind) < minNodes :
                minNodes = len(ind)
        
        # Quit if the minimum number of nodes gets small enough.
        if minNodes < 10 * numP :
            print('{0:1d} x {1:1d} subdomains'.format(nSubd, mSubd))
            return xmc, ymc, w, ell
        else :
            nSubd *= 2
            mSubd *= 2

################################################################################

def polymat(x, y, deg, kind="i") :
    """
	Make a polynomial matrix with basis functions arranged in rows.
	"""
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # deg                                largest poly degree to include in basis

    if (deg < 0) or (deg > 4) :
        s = "Use a polynomial degree from 0 (constant) up to 4, please."
        raise ValueError(s)
    
    numPoly = int(round((deg + 1) * (deg + 2) / 2))
    p = np.zeros((numPoly, len(x)))
    
    if kind == "i" :
        p[0,:] = 1
        if deg >= 1 :
            p[1,:] = x
            p[2,:] = y
        if deg >= 2 :
            p[3,:] = x**2
            p[4,:] = x*y
            p[5,:] = y**2
        if deg >= 3 :
            p[6,:] = x**3
            p[7,:] = x**2*y
            p[8,:] = x*y**2
            p[9,:] = y**3
        if deg >= 4 :
            p[10,:] = x**4
            p[11,:] = x**3*y
            p[12,:] = x**2*y**2
            p[13,:] = x*y**3
            p[14,:] = y**4
    elif kind == "x" :
        p[0,:] = 0
        if deg >= 1 :
            p[1,:] = 1
            p[2,:] = 0
        if deg >= 2 :
            p[3,:] = 2*x
            p[4,:] = y
            p[5,:] = 0
        if deg >= 3 :
            p[6,:] = 3*x**2
            p[7,:] = 2*x*y
            p[8,:] = y**2
            p[9,:] = 0
        if deg >= 4 :
            p[10,:] = 4*x**3
            p[11,:] = 3*x**2*y
            p[12,:] = 2*x*y**2
            p[13,:] = y**3
            p[14,:] = 0
    elif kind == "y" :
        p[0,:] = 0
        if deg >= 1 :
            p[1,:] = 0
            p[2,:] = 1
        if deg >= 2 :
            p[3,:] = 0
            p[4,:] = x
            p[5,:] = 2*y
        if deg >= 3 :
            p[6,:] = 0
            p[7,:] = x**2
            p[8,:] = x*2*y
            p[9,:] = 3*y**2
        if deg >= 4 :
            p[10,:] = 0
            p[11,:] = x**3
            p[12,:] = x**2*2*y
            p[13,:] = x*3*y**2
            p[14,:] = 4*y**3
    else :
        s = "Optional variable \"kind\" should be \"i\", \"x\", or \"y\"."
        raise ValueError(s)
    
    return p

################################################################################

def phs(x, y, rbfPow) :
    """
	Evaluate a polyharmonic spline basis function.
	"""
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # rbfPow                                             exponent in the phs rbf

    return (x**2 + y**2) ** (rbfPow/2)

################################################################################

def phs_x(x, y, rbfPow) :
    """
	Evaluate the derivative of a PHS basis function with respect to x.
	"""
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # rbfPow                                             exponent in the phs rbf

    return (rbfPow*x) * (x**2 + y**2) ** ((rbfPow - 2)/2)

################################################################################

def phs_y(x, y, rbfPow) :
    """
	Evaluate the derivative of a PHS basis function with respect to y.
	"""
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # rbfPow                                             exponent in the phs rbf

    return (rbfPow*y) * (x**2 + y**2) ** ((rbfPow - 2)/2)

################################################################################

def rbfmat(x, y, xc, yc, rbfPow, func=phs) :
    """
    RBF matrix with basis functions arranged in columns.
    """
    # x                                                     x-coords of eval pts
    # y                                                     y-coords of eval pts
    # xc                                                 x-coords of rbf centers
    # yc                                                 y-coords of rbf centers
    # rbfPow                                             exponent in the phs rbf

    nRows = len(x)
    nCols = len(xc)
    
    A = np.zeros((nRows, nCols))
    for i in range(nRows) :
        for j in range(nCols) :
            A[i,j] = func(x[i] - xc[j], y[i] - yc[j], rbfPow)
    
    return A

################################################################################

def interp(x, y, f, xe, ye, rbfPow=-1, deg=-1, nSubd=-1, mSubd=-1) :
    """
    Interpolate (x,y,f) to (xe,ye,fe_approx) using PHS RBFs and polynomials.
    """
    # x                                     x-coords where you KNOW the function
    # y                                     y-coords where you KNOW the function
    # f                                        known values of function on nodes
    # xe                                    x-coords where you WANT the function
    # ye                                    y-coords where you WANT the function
    # OPTIONAL:
    # rbfPow                                                 exponent of phs rbf
    # deg                                     largest polynomial degree in basis
    # nSubd                                    number of subdomains horizontally
    # mSubd                                      number of subdomains vertically
    if (rbfPow == -1) and (deg == -1) :
        rbfPow = 3
        deg = 1
    
    if rbfPow == -1 :
        print('RBF = NONE, polynomials up to degree {0:1d} are included.'.format(deg))
    else :
        print('RBF = r**{0:1d}, polynomials up to degree {1:1d} are included.'.format(rbfPow, deg))
    
    # Normalize coordinates for good conditioning.
    x, y, xe, ye = normalize(x, y, xe, ye)
    
    # Info (coords, half-width, half-length) about the rectangular subdomains.
    if (nSubd != -1) and (mSubd != -1) :
        xmc, ymc, w, ell = rectangles(x, y, xe, ye, nSubd=nSubd, mSubd=mSubd)
    elif deg != -1 :
        xmc, ymc, w, ell = rectangles(x, y, xe, ye, deg=deg)
    else :
        s = "Need either (nSubd,mSubd) or deg, or both."
        raise ValueError(s)
    
    # Set up a few helper variables.
    numP = int(round((deg + 1) * (deg + 2) / 2))
    zp1 = np.zeros((numP, 1))
    zp2 = np.zeros((numP, numP))
    fe_approx = np.zeros(len(xe))
    
    for i in range(len(xmc)) :

        # Get all nodes in the rectangular subdomain or adjacent subdomains.
        ind = inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
        if len(ind) < round(1.5 * numP) :
            print('numLocalNodes = {0:2d}'.format(len(ind)))
            s = "Not enough data for this polynomial degree."
            raise ValueError(s)
        xind = x[ind]
        yind = y[ind]

        # Make the polynomial matrix.
        p = polymat(xind, yind, deg)

        # Find evaluation points in the rectangular subdomain.
        IND = inrectangle(xe, ye, xmc[i], ymc[i], ell, w)
        if len(IND) == 0 :
            continue
        xeIND = xe[IND]
        yeIND = ye[IND]

        # Put together the RBF-poly approximant at the evaluation points.
        if (rbfPow == -1) :
            # Just do regular polynomial least squares.
            lam = np.linalg.lstsq(p.T, f[ind], rcond=None)[0]
            p = polymat(xeIND, yeIND, deg, kind="i")
            fe_approx[IND] = p.T.dot(lam).flatten()
        else :
            # Make the rbf matrix (square).
            A = rbfmat(xind, yind, xind, yind, rbfPow)
            # Put them together to create the combined rbf-poly matrix (square).
            A = np.hstack((A, p.T))
            p = np.hstack((p, zp2))
            A = np.vstack((A, p))
            # Get function values and solve for coefficients, $lam.
            lam = np.zeros((len(ind), 1))
            lam[:,0] = f[ind]
            lam = np.vstack((lam, zp1))
            lam = np.linalg.solve(A, lam)
            # Get rbf-poly evaluation matrix.
            A = rbfmat(xeIND, yeIND, xind, yind, rbfPow, func=phs)
            p = polymat(xeIND, yeIND, deg, kind="i").T
            # Evaluate the interpolant at the evaluation points in the subdomain.
            fe_approx[IND] = np.hstack((A, p)).dot(lam).flatten()
    
    return fe_approx

