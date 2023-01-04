# windows

# This package contains subroutines for two-dimensional RBF interpolation using
# polyharmonic spline (PHS) radial basis functions (RBFs), together with
# polynomial functions up to a specified degree.  The main subroutine is called
# "interp", but it uses the other subroutines to achieve a local, yet accurate
# approximation.  In the typical situation, the user will have nodes (x, y), and
# corresponding known function values f.  The goal is to use these known values
# to predict the value of the function on some other set of points (xe, ye),
# called "evaluation points".  If the user is not interested in adjusting
# parameters, then they can use the interp function by passing in five inputs:
# x, y, f, xe, ye.
# The smallest rectangle that contains all nodes (x,y) will become the overall
# computational domain.  If the user does not supply extra inputs, then the
# domain will automatically be broken into a number of rectangular subdomains,
# so that many small "local" interpolation problems can be solved rather than
# one large "global" problem.  The main computational effort is solving for the
# coefficients that determine how much of each basis function are needed to
# match the function values at the nodes.

# Greg Barnett
# January 2023

################################################################################

import numpy as np

################################################################################

def normalize(x, y, xe, ye) :
    # Shift and scale coordinates for a well-conditioned linear system.
    
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
    # Create "jostled" (not corners) Cartesian nodes, which are randomly moved.
    
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
    yy = yy.flatten
    
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
    # Find index of all points that lie in a particular rectangular subdomain.

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
    # Find the center (xmc, ymc) and dimensions of each rectangular subdomain.
    
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
        raise ValueError("Must supply either (nSubd, mSubd) or deg, but not both.\n")

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
            print('{0:1d} x {1:1d} subdomains\n'.format(nSubd, mSubd))
            return xmc, ymc, w, ell
        
        # Find minimum number of nodes in a subdomain or adjacent subdomains.
        minNodes = 99
        for i in range(len(xmc)) :
            ind = inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
            if len(ind) < minNodes :
                minNodes = len(ind)
        
        # Quit if the minimum number of nodes gets small enough.
        if minNodes < 4 * numP :
            print('{0:1d} x {1:1d} subdomains\n'.format(nSubd, mSubd))
            return xmc, ymc, w, ell
        else :
            nSubd += 1
            mSubd += 1

################################################################################

def polymat(x, y, deg) :
	# Make a polynomial matrix with basis functions arranged in rows.
	
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # deg                                largest poly degree to include in basis

    if (deg < 0) or (deg > 4) :
        raise ValueError("Use a polynomial degree from 0 (constant) up to 4, please.\n")
    
    p = np.ones(np.shape(x))
    
    if deg >= 1 :
        p = np.vstack((p, x, y))
    if deg >= 2 :
        p = np.vstack((p, x**2, x*y, y**2))
    if deg >= 3 :
        p = np.vstack((p, x**3, x**2*y, x*y**2, y**3))
    if deg >= 4 :
        p = np.vstack((p, x**4, x**3*y, x**2*y**2, x*y**3, y**4))
    
    return p

################################################################################

def phs(x, y, rbfPow) :
	# Evaluate a polyharmonic spline function.
	
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # rbfPow                                             exponent in the phs rbf

    return (x**2 + y**2) ** (rbfPow/2)

################################################################################

def rbfmat(x, y, xc, yc, rbfPow) :
    # RBF matrix with basis functions arranged in columns.

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
            A[i,j] = phs(x[i] - xc[j], y[i] - yc[j], rbfPow)
    
    return A

################################################################################

def interp(x, y, f, xe, ye, rbfPow=-1, deg=-1, nSubd=-1, mSubd=-1) :
    # Interpolate (x,y,f) to (xe,ye,fe_approx) using PHS RBFs and polynomials.
    
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
    print('RBF is r**{0:1d}, polynomials up to degree {1:1d} are included.\n'.format(rbfPow, deg))
    
    # Normalize coordinates for good conditioning.
    x, y, xe, ye = normalize(x, y, xe, ye)
    
    # Info (coords, half-width, half-length) about the rectangular subdomains.
    if (nSubd != -1) and (mSubd != -1) :
        xmc, ymc, w, ell = rectangles(x, y, xe, ye, nSubd=nSubd, mSubd=mSubd)
    elif deg != -1 :
        xmc, ymc, w, ell = rectangles(x, y, xe, ye, deg=deg)
    else :
        raise ValueError("Need either (nSubd,mSubd) or deg, or both.\n")
    
    # Set up a few helper variables.
    numP = int(round((deg + 1) * (deg + 2) / 2))
    zp1 = np.zeros((numP, 1))
    zp2 = np.zeros((numP, numP))
    fe_approx = np.zeros(len(xe))
    
    for i in range(len(xmc)) :
         
         # Get all nodes in the rectangular subdomain or adjacent subdomains.
         ind = inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
         if len(ind) < int(round(1.5 * numP)) :
             print('numLocalNodes = {0:2d}\n'.format(len(ind)))
             raise ValueError("Not enough data for this polynomial degree.\n")
         xind = x[ind]
         yind = y[ind]
         
         # Make the polynomial matrix.
         p = polymat(xind, yind, deg)
         
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
         
         # Find evaluation points in the rectangular subdomain.
         IND = inrectangle(xe, ye, xmc[i], ymc[i], ell, w)
         if len(IND) == 0 :
             next
         xeIND = xe[IND]
         yeIND = ye[IND]
         
         # Get rbf-poly evaluation matrix.
         A = rbfmat(xeIND, yeIND, xind, yind, rbfPow)
         p = polymat(xeIND, yeIND, deg).T
         
         # Evaluate the interpolant at the evaluation points in the subdomain.
         fe_approx[IND] = np.hstack((A, p)).dot(lam).flatten()
    
    return fe_approx
