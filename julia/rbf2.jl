# windows

using Printf
using LinearAlgebra

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

function rbf2_meshgrid(x, y)
    # Create matrices xx and yy based on arrays x and y.  x is repeated and
    # stacked vertically, while y is repeated and stacked horizontally.  In the
    # end, the output matrices xx and yy should have numRows = len(y) and
    # numCols = len(x).
    
    # x                                                          array of x-vals
    # y                                                          array of y-vals
    
    xx = repeat(x', length(y), 1)
    yy = repeat(y, 1, length(x))
    
    return (xx, yy)
end

################################################################################

function rbf2_normalize(x, y, xe, ye)
    # Shift and scale coordinates for a well-conditioned linear system.
    
    # x                                                        x-coords of nodes
    # y                                                        y-coords of nodes
    # xe                                           x-coords of evaluation points
    # ye                                           y-coords of evaluation points
    
    # Shift so that (0,0) is the center of the computational domain.
    xavg = sum(x) / length(x)
    yavg = sum(y) / length(y)
    x .= x .- xavg
    y .= y .- yavg
    xe .= xe .- xavg
    ye .= ye .- yavg

    # Re-scale the x and y coordinates (needed for high order poly).
    alp = (norm(x, Inf) + norm(y, Inf)) / 2
    x .= x ./ alp
    y .= y ./ alp
    xe .= xe ./ alp
	ye .= ye ./ alp
    
    return (x, y, xe, ye)
end

################################################################################

function rbf2_jostle(nx, ny, alp, a, b, c, d)
	# Create "jostled" (not corners) Cartesian nodes, which are randomly moved.
	
	# nx                         number of nodes going across (columns of nodes)
	# ny                              number of nodes going down (rows of nodes)
	# alp                                max proportion of node spacing for move
	# a                                                            left boundary
	# b                                                           right boundary
	# c                                                          bottom boundary
	# d                                                             top boundary
	
	xx = linspace(a, b, nx)
	yy = linspace(c, d, ny)
	(xx, yy) = rbf2_meshgrid(xx, yy)
	xx = xx'[:]
	yy = yy'[:]
	
	eps = .001
	alp = ((b - a) / (nx - 1) * alp + (d - c) / (ny - 1) * alp) / 2
	
	ne = length(xx)
	
	for i in 1 : ne
		x = xx[i]
		y = yy[i]
		if (x > (a + eps)) && (x < (b - eps)) && (y > (c + eps)) && (y < (d - eps))
			# Interior
			r = -alp + 2 * alp * rand()
			xx[i] = xx[i] + r
			r = -alp + 2 * alp * rand()
			yy[i] = yy[i] + r
		elseif ((x < (a + eps)) || (x > (b - eps))) && (y > (c + eps)) && (y < (d - eps))
			# Left and right boundaries, but no corners.
			r = -alp + 2 * alp * rand()
			yy[i] = yy[i] + r
		elseif ((y < (c + eps)) || (y > (d - eps))) && (x > (a + eps)) && (x < (b - eps))
			# Top and bottom boundaries, but no corners.
			r = -alp + 2 * alp * rand()
			xx[i] = xx[i] + r
		end
	end
	
	return (xx, yy)
end

################################################################################

function rbf2_inrectangle(x, y, xmci, ymci, ell, w)
    # Find index of all points that lie in a particular rectangular subdomain.

    # x                                                        array of x-coords
    # y                                                        array of y-coords
    # xmci                     single x-coord of center of rectangular subdomain
    # ymci                     single y-coord of center of rectangular subdomain
    # ell                                        length of rectangular subdomain
    # w                                           width of rectangular subdomain

    ind = []
    n = length(x)

    for j in 1 : length(x)
        if abs(x[j] - xmci) <= w
            if abs(y[j] - ymci) <= ell
                ind = vcat(ind, [j])
            end
        end
    end
    
    return ind
end

################################################################################

function rbf2_rectangles(x, y, xe, ye; nSubd=-1, mSubd=-1, deg=-1)
    # Find the center (xmc, ymc) and dimensions of each rectangular subdomain.
    
    # x                                                        x-coords of nodes
    # y                                                        y-coords of nodes
    # xe                                                    x-coords of eval pts
    # ye                                                    y-coords of eval pts
    if (nSubd != -1) && (mSubd != -1)
        deg = -1                                             # polynomial degree
    elseif deg != -1 
        nSubd = -1                           # number of subdomains going across
        mSubd = -1                             # number of subdomains going down
    else
        error("Must supply either (nSubd, mSubd) or deg, but not both.\n")
    end

    # The rectangular computational domain, [a,b] x [c,d].
    ab = vcat(x, xe)
    a = minimum(ab)
    b = maximum(ab)
    cd = vcat(y, ye)
    c = minimum(cd)
    d = maximum(cd)

    # @printf("%f, %f, %f, %f\n", a, b, c, d)
    
    # Initial number of subdomains across and down is small.
    if deg != -1
        if (norm(cd, Inf) > norm(ab, Inf))
            nSubd = 2
            mSubd = Int(round(norm(cd, Inf) / norm(ab, Inf) * 2))
        else
            mSubd = 2
            nSubd = Int(round(norm(ab, Inf) / norm(cd, Inf) * 2))
        end
        # Number of polynomial functions.
        np = Int(round((deg + 1) * (deg + 2) / 2))
    end
    
    while true

        # (xmc,ymc) are coordinates of the center of each rectangular subdomain.
        eps = 0.001
        dx = (b - a + 2*eps) / nSubd
        dy = (d - c + 2*eps) / mSubd
        xmc = range(a - eps + dx/2, b + eps - dx/2, length=nSubd)
        ymc = range(c - eps + dy/2, d + eps - dy/2, length=mSubd)
        (xmc, ymc) = rbf2_meshgrid(xmc, ymc)
        xmc = xmc'[:]
        ymc = ymc'[:]

        # Half-width and half-length of each rectangular subdomain.
        w = (b - a + 2*eps) / nSubd / 2
        ell = (d - c + 2*eps) / mSubd / 2
        
        # If nSubd and mSubd are function inputs, then return the results now.
        if deg == -1
            @printf("%i x %i subdomains\n", nSubd, mSubd)
            return (xmc, ymc, w, ell)
        end
        
        # Find minimum number of nodes in a subdomain or adjacent subdomains.
        minNodes = 99
        for i in 1 : length(xmc)
            ind = rbf2_inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
            if length(ind) < minNodes
                minNodes = length(ind)
            end
        end
        
        # Quit if the minimum number of nodes gets small enough.
        if minNodes < 4 * np
            @printf("%i x %i subdomains\n", nSubd, mSubd)
            return (xmc, ymc, w, ell)
        else
            nSubd += 1
            mSubd += 1
        end
    end
end

################################################################################

function rbf2_polymat(x, y, deg)
	# Make a polynomial matrix with basis functions arranged in columns.
	
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # deg                                largest poly degree to include in basis

    if (deg < 0) || (deg > 4)
        error("Use a polynomial degree from 0 (constant) up to 4, please.\n")
    end
    
    p = []

    if deg >= 0
        p = ones(size(x))
    end
    if deg >= 1
        p = hcat(p, x, y)
    end
    if deg >= 2
        p = hcat(p, x .^ 2, x .* y, y .^ 2)
    end
    if deg >= 3
        p = hcat(p, x .^ 3, x .^ 2 .* y, x .* y .^ 2, y .^ 3)
    end
    if deg >= 4
        p = hcat(p, x .^ 4, x .^ 3 .* y, x .^ 2 .* y .^ 2, x .* y .^ 3, y .^ 4)
    end
    
    return p
end

################################################################################

function rbf2_phs(x, y, rbfPow)
	# Evaluate a polyharmonic spline function.
	
    # x                                                        x-coords of input
    # y                                                        y-coords of input
    # rbfPow                                             exponent in the phs rbf

    return (x .^ 2 + y .^ 2) .^ (rbfPow/2)
end

################################################################################

function rbf2_rbfmat(x, y, xc, yc, rbfPow)
    # RBF matrix with basis functions arranged in columns.

    # x                                                     x-coords of eval pts
    # y                                                     y-coords of eval pts
    # xc                                                 x-coords of rbf centers
    # yc                                                 y-coords of rbf centers
    # rbfPow                                             exponent in the phs rbf

    nRows = length(x)
    nCols = length(xc)
    
    A = zeros(nRows, nCols)
    for i in 1 : nRows
        for j in 1 : nCols
            A[i,j] = rbf2_phs(x[i] - xc[j], y[i] - yc[j], rbfPow)
        end
    end
    
    return A
end

################################################################################

function rbf2_interp(x, y, f, xe, ye; rbfPow=-1, deg=-1, nSubd=-1, mSubd=-1)
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
    if (rbfPow == -1) && (deg == -1)
        rbfPow = 3
        deg = 1
    end
    @printf("RBF is r^%i, polynomials up to degree %i are included.\n", rbfPow, deg)
    
    # Normalize coordinates for good conditioning.
    (x, y, xe, ye) = rbf2_normalize(x, y, xe, ye)
    
    # Info (coords, half-width, half-length) about the rectangular subdomains.
    if (nSubd != -1) && (mSubd != -1)
        (xmc, ymc, w, ell) = rbf2_rectangles(x, y, xe, ye; nSubd=nSubd, mSubd=mSubd)
    else
        (xmc, ymc, w, ell) = rbf2_rectangles(x, y, xe, ye; deg=deg)
    end

    # Set up a few helper variables.
    numP = Int(round((deg + 1) * (deg + 2) / 2))
    zp1 = zeros(numP, 1)
    zp2 = zeros(numP, numP)
    fe_approx = zeros(length(xe))
    nSubdomains = length(xmc)

    for i in 1 : nSubdomains
        
        # Get all nodes in the rectangular subdomain or adjacent subdomains.
        ind = rbf2_inrectangle(x, y, xmc[i], ymc[i], 3*ell, 3*w)
        if length(ind) < Int(round(1.5 * numP))
            @printf("numLocalNodes = %i\n", length(ind))
            error("Not enough data for this polynomial degree.\n")
        end
        xind = x[ind]
        yind = y[ind]
        
        # Make the polynomial matrix (tall and skinny).
        p = rbf2_polymat(xind, yind, deg)
        
        # Make the rbf matrix (square).
        A = rbf2_rbfmat(xind, yind, xind, yind, rbfPow)
        
        # Put them together to create the combined rbf-poly matrix (square).
        A = hcat(A, p)
        p = hcat(p', zp2)
        A = vcat(A, p)
        
        # Get function values and solve for coefficients, $lam.
        lam = f[ind]
		lam = vcat(lam, zp1)
        lam = A \ lam
        
        # Find evaluation points in the rectangular subdomain.
        IND = rbf2_inrectangle(xe, ye, xmc[i], ymc[i], ell, w)
		if length(ind) == 0;  next;  end;
        xeIND = xe[IND]
        yeIND = ye[IND]
        
        # Get rbf-poly evaluation matrix.
        A = rbf2_rbfmat(xeIND, yeIND, xind, yind, rbfPow)
        p = rbf2_polymat(xeIND, yeIND, deg)
        
        # Evaluate the interpolant at the evaluation points in the subdomain.
        fe_approx[IND] = hcat(A, p) * lam
    end
    
	return fe_approx
end

