# windows

package rbf2;

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
# December 2022

################################################################################

use strict;
use warnings;

use lib ".";
use linalg;

################################################################################

sub normalize {
    # Shift and scale coordinates for a well-conditioned linear system.
    
    my $x = linalg::copy(shift);                             # x-coords of nodes
    my $y = linalg::copy(shift);                             # y-coords of nodes
    my $xe = linalg::copy(shift);                         # x-coords of eval pts
    my $ye = linalg::copy(shift);                         # y-coords of eval pts
    
    # Shift so that (0,0) is the center of the computational domain.
    my $xavg = linalg::avg($x);
    my $yavg = linalg::avg($y);
    $x = linalg::scalaradd($x, -$xavg);
    $y = linalg::scalaradd($y, -$yavg);
    $xe = linalg::scalaradd($xe, -$xavg);
    $ye = linalg::scalaradd($ye, -$yavg);

    # Re-scale the x and y coordinates (needed for high order poly).
    my $alp = (linalg::norm($x, "inf") + linalg::norm($y, "inf")) / 2;
    $x = linalg::scalarmul($x, (1/$alp));
    $y = linalg::scalarmul($y, (1/$alp));
    $xe = linalg::scalarmul($xe, (1/$alp));
    $ye = linalg::scalarmul($ye, (1/$alp));
    
    return ($x, $y, $xe, $ye);
}

################################################################################

sub jostle {
	# Create "jostled" (not corners) Cartesian nodes, which are randomly moved.
	
	my $nx = shift;            # number of nodes going across (columns of nodes)
	my $ny = shift;                 # number of nodes going down (rows of nodes)
	my $alp = shift;                   # max proportion of node spacing for move
	my $a = shift;                                               # left boundary
	my $b = shift;                                              # right boundary
	my $c = shift;                                             # bottom boundary
	my $d = shift;                                                # top boundary
	
	my $xx = linalg::linspace($a, $b, $nx);
	my $yy = linalg::linspace($c, $d, $ny);
	($xx, $yy) = linalg::meshgrid($xx, $yy);
	$xx = linalg::flatten($xx);
	$yy = linalg::flatten($yy);
	
	my $eps = .0001 * (($b - $a) + ($d - $c)) / 2;
	$alp = (($b - $a) / ($nx - 1) * $alp + ($d - $c) / ($ny - 1) * $alp) / 2;
	
	my $ne = scalar @{$xx};
	
	my ($i, $x, $y, $r);
	
	for ($i = 0; $i < $ne; $i++) {
		$x = @{$xx}[$i];
		$y = @{$yy}[$i];
		if (($x > ($a + $eps)) && ($x < ($b - $eps)) && ($y > ($c + $eps))
		&& ($y < ($d - $eps))) {
			# Interior
			$r = (-$alp + 2 * $alp * rand);
			@{$xx}[$i] += $r;
			$r = (-$alp + 2 * $alp * rand);
			@{$yy}[$i] += $r;
		} elsif ((($x < ($a + $eps)) || ($x > ($b - $eps))) && ($y > ($c + $eps))
		&& ($y < ($d - $eps))) {
			# Left and right boundaries, but no corners.
			$r = (-$alp + 2 * $alp * rand);
			@{$yy}[$i] += $r;
		} elsif ((($y < ($c + $eps)) || ($y > ($d - $eps))) && ($x > ($a + $eps))
		&& ($x < ($b - $eps))) {
			# Top and bottom boundaries, but no corners.
			$r = (-$alp + 2 * $alp * rand);
			@{$xx}[$i] += $r;
		}
	}
	
	return ($xx, $yy);
}

################################################################################

sub inrectangle {
    # Find index of all points that lie in a particular rectangular subdomain.

    my $x = shift;                                # pointer to array of x-coords
    my $y = shift;                                # pointer to array of y-coords
    my $xmci = shift;        # single x-coord of center of rectangular subdomain
    my $ymci = shift;        # single y-coord of center of rectangular subdomain
    my $ell = shift;                           # length of rectangular subdomain
    my $w = shift;                              # width of rectangular subdomain

    my @ind = ();
    my $n = scalar @{$x};

    for (my $j = 0; $j < $n; $j++) {
        if ((abs (@{$x}[$j] - $xmci)) <= $w) {
            if ((abs (@{$y}[$j] - $ymci)) <= $ell) {
                push @ind, $j;
            }
        }
    }
    
    return \@ind;
}

################################################################################

sub rectangles {
    # Find the center (xmc, ymc) and dimensions of each rectangular subdomain.
    
    my $x = shift;                                           # x-coords of nodes
    my $y = shift;                                           # y-coords of nodes
    my $xe = shift;                                       # x-coords of eval pts
    my $ye = shift;                                       # y-coords of eval pts
    my ($nSubd, $mSubd, $deg);
    $nSubd = $mSubd = $deg = -1;
    if ((scalar @_) == 2) {
        $nSubd = shift;                      # number of subdomains horizontally
        $mSubd = shift;                        # number of subdomains vertically
    } elsif ((scalar @_) == 1) {
        $deg = shift;                       # highest degree polynomial included
    } else {
        print STDERR "Must supply either (nSubd, mSubd) or deg, but not both.\n"; die;
    }

    # The rectangular computational domain, [a,b] x [c,d].
    my @ab = (@{$x}, @{$xe});
    my $a = linalg::min(\@ab);
    my $b = linalg::max(\@ab);
    my @cd = (@{$y}, @{$ye});
    my $c = linalg::min(\@cd);
    my $d = linalg::max(\@cd);

    # print ($a . " " . $b . " " . $c . " " . $d . "\n");
    
    # Initial number of subdomains across and down is small.
    my $np;
    if ($deg != -1) {
        if (linalg::norm(\@cd, "inf") > linalg::norm(\@ab, "inf")) {
            $nSubd = 2;
            $mSubd = int (linalg::norm(\@cd, "inf") / linalg::norm(\@ab, "inf") * 2 + .5);
        } else {
            $mSubd = 2;
            $nSubd = int (linalg::norm(\@ab, "inf") / linalg::norm(\@cd, "inf") * 2 + .5);
        }
        # Number of polynomial functions.
        $np = int (($deg + 1) * ($deg + 2) / 2 + .5);
    }
    
    while ("true") {

        # (xmc,ymc) are coordinates of the center of each rectangular subdomain.
        my $eps = .0001 * (($b - $a) + ($d - $c)) / 2;
        my $dx = ($b - $a + 2*$eps) / $nSubd;
        my $dy = ($d - $c + 2*$eps) / $mSubd;
        my $xmc = linalg::linspace($a - $eps + $dx/2, $b + $eps - $dx/2, $nSubd);
        my $ymc = linalg::linspace($c - $eps + $dy/2, $d + $eps - $dy/2, $mSubd);
        ($xmc, $ymc) = linalg::meshgrid($xmc, $ymc);
        $xmc = linalg::flatten($xmc);
        $ymc = linalg::flatten($ymc);

        # Half-width and half-length of each rectangular subdomain.
        my $w = ($b - $a + 2*$eps) / $nSubd / 2;
        my $ell = ($d - $c + 2*$eps) / $mSubd / 2;
        
        # If nSubd and mSubd are function inputs, then return the results now.
        if ($deg == -1) {
            print ("$nSubd x $mSubd subdomains\n");
            return ($xmc, $ymc, $w, $ell);
        }
        
        # Find minimum number of nodes in a subdomain or adjacent subdomains.
        my $minNodes = 99;
        for (my $i = 0; $i < (scalar @{$xmc}); $i++) {
            my $ind = rbf2::inrectangle($x, $y, @{$xmc}[$i], @{$ymc}[$i], 3*$ell, 3*$w);
            if ((scalar @{$ind}) < $minNodes) {
                $minNodes = (scalar @{$ind});
            }
        }
        
        # Quit if the minimum number of nodes gets small enough.
        if ($minNodes < (4 * $np)) {
            print ("$nSubd x $mSubd subdomains\n");
            return ($xmc, $ymc, $w, $ell);
        } else {
            $nSubd += 1;
            $mSubd += 1;
        }
    }
}

################################################################################

sub polymat {
	# Make a polynomial matrix with basis functions arranged in columns.
	
    my $x = shift;                                           # x-coords of input
    my $y = shift;                                           # y-coords of input
    my $deg = shift;                   # largest poly degree to include in basis

    my $n = scalar @{$x};
    my @p = ();
    for (my $i = 0; $i < $n; $i++) {
        my @tmp = (1);
        push @p, \@tmp;
    }

    if ($deg < 0 || $deg > 4) {
        print STDERR "Use a polynomial degree from 0 (constant) up to 4, please.\n"; die;
    }

    for (my $i = 0; $i < $n; $i++) {
        my $xi = @{$x}[$i];
        my $yi = @{$y}[$i];
        if ($deg >= 1) {
            push @{$p[$i]}, $xi, $yi;
        }
        if ($deg >= 2) {
            push @{$p[$i]}, $xi**2, $xi*$yi, $yi**2;
        }
        if ($deg >= 3) {
            push @{$p[$i]}, $xi**3, $xi**2*$yi, $xi*$yi**2, $yi**3;
        }
        if ($deg >= 4) {
            push @{$p[$i]}, $xi**4, $xi**3*$yi, $xi**2*$yi**2, $xi*$yi**3, $yi**4;
        }
    }
    
    return \@p;
}

################################################################################

sub phs {
	# Evaluate a polyharmonic spline function.
	
    my $x = shift;                                            # x-coord of input
    my $y = shift;                                            # y-coord of input
    my $rbfPow = shift;                                # exponent in the phs rbf

    if (! ref $x && ! ref $y) {
        return ($x**2 + $y**2) ** ($rbfPow/2);
    } elsif (ref $x && ref $y && (scalar @{$x} == scalar @{$y})
    && ! ref @{$x}[0] && ! ref @{$y}[0]) {
        my $nE = scalar @{$x};
        my @z = linalg::alloc($nE);
        for (my $i = 0; $i < $nE; $i++) {
            $z[$i] = linalg::phs(@{$x}[$i], @{$y}[$i], $rbfPow);
        }
        return \@z;
    } else {
        print STDERR "Bad input.  Two equal-length arrays or two scalars, please.\n"; die;
    }
}

################################################################################

sub rbfmat {
    # RBF matrix with basis functions arranged in columns.

    my $x = shift;                             # pointer to x-coords of eval pts
    my $y = shift;                             # pointer to y-coords of eval pts
    my $xc = shift;                         # pointer to x-coords of rbf centers
    my $yc = shift;                         # pointer to y-coords of rbf centers
    my $rbfPow = shift;                                # exponent in the phs rbf

    my $nRows = scalar @{$x};
    my $nCols = scalar @{$xc};
    
    my $A = linalg::alloc($nRows, $nCols);
    for (my $i = 0; $i < $nRows; $i++) {
        for (my $j = 0; $j < $nCols; $j++) {
            @{@{$A}[$i]}[$j] = rbf2::phs(@{$x}[$i] - @{$xc}[$j]
			, @{$y}[$i] - @{$yc}[$j], $rbfPow);
        }
    }
    
    return $A;
}

################################################################################

sub interp {
    # Interpolate (x,y,f) to (X,Y,F) using PHS RBFs and polynomials.

    my $x = shift;             # pointer to x-coords where you KNOW the function
    my $y = shift;             # pointer to y-coords where you KNOW the function
    my $f = shift;                # pointer to known values of function on nodes
    my $X = shift;             # pointer to x-coords where you WANT the function
    my $Y = shift;             # pointer to y-coords where you WANT the function
    # Optional:
    my ($rbfPow, $deg, $nSubd, $mSubd);
    if ((scalar @_) == 2 || (scalar @_) == 4) {
        $rbfPow = shift;                                   # exponent of phs rbf
        $deg = shift;                       # largest polynomial degree in basis
        if (scalar @_) {
            $nSubd = shift;                  # number of subdomains horizontally
            $mSubd = shift;                    # number of subdomains vertically
        }
    } elsif (! (scalar @_)) {
        $rbfPow = 3;
        $deg = 1;
    } else {
        print STDERR "Please use either 0, 2, or 4 optional inputs.\n"; die;
    }

    print ("RBF = r**$rbfPow, polynomials up to degree $deg are included.\n");
    
    # Normalize coordinates for good conditioning.
    ($x, $y, $X, $Y) = rbf2::normalize($x, $y, $X, $Y);
    
    # Info (coords, half-width, half-length) about the rectangular subdomains.
    my ($xmc, $ymc, $w, $ell);
    if (($nSubd != -1) && ($mSubd != -1)) {
        ($xmc, $ymc, $w, $ell) = rbf2::rectangles($x, $y, $X, $Y, $nSubd, $mSubd);
    } else {
        ($xmc, $ymc, $w, $ell) = rbf2::rectangles($x, $y, $X, $Y, $deg);
    }

    # Set up a few helper variables.
    my $numP = int (($deg + 1) * ($deg + 2) / 2 + .5);
    my $zp1 = linalg::zeros($numP);
    my $zp2 = linalg::zeros($numP, $numP);
    my $F = linalg::alloc(scalar @{$X});
    my $nSubdomains = scalar @{$xmc};

    for (my $i = 0; $i < $nSubdomains; $i++) {
        
        # Get all nodes in the rectangular subdomain or adjacent subdomains.
        my $ind = rbf2::inrectangle($x, $y, @{$xmc}[$i], @{$ymc}[$i], 3*$ell, 3*$w);
        if (scalar @{$ind} < int (1.5 * $numP + .5)) {
            print ("numLocalNodes = " . (scalar @{$ind}) . "\n");
            print STDERR "Not enough data for this polynomial degree.\n"; die;
        }
        my $xind = linalg::get($x, $ind);
        my $yind = linalg::get($y, $ind);
        
        # Make the polynomial matrix (tall and skinny).
        my $p = rbf2::polymat($xind, $yind, $deg);
        
        # Make the rbf matrix (square).
        my $A = rbf2::rbfmat($xind, $yind, $xind, $yind, $rbfPow);
        
        # Put them together to create the combined rbf-poly matrix (square).
        $A = linalg::hstack($A, $p);
        $p = linalg::hstack(linalg::transpose($p), $zp2);
        $A = linalg::vstack($A, $p);
        
        # Get function values and solve for coefficients, $lam.
        my $lam = linalg::get($f, $ind);
		push @{$lam}, @{$zp1};
        $lam = linalg::solve($A, $lam);
        
        # Find evaluation points in the rectangular subdomain.
        my $IND = rbf2::inrectangle($X, $Y, @{$xmc}[$i], @{$ymc}[$i], $ell, $w);
		if (! (scalar @{$IND})) { next; }
        my $XIND = linalg::get($X, $IND);
        my $YIND = linalg::get($Y, $IND);
        
        # Get rbf-poly evaluation matrix.
        $A = rbf2::rbfmat($XIND, $YIND, $xind, $yind, $rbfPow);
        $p = rbf2::polymat($XIND, $YIND, $deg);
        
        # Evaluate the interpolant at the evaluation points in the subdomain.
        linalg::set($F, $IND, linalg::dot(linalg::hstack($A, $p), $lam));
    }
    
	return $F;
}

################################################################################

return 1;

