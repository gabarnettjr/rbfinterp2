# windows

package rbf;

use strict;
use warnings;

use lib ".";
use linalg;

################################################################################

sub inrectangle {
    # Find index of all points that lie in a rectangular subdomain.

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

sub polymat {
	# Make a polynomial matrix with basis functions arranged in columns.
	
    my $x = shift;                                           # x-coords of input
    my $y = shift;                                           # y-coords of input
    my $deg = shift;                   # largest poly degree to include in basis

    my $n = scalar @{$x};
    my @p = ();

    for (my $i = 0; $i < $n; $i++) {
        my @tmp = ();
        push @p, \@tmp;
    }

    if ($deg < -1 || $deg > 4) {
        die "\nUse a polynomial degree from -1 (no poly) up to 4, please.\n";
    }

    for (my $i = 0; $i < $n; $i++) {
        if ($deg >= 0) {
            push @{$p[$i]}, 1;
        }
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
        my @z = linalg::zeros($nE);
        for (my $i = 0; $i < $nE; $i++) {
            $z[$i] = linalg::phs(@{$x}[$i], @{$y}[$i], $rbfPow);
        }
        return \@z;
    } else {
        die "\nBad input.  Two equal-length arrays or two scalars, please.\n";
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
    
    my $A = linalg::zeros($nRows, $nCols);
    for (my $i = 0; $i < $nRows; $i++) {
        for (my $j = 0; $j < $nCols; $j++) {
            @{@{$A}[$i]}[$j] = rbf::phs(@{$x}[$i] - @{$xc}[$j]
            , @{$y}[$i] - @{$yc}[$j], $rbfPow);
        }
    }
    return $A;
}

################################################################################

sub interp2 {
    # Interpolate (x,y,f) to (X,Y,F)

    my $rbfPow = shift;                                    # exponent of phs rbf
    my $deg = shift;                        # largest polynomial degree in basis
    my $x = shift;             # pointer to x-coords where you KNOW the function
    my $y = shift;             # pointer to y-coords where you KNOW the function
    my $f = shift;                # pointer to known values of function on nodes
    my $X = shift;             # pointer to x-coords where you WANT the function
    my $Y = shift;             # pointer to y-coords where you WANT the function
    my $xmc = shift;  # pointer to x-coords of centers of rectangular subdomains
    my $ymc = shift;  # pointer to y-coords of centers of rectangular subdomains
    my $ell = shift;                      # length of each rectangular subdomain
    my $w = shift;                         # width of each rectangular subdomain

    my $ELL = (3 * $ell);
    my $W = (3 * $w);
    my $numP = int (($deg + 1) * ($deg + 2) / 2);
    my $zp1 = linalg::zeros($numP);
    my $zp2 = linalg::zeros($numP, $numP);
    my $F = linalg::zeros(scalar @{$X});
    my $nSubdomains = scalar @{$xmc};

    for (my $i = 0; $i < $nSubdomains; $i++) {
        # Get all nodes in the rectangular subdomain or adjacent subdomains.
        my $ind = rbf::inrectangle($x, $y, @{$xmc}[$i], @{$ymc}[$i], $ELL, $W);
        if ((scalar @{$ind}) < int (1.5 * $numP)) {
            print ((scalar @{$ind}) . "\n");
            die "\nNot enough data for this polynomial degree.\n";
        }
        my $xind = linalg::get($x, $ind);
        my $yind = linalg::get($y, $ind);
        # Make the polynomial matrix.
        my $p = rbf::polymat($xind, $yind, $deg);
        # Make the rbf matrix.
        my $A = rbf::rbfmat($xind, $yind, $xind, $yind, $rbfPow);
        # Put them together to create the combined rbf-poly matrix.
        $A = linalg::hstack($A, $p);
        $p = linalg::hstack(linalg::transpose($p), $zp2);
        $A = linalg::vstack($A, $p);
        # Get function values and solve for coefficients.
        my $lam = linalg::get($f, $ind);
		push @{$lam}, @{$zp1};
        $lam = linalg::solve($A, $lam);
        # Find evaluation points in the rectangular subdomain.
        my $IND = rbf::inrectangle($X, $Y, @{$xmc}[$i], @{$ymc}[$i], $ell, $w);
        my $XIND = linalg::get($X, $IND);
        my $YIND = linalg::get($Y, $IND);
        # Get rbf-poly evaluation matrix.
        $A = rbf::rbfmat($XIND, $YIND, $xind, $yind, $rbfPow);
        $p = rbf::polymat($XIND, $YIND, $deg);
        # Evaluate the interpolant at the evaluation points in the subdomain.
        linalg::set($F, $IND, linalg::dot(linalg::hstack($A, $p), $lam));
    }
	return $F;
}

################################################################################

return 1;

