# windows

package linalg;

# A package for performing linear algebra operations on matrices and arrays.
# IMPORTANT: The inputs and outputs for the functions below are always POINTERS
# to an array/matrix, not the array itself.  The subroutines that begin with the
# prefix "test_" are a good place to look for examples showing how to use this
# package.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

################################################################################

sub copy {
    # Copy a matrix or array.  This only happens by defaut for arrays in perl.

    my $A = shift;                            # pointer to what you want to copy

    my @B;
    if (ref @{$A}[0]) {
        @B = ();
        foreach my $row (@{$A}) {
            my @tmp = @{$row};
            push @B, \@tmp;
        }
    } else {
        @B = @{$A};
    }
    return \@B;
}

################################################################################

sub sum {
    # Calculate the sum of all elements of an array.

    my $x = shift;                                         # pointer to an array

    my $s = 0;
    foreach my $r (@{$x}) {
        $s += $r;
    }
    return $s;
}

################################################################################

sub avg {
    # Calculate the average value of all elements of an array.
    return linalg::sum($_[0]) / (scalar @{$_[0]});
}

################################################################################

sub min {
    # Find the minimum value of all elements in an array.

    my $x = shift;                                  # pointer to the input array

    # Start by setting $m equal to the first element, then replace it if you
    # find something smaller.
    my $ell = scalar @{$x};
    my $m = @{$x}[0];
    for (my $i = 1; $i < $ell; $i++) {
        if (@{$x}[$i] < $m) {
            $m = @{$x}[$i];
        }
    }
    return $m;
}

################################################################################

sub indmin {
    # Find the index where the minimum value of an array occurs.

    my $x = shift;                                  # pointer to the input array

    # Start by setting $ind to the index of the first element, then replace it
    # if you find something smaller.
    my $ell = scalar @{$x};
    my $ind = 0;
    for (my $i = 1; $i < $ell; $i++) {
        if (@{$x}[$i] < @{$x}[$ind]) {
            $ind = $i;
        }
    }
    return $ind;
}

################################################################################

sub max {
    # Find the maximum value in an array.

    my $x = shift;                                  # pointer to the input array

    # Start out by setting the max $M equal to the first element, but replace it
    # if you find something larger.
    my $ell = scalar @{$x};
    my $M = @{$x}[0];
    for (my $i = 1; $i < $ell; $i++) {
        if (@{$x}[$i] > $M) {
            $M = @{$x}[$i];
        }
    }
    return $M;
}

################################################################################

sub indmax {
    # Find the index where the maximum value of an array occurs.

    my $x = shift;                                  # pointer to the input array

    # Start by choosing the first index, but replace with a new value if you
    # find a larger element of the array.
    my $ell = scalar @{$x};
    my $ind = 0;
    for (my $i = 1; $i < $ell; $i++) {
        if (@{$x}[$i] > @{$x}[$ind]) {
            $ind = $i;
        }
    }
    return $ind;
}

################################################################################

sub test_min_max {
    my $x = [4, 9, 7, 5, 8, 2, 1, 6];
    linalg::printmat($x);
    print ("min = " . linalg::min($x) . "\n");
    print ("max = " . linalg::max($x) . "\n");
    print ("indmin = " . linalg::indmin($x) . "\n");
    print ("indmax = " . linalg::indmax($x) . "\n");
    print "\n";
}

################################################################################

sub absval {
    # Get the absolute value of an array or matrix (abs val of each element).

    my @x = @{(shift)};                                # copy of the input array

    if (! ref $x[0]) {
        # Flip the sign of negative elements in the array, making them positive.
        my $n = scalar @x;
        for (my $i = 0; $i < $n; $i++) {
            if ($x[$i] < 0) {
                $x[$i] = -$x[$i];
            }
        }
        return \@x;
    } else {
        # Get the absolute value of each row in the matrix.
        my $n = scalar @x;
        for (my $i = 0; $i < $n; $i++) {
            $x[$i] = linalg::absval($x[$i]);
        }
        return \@x;
    }
}

################################################################################

sub test_absval {
    # Get absval of an array.
    my $x = [1, 2, -3];
    my $y = linalg::absval($x);
    linalg::printmat($y);
    # Get absval of a matrix.
    my $z = [[-1,2,-3], [5,4,-2]];
    linalg::printmat(linalg::absval($z));
}

################################################################################

sub norm {
    # Calculate the norm of a simple array (vector).  Default is 2-norm.
    
    my $x = shift;                                   # array to find the norm of
    my $type = shift;                             # type of norm (1,2,...,"inf")
    
    if (! $type) {
        return linalg::norm($x, 2);
    } elsif ($type eq "inf") {
        return linalg::max(linalg::absval($x));
    } else {
        my $n = 0;
        foreach my $r (@{$x}) {
            $n += ($r**$type);
        }
        $n = ($n ** (1/$type));
        return $n;
    }
}

################################################################################

sub test_norm {
    my $x = [1,2,3];
    linalg::printmat($x);
    print ("1-norm   = " . linalg::norm($x, 1) . "\n");
    print ("2-norm   = " . linalg::norm($x) . "\n");
    print ("3-norm   = " . linalg::norm($x, 3) . "\n");
    print ("inf-norm = " . linalg::norm($x, "inf") . "\n");
    print "\n";
}

################################################################################

sub linspace {
    # Get an array that starts at $a, stops at $b, and contains a total of $n
    # equally spaced values.  If $n = 1, it just gives the avg of $a and $b.

    my $a = shift;                                               # left endpoint
    my $b = shift;                                              # right endpoint
    my $n = shift;                                      # total number of points

    # Check that the input is of the correct form.
    if (ref $a || ref $b || ref $n) {
        print STDERR "Inputs should all be scalar values, not references.\n"; die;
    } elsif ($n == 1) {
        my @x = ();
        push @x, (($a + $b) / 2);
        return \@x;
    }

    # Tolerance to be used when comparing real number values for equality.
    my $eps = 1e-10;
    
    # Constant distance between adjacent points in the array.
    my $dx = ($b - $a) / ($n - 1);

    # Push values onto the array, one at a time.
    my @x = ();
    my $x = $a;
    while (abs ($x - $b) > $eps) {
        push @x, $x;
        $x += $dx;
    }
    push @x, $b;

    # Check that you ended up with the right number of elements, then return @x.
    if ((scalar @x) != $n) {
        print ((scalar @x) . "\n");
        print "$n\n";
        print (abs ($x - $b) . "\n");
        print STDERR "Something went wrong.  These should be equal.\n"; die;
    }
    return \@x;
}

################################################################################

sub test_linspace {
    my $x = linalg::linspace(0, 1, 11);
    my $y = linalg::linspace(-10, 10, 21);
    linalg::printmat($x);
    linalg::printmat($y);
}

################################################################################

sub get {
    # Get @x at the indices contained in @ind.

    my $x = shift;                                             # array of values
    my $ind = shift;                                          # array of indices

    # Check input.
    if (! ref $x || ! ref $ind) {
        print STDERR "Both inputs should be references to arrays.\n"; die;
    } elsif (ref @{$x}[0]) {
        print STDERR "Only works for simple arrays, not matrices.\n"; die;
    }

    # Get values and return them in an array.
    my @xind = ();
    foreach my $i (@{$ind}) {
        push @xind, @{$x}[$i];
    }
    return \@xind;
}

################################################################################

sub set {
    # Set the values of @x at indices @ind equal to @y.

    my $x = shift;                                  # pointer to array of values
    my $ind = shift;                               # pointer to array of indices
    my $y = shift;                            # pointer to other array of values

    # Check that the input is okay.
    if ((scalar @{$ind}) != (scalar @{$y})) {
        print STDERR "\@ind and \@y must be the same length.\n"; die;
    } elsif (! ref $x || ! ref $ind || ! ref $y) {
        print STDERR "All inputs should be references to arrays.\n"; die;
    } elsif (ref @{$x}[0] || ref @{$ind}[0] || ref @{$y}[0]) {
        print STDERR "Inputs should be simple arrays, not matrices.\n"; die;
    }
    
    # Set the array elements to the values.
    my $j = 0;
    foreach my $i (@{$ind}) {
        @{$x}[$i] = @{$y}[$j];
        $j++;
    }
}

################################################################################

sub test_get_set {
    my $x = [0, 1, 2, 3, 4];
    my $ind = [2, 4];
    my $y = linalg::get($x, $ind);
    my $z = [-12, 99];
    linalg::printmat($x);
    linalg::set($x, [1,3], $z);
    linalg::printmat($x);
    linalg::printmat($y);
}

################################################################################

sub meshgrid {
    # Create matrices @xx and @yy based on arrays @x and @y.  @x is repeated and
	# stacked vertically, while @y is repeated and stacked horizontally.  In the
	# end, the output matrices xx and yy should have numRows = len(y) and
	# numCols = len(x).

    my $x = shift;                                # pointer to array of x-coords
    my $y = shift;                                # pointer to array of y-coords

    my $nRows = scalar @{$y};
    my $nCols = scalar @{$x};
    
    # Check that the input is appropriate.
    if (ref @{$x}[0] || ref @{$y}[0]) {
        print STDERR "\@x and \@y must both be simple arrays.\n"; die;
    }

    # Repeat @x vertically, @y horizontally, and return pointers to matrices.
    my @xx = ();
    my @yy = ();
    for (my $i = 0; $i < $nRows; $i++) {
        # append new row of @xx.
        my @tmp1 = @{$x};
        push @xx, \@tmp1;
        # append new row of @yy.
        my @tmp2 = ();
        for (my $j = 0; $j < $nCols; $j++) {
            push @tmp2, @{$y}[$i];
        }
        push @yy, \@tmp2;
    }
    return (\@xx, \@yy);
}

################################################################################

sub test_meshgrid {
    my $x = [1,2,3,4];
    my $y = [5,6,7];
    my ($xx, $yy) = linalg::meshgrid($x, $y);
    linalg::printmat($xx);
    linalg::printmat($yy);
}

################################################################################

sub flatten {
    # Flatten the matrix @A to the array @x.

    my $A = shift;                 # matrix (array of pointers to arrays (rows))
    
    # Make sure that the input is a matrix.
    if (! ref @{$A}[0]) {
        print STDERR "Only a matrix can be flattened, not an array.\n"; die;
    }

    # Push each row of matrix @A onto array @x, then return @x.
    my $nRows = scalar @{$A};
    my @x = ();
    for (my $i = 0; $i < $nRows; $i++) {
        push @x, @{@{$A}[$i]};
    }
    return \@x;
}

################################################################################

sub scalarmul {
    # Multiply an array or matrix by a scalar.

    my $x;                                                               # array
    my $r;                                                              # scalar

    # Determine the ordering of the inputs.
    if (ref $_[0] && ! ref $_[1]) {
        $x = linalg::copy(shift);
        $r = shift;
    } elsif (ref $_[1] && ! ref $_[0]) {
        $r = shift;
        $x = linalg::copy(shift);
    } else {
        print STDERR "One input should be an array, the other a scalar.\n";  die;
    }
    
    # Multiply each element of @x by $r, then return a pointer to @x.
    if (ref @{$x}[0]) {
        for (my $i = 0; $i < (scalar @{$x}); $i++) {
            @{$x}[$i] = linalg::scalarmul($r, @{$x}[$i]);
        }
        return $x;
    } else {
        for (my $i = 0; $i < (scalar @{$x}); $i++) {
            @{$x}[$i] *= $r;
        }
        return $x;
    }
}

################################################################################

sub scalaradd {
    # Add an array or matrix to a scalar.

    my $x;                                                               # array
    my $r;                                                              # scalar

    # Determine the ordering of the inputs.
    if (ref $_[0] && ! ref $_[1]) {
        $x = linalg::copy(shift);
        $r = shift;
    } elsif (ref $_[1] && ! ref $_[0]) {
        $r = shift;
        $x = linalg::copy(shift);
    } else {
        print STDERR "One input should be an array, the other a scalar.\n";  die;
    }
    
    # Add each element of @x to $r, then return a pointer to @x.
    if (ref @{$x}[0]) {
        for (my $i = 0; $i < (scalar @{$x}); $i++) {
            @{$x}[$i] = linalg::scalaradd($r, @{$x}[$i]);
        }
        return $x;
    } else {
        for (my $i = 0; $i < (scalar @{$x}); $i++) {
            @{$x}[$i] += $r;
        }
        return $x;
    }
}

################################################################################

sub test_scalaradd_scalarmul {
    my $A = [[1,2,3], [4,5,6]];
    linalg::printmat($A);
    linalg::printmat(linalg::scalaradd($A, -1));
    linalg::printmat(linalg::scalarmul($A, -1));
}

################################################################################

sub hstack {
    # Stack two matrices horizontally.

    my $a = shift;                                     # pointer to first matrix
    my $b = shift;                                    # pointer to second matrix
    
    # Make sure the operation is possible, given inputs $a and $b.
    my $numRows = scalar @{$a};
    if (! ref @{$a}[0] || ! ref @{$b}[0]) {
        print STDERR "This only works if both inputs are matrices.\n"; die;
    } elsif ($numRows != (scalar @{$b})) {
        print STDERR "Matrices must have same number of rows to stack horizontally.\n"; die;
    }
    
    # Make and return the new stacked matrix.
    my @c = ();
    for (my $i = 0; $i < $numRows; $i++) {
        my @tmp = @{@{$a}[$i]};
        push @tmp, @{@{$b}[$i]};
        push @c, \@tmp;
    }
    return \@c;
}

################################################################################

sub vstack {
    # Stack two matrices vertically.

    my $a = linalg::copy(shift);                       # pointer to first matrix
    my $b = linalg::copy(shift);                      # pointer to second matrix
    
    # Check that the two inputs are valid.
    if (! ref @{$a}[0] || ! ref @{$b}[0]) {
        print STDERR "This only works if both inputs are matrices.\n"; die;
    }
    
    # Make and return the new stacked matrix.
    push @{$a}, @{$b};
    return $a;
}

################################################################################

sub test_hstack_vstack {
    # Make a matrix @A.
    my $A = [[1,2,3], [4,5,6]];
    # Make a copy @B of @A, and change one element.
    my $B = linalg::copy($A);
    @{@{$B}[1]}[1] = 23;
    # Stack @A next to @B and call it @C.
    my $C = linalg::hstack($A, $B);
    # Stack @A on top of @B and call it @D.
    my $D = linalg::vstack($A, $B);
    # Print all three matrices.
    linalg::printmat($A);
    linalg::printmat($B);
    linalg::printmat($C);
    linalg::printmat($D);
}

################################################################################

sub zeros {
    # Create an array or matrix of zeros.

    my $nRows = shift;                                  # desired number of rows
    my $nCols = shift;                               # desired number of columns

    # Check if it is an array or matrix, based on the number of inputs,
    # then create the array/matrix of zeros, @z.
    my @z = ();
    if ($nRows && ! $nCols) {
        for (my $i = 0; $i < $nRows; $i++) {
            push @z, 0;
        }
    } elsif ($nRows && $nCols) {
        for (my $i = 0; $i < $nRows; $i++) {
            my @tmp = ();
            for (my $j = 0; $j < $nCols; $j++) {
                push @tmp, 0;
            }
            push @z, \@tmp;
        }
    } elsif ($nRows == 0 || $nCols == 0) {
        # Do nothing.
    } else {
        print STDERR "Bad input.  Please try again.\n"; die;
    }
    return \@z;
}

################################################################################

sub test_zeros {
	my $z0 = linalg::zeros(0);
    my $z1 = linalg::zeros(5);
    my $z2 = linalg::zeros(3, 6);
	linalg::printmat($z0);
    linalg::printmat($z1);
    linalg::printmat($z2);
}

################################################################################

sub eye {
    # Get a square identity matrix, which is all zeros except ones on diagonal.
    
    my $m = shift;               # number of rows and columns of identity matrix
    
    my $I = linalg::zeros($m, $m);
    for (my $j = 0; $j < $m; $j++) {
        @{@{$I}[$j]}[$j] = 1;
    }
    return $I;
}

################################################################################

sub test_eye {
    my $A = [[1,2,3], [4,5,6], [7,8,9]];
    my $I = linalg::eye(scalar @{$A});
    linalg::printmat(linalg::add($A, linalg::scalarmul(-1, $I)));
}

################################################################################

sub transpose {
    # Transpose a matrix (swap rows and columns).

    my $A = shift;                                 # the matrix to be transposed

    # Number of rows and columns in the original matrix.
    my $m = scalar @{$A};
    my $n = scalar @{@{$A}[0]};
    
    # Assign elements to the transpose matrix, then return it.
    my $At = linalg::zeros($n, $m);
    for (my $i = 0; $i < $n; $i++) {
        for (my $j = 0; $j < $m; $j++) {
            @{@{$At}[$i]}[$j] = @{@{$A}[$j]}[$i];
        }
    }
    return $At;
}

################################################################################

sub test_transpose {
    # Make a matrix @A.
    my $A = [[1,2,3], [4,5,6]];
    # Make the transpose of @A.
    my $B = linalg::transpose($A);
    # Print both matrices.
    linalg::printmat($A);
    linalg::printmat($B);
}

################################################################################

sub add {
    # Add two matrices or arrays together.

    my $x = shift;                                 # pointer to the first matrix
    my $y = shift;                                # pointer to the second matrix

    if (! ref @{$x}[0] && ! ref @{$y}[0]) {
        if ((scalar @{$x}) != (scalar @{$y})) {
            print STDERR "Arrays must be the same length to be added together.\n"; die;
        }
        my $ne = scalar @{$x};
        my @z = @{$x};
        for (my $i = 0; $i < $ne; $i++) {
            $z[$i] += @{$y}[$i];
        }
        return \@z;
    } elsif (ref @{$x}[0] && ref @{$y}[0]) {
        if ((scalar @{$x}) != (scalar @{$y})) {
            print STDERR "Matrices must have same number of rows to be added together.\n"; die;
        } elsif ((scalar @{@{$x}[0]}) != (scalar @{@{$y}[0]})) {
            print STDERR "Matrices must have same number of columns to be added together.\n"; die;
        }
        my $nRows = scalar @{$x};
        my @z = @{$x};
        for (my $i = 0; $i < $nRows; $i++) {
            $z[$i] = linalg::add(@{$x}[$i], @{$y}[$i]);
        }
        return \@z;
    } else {
        print STDERR "Inputs must either both be arrays or both be matrices.\n"; die;
    }
}

################################################################################

sub test_add {
	# Make a 2x3 matrix @A.
    my $A = [[1,2,3], [4,5,6]];
	# Make another 2x3 matrix @B.
    my $B = [[2,3,4], [3,4,5]];
	# See what @A and @B look like.
	linalg::printmat($A);
	linalg::printmat($B);
	# Find the sum of @A and @B.
	my $C = linalg::add($A, $B);
	linalg::printmat($C);
    # Add first row of @A to first row of @B.
    my $c = linalg::add(@{$A}[0], @{$B}[0]);
    linalg::printmat($c);
}

################################################################################

sub dot {
    # Compute the dot product of arrays or matrices.  Note that simple arrays
	# are considered to be column vectors.  To calculate the dot product of a
	# row vector x and a matrix A, instead get the dot product of A^t with x^t.

    my $x = shift;                                       # first array or matrix
    my $y = shift;                                      # second array or matrix

    if (! ref @{$x}[0] && ! ref @{$y}[0]) {
        if ((scalar @{$x}) != (scalar @{$y})) {
			print (scalar @{$x});
			print "\n";
			print (scalar @{$y});
			print "\n";
            print STDERR "Arrays must be the same length to be dotted.\n"; die;
        }
        my $ne = scalar @{$x};
        my $dot = 0;
        for (my $i = 0; $i < $ne; $i++) {
            $dot += (@{$x}[$i] * @{$y}[$i]);
        }
        return $dot;
    } elsif (ref @{$x}[0] && ! ref @{$y}[0]) {
        my $colsX = scalar @{@{$x}[0]};
        my $rowsY = scalar @{$y};
        if ($colsX != $rowsY) {
            print "\n\$colsX = $colsX\n";
            print "\n\$rowsY = $rowsY\n";
            print STDERR "nCols of LHS matrix must equal length of RHS array..\n"; die;
        }
		my $nRows = scalar @{$x};
        my $z = linalg::zeros($nRows);
        for (my $i = 0; $i < $nRows; $i++) {
            @{$z}[$i] = linalg::dot(@{$x}[$i], $y);
        }
        return $z;
    } elsif (ref @{$x}[0] && ref @{$y}[0]) {
        my $colsX = scalar @{@{$x}[0]};
        my $rowsY = scalar @{$y};
        if ($colsX != $rowsY) {
            print "\n\$colsX = $colsX\n";
            print "\n\$rowsY = $rowsY\n";
            print STDERR "nCols of first matrix must equal nRows of second matrix.\n"; die;
        }
        my $yT = linalg::transpose($y);
        my $nRows = scalar @{$x};
		my $nCols = scalar @{$yT};
        my $z = linalg::zeros($nRows, $nCols);
        for (my $i = 0; $i < $nRows; $i++) {
            for (my $j = 0; $j < $nCols; $j++) {
                @{@{$z}[$i]}[$j] = linalg::dot(@{$x}[$i], @{$yT}[$j]);
            }
        }
        return $z;
    } else {
        print STDERR "Inputs must matrix-matrix, array-array, or matrix-array.\n"; die;
    }
}

################################################################################

sub test_dot {
	# Make a 2x3 matrix @A.
    my $A = [[1,2,3], [4,5,6]];
	# Make a 3x4 matrix @B.
    my $B = [[2,3,4,5], [3,4,5,6], [4,5,6,7]];
	# See what @A and @B look like.
	linalg::printmat($A);
	linalg::printmat($B);
	# Find the product of @A and @B.
	my $C = linalg::dot($A, $B);
	linalg::printmat($C);
}

################################################################################

sub printmat {
    # Print an array or matrix to standard output.

    my $a = shift;                                 # pointer to the input matrix

    # Print each row, followed by a new line character.
    foreach my $row (@{$a}) {
        if (ref @{$a}[0]) {
            foreach my $elem (@{$row}) {
                printf "%10.7f ", $elem;
            }
            print "\n";
        } else {
            printf "%8.5f ", $row;
        }
    }
    print "\n";
    if (! ref @{$a}[0]) {
        print "\n";
    }
}

################################################################################

sub swaprows {
    # Swap rows $i and $j of the matrix @A.

    my $A = shift;                                    # pointer to the matrix @A
    my $i = shift;                                            # index of one row
    my $j = shift;                                          # index of other row
    
    # Swap pointers to the rows.
    my $tmp = @{$A}[$i];
    @{$A}[$i] = @{$A}[$j];
    @{$A}[$j] = $tmp;
}

################################################################################

sub solve {
    # Solve the linear system @A*@x=@b, with matrix @A and array @b known.

    my $A = linalg::copy(shift);                          # the square matrix @A
    my $b = linalg::copy(shift);                                  # the array @b

    my $nRows = scalar @{$A};
    my $nCols = scalar @{@{$A}[0]};
    if ($nRows != $nCols || $nRows != (scalar @{$b})) {
        print STDERR "A should be square.  Rows(A) should equal len(b).\n"; die;
    }

    # Apply row operations to transform A to upper triangular form.
    for (my $j = 0; $j < $nCols - 1; $j++) {
        # Find largest element in leftmost column, and make this the pivot row.
        # Apparently, this is important, because when I was not doing this, I
        # was getting noticeably different answers from python when the number
        # of subdomains was small.
        my $indmax = $j;
        for (my $i = $j+1; $i < $nRows; $i++) {
            if ((abs @{@{$A}[$i]}[$j]) > (abs @{@{$A}[$indmax]}[$j])) {
                $indmax = $i;
            }
        }
        # Swap row $indmax and row $j.
        linalg::swaprows($A, $indmax, $j);
        my $tmp = @{$b}[$j];
        @{$b}[$j] = @{$b}[$indmax];
        @{$b}[$indmax] = $tmp;
        # Zero out the $jth column.
        for (my $i = $j + 1; $i < $nRows; $i++) {
            my $factor = -@{@{$A}[$i]}[$j] / @{@{$A}[$j]}[$j];
            for (my $k = $j; $k < $nCols; $k++) {
                @{@{$A}[$i]}[$k] += ($factor * @{@{$A}[$j]}[$k]);
            }
            @{$b}[$i] += ($factor * @{$b}[$j]);
        }
    }

    # Use back substitution to finish solving for @x.
    my $x = $b;
    @{$x}[$nRows-1] = @{$b}[$nRows-1] / @{@{$A}[$nRows-1]}[$nRows-1];
    my $k = $nRows;
    for (my $i = $nRows - 2; $i >= 0; $i--) {
        $k--;
        my $dot = 0;
        for (my $j = $k; $j < $nRows; $j++) {
            $dot += (@{@{$A}[$i]}[$j] * @{$x}[$j]);
        }
        @{$x}[$i] = (@{$b}[$i] - $dot) / @{@{$A}[$i]}[$i];
    }
    return $x;
}

################################################################################

sub test_solve {
    # Make up a matrix @A.
    my $A = [[1,2,1], [4,5,2], [1,8,3]];
    # Make up a vector @b.
    my $b = [9, 3, 2];
    # Solve the system.
    my $x = linalg::solve($A, $b);
    # Test if it worked.
    print "\@A =\n";
    linalg::printmat($A);
    print "\@x = ";
    linalg::printmat($x);
    my $prod = dot($A, $x);
    print "\@A\*\@x = ";
    linalg::printmat($prod);
    print "\@b    = ";
    linalg::printmat($b);
}

################################################################################

return 1;

