# windows

package io;

# Simple functions for loading or saving a single array as a *.txt file.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

################################################################################

sub loadArray {
    my $file = shift;

    unless (open FILE, "<$file") {
        print STDERR "Unable to open file for reading.\n"; die;
    }
    my @values = <FILE>;
    close FILE;
    return \@values;
}

################################################################################

sub saveArray {
    my $fileName = shift;
    my $values = shift;

    unless (open FILE, ">$fileName") {
        print STDERR "Unable to open file for writing.\n"; die;
    }
    foreach my $val (@{$values}) {
		printf FILE "%1.15e\n", $val;
    }
    close FILE;
}

################################################################################

return 1;

