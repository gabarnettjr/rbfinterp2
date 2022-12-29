# windows

package io;

use strict;
use warnings;

################################################################################

sub loadArray {
    my $file = shift;

    unless (open FILE, "<$file") {
        die "\nUnable to open file for reading.\n";
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
        die "\nUnable to open file for writing.\n";
    }
    foreach my $val (@{$values}) {
        printf FILE "%-17.14f\n", $val;
    }
    close FILE;
}

################################################################################

return 1;

