# windows

# Sometimes you already have some nodes that you are working with, (x, y), and
# what you need are some evaluation points (xe, ye) to interpolate to.  This
# script will read the files x.txt and y.txt, and use their values to generate
# a nice collection of corresponding evaluation points for testing.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

use lib ".";
use linalg;
use io;
use rbf2;

################################################################################

# Process the input, if there is any.

my ($coordsDir, $nx, $ny, $alp);

$coordsDir = "randomCoords";
$nx = 32;
$ny = 32;
$alp = 0;

if (scalar @ARGV) {
    $coordsDir = shift;
}
if (scalar @ARGV) {
    $nx = shift;
}
if (scalar @ARGV) {
    $ny = shift;
}
if (scalar @ARGV) {
	$alp = shift;
}

################################################################################

my $x = io::loadArray("$coordsDir\\x.txt");
my $y = io::loadArray("$coordsDir\\y.txt");

my $a = linalg::min($x);
my $b = linalg::max($x);
my $c = linalg::min($y);
my $d = linalg::max($y);

my ($xe, $ye) = rbf2::jostle($nx, $ny, $alp, $a, $b, $c, $d);

# my $xe = linalg::linspace($a, $b, $nx);
# my $ye = linalg::linspace($c, $d, $ny);

# ($xe, $ye) = linalg::meshgrid($xe, $ye);
# $xe = linalg::flatten($xe);
# $ye = linalg::flatten($ye);

# my $eps = .001;
# $alp = ($b - $a) / ($ne - 1) * $alp;

# $ne = scalar @{$xe};

# # Jostle the Cartesian nodes (except corners) by a fraction of the node spacing.
# for (my $i = 0; $i < $ne; $i++) {
    # my $x = @{$xe}[$i];
    # my $y = @{$ye}[$i];
    # if (($x > ($a + $eps)) && ($x < ($b - $eps)) && ($y > ($c + $eps)) && ($y < ($d - $eps))) {
        # my $r = -$alp + 2 * $alp * rand;
        # @{$xe}[$i] += $r;
        # $r = -$alp + 2 * $alp * rand;
        # @{$ye}[$i] += $r;
    # } elsif ((($x < ($a + $eps)) || ($x > ($b - $eps))) && ($y > ($c + $eps)) && ($y < ($d - $eps))) {
		# my $r = -$alp + 2 * $alp * rand;
		# @{$ye}[$i] += $r;
	# } elsif ((($y < ($c + $eps)) || ($y > ($d - $eps))) && ($x > ($a + $eps)) && ($x < ($b - $eps))) {
		# my $r = -$alp + 2 * $alp * rand;
		# @{$xe}[$i] += $r;
	# }
# }

io::saveArray("$coordsDir\\xe.txt", $xe);
io::saveArray("$coordsDir\\ye.txt", $ye);

