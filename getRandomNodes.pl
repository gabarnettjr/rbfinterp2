# windows

# Starting from a Cartesian grid, this script will jostle the nodes by some
# chosen amount to create "random" nodes for testing.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

use lib ".";
use linalg;
use io;

################################################################################

# Process the input, if there is any.

my ($coordsDir, $alp, $nx, $ny, $a, $b, $c, $d);

$coordsDir = "randomCoords";
$alp = .25;
$nx = 16;
$ny = 16;
$a = 0;
$b = 50000;
$c = 0;
$d = 50000;

if (scalar @ARGV) {
    $coordsDir = shift;
}
if (scalar @ARGV) {
    $alp = shift;
}
if (scalar @ARGV) {
    $nx = shift;
}
if (scalar @ARGV) {
	$ny = shift;
}
if (scalar @ARGV) {
	$a = shift;
}
if (scalar @ARGV) {
	$b = shift;
}
if (scalar @ARGV) {
	$c = shift;
}
if (scalar @ARGV) {
	$d = shift;
}

################################################################################

my $xx = linalg::linspace($a, $b, $nx);
my $yy = linalg::linspace($c, $d, $ny);

($xx, $yy) = linalg::meshgrid($xx, $yy);
$xx = linalg::flatten($xx);
$yy = linalg::flatten($yy);

my $eps = .001;
$alp = (($b - $a) / ($nx - 1) * $alp + ($d - $c) / ($ny - 1) * $alp) / 2;

my $ne = scalar @{$xx};

# Jostle the Cartesian nodes (except corners) by a fraction of the node spacing.
for (my $i = 0; $i < $ne; $i++) {
    my $x = @{$xx}[$i];
    my $y = @{$yy}[$i];
    if (($x > ($a + $eps)) && ($x < ($b - $eps)) && ($y > ($c + $eps))
	&& ($y < ($d - $eps))) {
		# Interior
        my $r = -$alp + 2 * $alp * rand;
        @{$xx}[$i] += $r;
        $r = -$alp + 2 * $alp * rand;
        @{$yy}[$i] += $r;
    } elsif ((($x < ($a + $eps)) || ($x > ($b - $eps))) && ($y > ($c + $eps))
	&& ($y < ($d - $eps))) {
		# Left and right boundaries, but no corners.
		my $r = -$alp + 2 * $alp * rand;
		@{$yy}[$i] += $r;
	} elsif ((($y < ($c + $eps)) || ($y > ($d - $eps))) && ($x > ($a + $eps))
	&& ($x < ($b - $eps))) {
		# Top and bottom boundaries, but no corners.
		my $r = -$alp + 2 * $alp * rand;
		@{$xx}[$i] += $r;
	}
}

io::saveArray("$coordsDir\\x.txt", $xx);
io::saveArray("$coordsDir\\y.txt", $yy);

