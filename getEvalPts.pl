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

my ($coordsDir, $nx, $ny, $alp, $a, $b, $c, $d);

$coordsDir = "randomCoords";
$nx = 32;
$ny = 32;
$alp = 0;
$a = $b = $c = $d = "";

if (scalar @ARGV) {
	$coordsDir = shift;
	if (! -d $coordsDir) {
		print STDERR "First input must be a coordinates directory.\n"; die;
	} elsif ((! -e "$coordsDir\\x.txt") || (! -e "$coordsDir\\y.txt")) {
		print STDERR "Coordinates directory must contain nodes.\n"; die;
	}
}
if (scalar @ARGV) { $nx = shift; }
if (scalar @ARGV) { $ny = shift; }
if (scalar @ARGV) { $alp = shift; }
if (scalar @ARGV) { $a = shift; }
if (scalar @ARGV) { $b = shift; }
if (scalar @ARGV) { $c = shift; }
if (scalar @ARGV) { $d = shift; }

################################################################################

if ($a eq "" || $b eq "" || $c eq "" || $d eq "") {
	# Use nodes to get boundaries of eval pts.
	my $x = io::loadArray("$coordsDir\\x.txt");
	my $y = io::loadArray("$coordsDir\\y.txt");
	$a = linalg::min($x);
	$b = linalg::max($x);
	$c = linalg::min($y);
	$d = linalg::max($y);
}

my ($xe, $ye) = rbf2::jostle($nx, $ny, $alp, $a, $b, $c, $d);

io::saveArray("$coordsDir\\xe.txt", $xe);
io::saveArray("$coordsDir\\ye.txt", $ye);

