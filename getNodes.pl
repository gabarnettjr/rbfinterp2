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
use rbf2;

################################################################################

# Process the input, if there is any.

my ($coordsDir, $nx, $ny, $alp, $a, $b, $c, $d);

$coordsDir = "randomCoords";
$nx = 16;
$ny = 16;
$alp = .3;
$a = 0;
$b = 50000;
$c = 0;
$d = 50000;

if (scalar @ARGV) {
	$coordsDir = shift;
	if (! -d $coordsDir) {
		print STDERR "First input must be a coordinates directory.\n"; die;
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

my ($xx, $yy) = rbf2::jostle($nx, $ny, $alp, $a, $b, $c, $d);

io::saveArray("$coordsDir\\x.txt", $xx);
io::saveArray("$coordsDir\\y.txt", $yy);

