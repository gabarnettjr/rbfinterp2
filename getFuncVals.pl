# windows

# If you already have defined nodes (x, y) and evaluation points (xe, ye), then
# this script will use a smooth function to define data at these points.
# The script generates f.txt to go with (x, y), and fe.txt to go with (xe, ye).
# In order for this to work, you need to already have x.txt, y.txt, xe.txt, and
# ye.txt saved in the folder at location $dataDir.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

use Math::Trig;

use lib ".";
use io;
use linalg;

################################################################################

my ($dataDir, $f);

$dataDir = "randomCoords\\smoothData";
$f = "f";

if (scalar @ARGV) {
	$dataDir = shift;
	if (! -d $dataDir) {
		print STDERR "First input must be a data directory.\n"; die;
	} elsif (! -e "$dataDir\\f.txt") {
		print STDERR "Data directory must contain function values.\n"; die;
	}
}
if (scalar @ARGV) { $f = shift; }

################################################################################

my $x = io::loadArray("$dataDir\\..\\x.txt");
my $y = io::loadArray("$dataDir\\..\\y.txt");

my $xe = io::loadArray("$dataDir\\..\\xe.txt");
my $ye = io::loadArray("$dataDir\\..\\ye.txt");

################################################################################

sub f {
	my ($x, $y) = @_;
	my $z;
	if ($f eq "f") {
		$z = (cos 2*pi * $x / 25000) * (sin 2*pi * $y / 25000);
	} elsif ($f eq "g") {
		$z = (cos 2*pi * $x / 25000) * (cos 2*pi * $y / 25000)
	       + (sin 2*pi * $x / 25000) * (sin 2*pi * $y / 25000);
	} elsif ($f eq "h") {
		$z = (cos 2*pi * $x / 50000) * (cos 2*pi * $y / 50000)
	       + (sin 2*pi * $x / 50000) * (sin 2*pi * $y / 50000);
	} elsif ($f eq "i") {
		$z = (exp -(1/8000)**2 * (($x - 10000)**2 + ($y - 35000)**2))
		   + (exp -(1/8000)**2 * (($x - 32000)**2 + ($y - 5000)**2))
		   + (exp -(1/5000)**2 * (($x - 41000)**2 + ($y - 38000)**2));
	} else {
		print STDERR "Invalid choice for \$f.  Choose f, g, or h.\n"; die;
	}
	return $z;
}

sub eff {
	my ($x, $y) = @_;
	my $z = linalg::zeros(scalar @{$x});
	for (my $i = 0; $i < (scalar @{$x}); $i++) {
		@{$z}[$i] = f(@{$x}[$i], @{$y}[$i]);
	}
	return $z;
}

################################################################################

io::saveArray("$dataDir\\f.txt", eff($x, $y));

io::saveArray("$dataDir\\fe.txt", eff($xe, $ye));

