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

my ($dataDir, $ftype);

$dataDir = "randomCoords\\smoothData";
$ftype = "peaks and valleys";

if (scalar @ARGV) {
	$dataDir = shift;
	if (! -d $dataDir) {
		print STDERR "First input must be a data directory.\n"; die;
	}
}
if (scalar @ARGV) { $ftype = (lc shift); }

################################################################################

my $x = io::loadArray("$dataDir\\..\\x.txt");
my $y = io::loadArray("$dataDir\\..\\y.txt");
my $xe = io::loadArray("$dataDir\\..\\xe.txt");
my $ye = io::loadArray("$dataDir\\..\\ye.txt");

my @ab = (@{$x}, @{$xe});
my $a = linalg::min(\@ab);
my $b = linalg::max(\@ab);
my @cd = (@{$y}, @{$ye});
my $c = linalg::min(\@cd);
my $d = linalg::max(\@cd);

################################################################################

sub f {
    my ($ftype, $x, $y, $a, $b, $c, $d) = @_;
    my $w = ($b - $a);
    my $ell = ($d - $c);
    my $s = ($w + $ell) / 2;
    my $z;
    if ($ftype eq "peaks and valleys" || $ftype eq "1") {
        $z = (cos 2*pi * $x / ($w/2)) * (sin 2*pi * $y / ($ell/2));
    } elsif ($ftype eq "narrow stripes" || $ftype eq "2") {
        $z = (cos 2*pi * $x / ($w/2)) * (cos 2*pi * $y / ($ell/2))
           + (sin 2*pi * $x / ($w/2)) * (sin 2*pi * $y / ($ell/2));
    } elsif ($ftype eq "wide stripes" || $ftype eq "3") {
        $z = (cos 2*pi * $x / $w) * (cos 2*pi * $y / $ell)
           + (sin 2*pi * $x / $w) * (sin 2*pi * $y / $ell);
    } elsif ($ftype eq "bells" || $ftype eq "4") {
        $z = (exp -(1/(.16*$s))**2 * (($x - ($a + .20*$w))**2 + ($y - ($c + .70*$ell))**2))
           + (exp -(1/(.16*$s))**2 * (($x - ($a + .64*$w))**2 + ($y - ($c + .10*$ell))**2))
           + (exp -(1/(.10*$s))**2 * (($x - ($a + .82*$w))**2 + ($y - ($c + .76*$ell))**2));
    } else {
        my $s = "Invalid choice for \$ftype.  Try again using one of these:\n";
        $s .= "(1) \"peaks and valleys\"\n";
        $s .= "(2) \"narrow stripes\"\n";
        $s .= "(3) \"wide stripes\"\n";
        $s .= "(4) \"bells\"\n";
        print STDERR $s; die;
    }
    return $z;
}

sub eff {
    my ($ftype, $x, $y, $a, $b, $c, $d) = @_;
    my $z = linalg::zeros(scalar @{$x});
    for (my $i = 0; $i < (scalar @{$x}); $i++) {
        @{$z}[$i] = f($ftype, @{$x}[$i], @{$y}[$i], $a, $b, $c, $d);
    }
    return $z;
}

################################################################################

io::saveArray("$dataDir\\f.txt", eff($ftype, $x, $y, $a, $b, $c, $d));

io::saveArray("$dataDir\\fe.txt", eff($ftype, $xe, $ye, $a, $b, $c, $d));

