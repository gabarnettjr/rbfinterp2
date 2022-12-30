# windows

use strict;
use warnings;

use lib ".";
use linalg;
use io;

my $coordsDir = "myCoords";
my $ne = 51;
my $alp = 0;

################################################################################

my $x = io::loadArray("$coordsDir\\x.txt");
my $y = io::loadArray("$coordsDir\\y.txt");

my $a = linalg::min($x);
my $b = linalg::max($x);
my $c = linalg::min($y);
my $d = linalg::max($y);

my $xe = linalg::linspace($a, $b, $ne);
my $ye = linalg::linspace($c, $d, $ne);

($xe, $ye) = linalg::meshgrid($xe, $ye);
$xe = linalg::flatten($xe);
$ye = linalg::flatten($ye);

my $eps = .001;
$alp = ($b - $a) / ($ne - 1) * $alp;

$ne = scalar @{$xe};

for (my $i = 0; $i < $ne; $i++) {
    my $x = @{$xe}[$i];
    my $y = @{$ye}[$i];
    if (($x > ($a + $eps)) && ($x < ($b - $eps)) && ($y > ($c + $eps)) && ($y < ($d - $eps))) {
        my $r = -$alp + 2 * $alp * rand;
        @{$xe}[$i] += $r;
        $r = -$alp + 2 * $alp * rand;
        @{$ye}[$i] += $r;
    }
}

io::saveArray("$coordsDir\\xe.txt", $xe);
io::saveArray("$coordsDir\\ye.txt", $ye);

