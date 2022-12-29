# windows

use strict;
use warnings;

use lib ".";
use io;
use linalg;
use rbf;

my ($dataDir, $checkError);
if ((scalar @ARGV) == 0) {
    $dataDir = "cdCoords\\smoothData";
    $checkError = 1;
} else {
    $dataDir = (shift @ARGV);
    if (scalar @ARGV) {
        $checkError = (shift @ARGV);
    } else {
        $checkError = 0;
    }
}
my $rbfPow = 3;
my $deg = 1;
my $nSubd = 8;
my $mSubd = $nSubd;

################################################################################

unless (-e "$dataDir\\f.txt") {
    # Use python to generate some smooth data to interpolate.  If you want to
    # generate new data using a new function, then delete f.txt, re-define the
    # function in generateSmoothData.py, and then re-run rbfinterp2.pl.
    system "python generateSmoothData.py $dataDir";
}

# Start the timer for loading/computing/saving in perl.
my $computeTime = time;

# Load the data at the nodes (where you know the function).
my $x = io::loadArray("$dataDir\\..\\x.txt");
my $y = io::loadArray("$dataDir\\..\\y.txt");
my $f = io::loadArray("$dataDir\\f.txt");

# Load the data at the evaluation points (where you WANT to know the function).
my $xe = io::loadArray("$dataDir\\..\\xe.txt");
my $ye = io::loadArray("$dataDir\\..\\ye.txt");

################################################################################

# Translate and scale the coordinates for a well-conditioned problem.

# Shift so that (0,0) is the center of the computational domain.
my $xavg = linalg::avg($x);
my $yavg = linalg::avg($y);
$x = linalg::scalaradd($x, -$xavg);
$y = linalg::scalaradd($y, -$yavg);
$xe = linalg::scalaradd($xe, -$xavg);
$ye = linalg::scalaradd($ye, -$yavg);

# Constant re-scaling of the x and y coordinates (needed for high order poly).
my $alp = (linalg::max(linalg::absval($x)) + linalg::max(linalg::absval($y))) / 2;
$x = linalg::scalarmul($x, (1/$alp));
$y = linalg::scalarmul($y, (1/$alp));
$xe = linalg::scalarmul($xe, (1/$alp));
$ye = linalg::scalarmul($ye, (1/$alp));

# The rectangular computational domain, [a,b] x [c,d].
my @ab = (@{$x}, @{$xe});
my $a = linalg::min(\@ab);
my $b = linalg::max(\@ab);
my @cd = (@{$y}, @{$ye});
my $c = linalg::min(\@cd);
my $d = linalg::max(\@cd);

# print ($a . " " . $b . " " . $c . " " . $d . "\n");

################################################################################

# Find the center (xmc,ymc) and dimensions of each rectangular subdomain.

my $eps = 0.001;

my $dx = ($b - $a + 2*$eps) / $nSubd;
my $xmc = linalg::linspace($a - $eps + $dx/2, $b + $eps - $dx/2, $nSubd);

my $dy = ($d - $c + 2*$eps) / $mSubd;
my $ymc = linalg::linspace($c - $eps + $dy/2, $d + $eps - $dy/2, $mSubd);

($xmc, $ymc) = linalg::meshgrid($xmc, $ymc);
$xmc = linalg::flatten($xmc);
$ymc = linalg::flatten($ymc);

my $wSubd = ($b - $a + 2*$eps) / $nSubd / 2;
my $ellSubd = ($d - $c + 2*$eps) / $mSubd / 2;

################################################################################

# Interpolate using polyharmonic spline (phs) radial basis functions (rbf).

my $fe_approx = rbf::interp2($rbfPow, $deg, $x, $y, $f, $xe, $ye
, $xmc, $ymc, $ellSubd, $wSubd);

io::saveArray("$dataDir\\fe_approx.txt", $fe_approx);

$computeTime = time - $computeTime;
print "\$computeTime = $computeTime\n";

# Use python to visualize the results.
system "python plotResults.py $dataDir $checkError";

