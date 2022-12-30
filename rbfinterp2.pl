# windows

use strict;
use warnings;

use lib ".";
use io;
use rbf2;

################################################################################

my ($dataDir, $checkError, $rbfPow, $deg, $nSubd, $mSubd);

$dataDir = "cdCoords\\smoothData";
$checkError = 1;
$rbfPow = 3;
$deg = 1;
$nSubd = "";
$mSubd = "";

if (scalar @ARGV) {
    if ($ARGV[0] =~ /^\-\-help$|^\-h$/) {
        goto HELP;
    }
    $dataDir = shift;
    $checkError = 0;
}
if (scalar @ARGV) {
    $checkError = shift;
}
if (scalar @ARGV) {
    $rbfPow = shift;
}
if (scalar @ARGV) {
    $deg = shift;
}
if (scalar @ARGV) {
    $nSubd = shift;
}
if (scalar @ARGV) {
    $mSubd = shift;
}
if (scalar @ARGV) {
    HELP:
    my $s = "\nPlease use up to six inputs.\n";
    $s .= "(1) The first input is the folder containing input function data.\n";
    $s .= "(2) The second input is whether or not to calculate the error,\n";
    $s .= "    where 1 means calculate the error, and 0 means don\'t.\n";
    $s .= "(3) The third input is the rbf exponent (default 3).\n";
    $s .= "(4) The fourth input is the polynomial degree (default 1).\n";
    $s .= "(5) The fifth input is number of subdomains across (default auto calculate).\n";
    $s .= "(6) The sixth input is number of subdomains going down (default auto calculate).\n\n";
    print STDERR $s; die;
}

################################################################################

# Start the timer for loading/computing/saving in perl.
my $computeTime = time;

# Load the data at the nodes (where you know the function).
my $x = io::loadArray("$dataDir\\..\\x.txt");
my $y = io::loadArray("$dataDir\\..\\y.txt");
my $f = io::loadArray("$dataDir\\f.txt");

# Load the evaluation points (where you WANT to know the function).
my $xe = io::loadArray("$dataDir\\..\\xe.txt");
my $ye = io::loadArray("$dataDir\\..\\ye.txt");

################################################################################

# Interpolate using polyharmonic spline (PHS) radial basis functions (RBFs).
my $fe_approx = rbf2::interp($x, $y, $f, $xe, $ye
, $rbfPow, $deg, $nSubd, $mSubd);

# Save the interpolated values.
io::saveArray("$dataDir\\fe_approx.txt", $fe_approx);

# Print the compute time.
$computeTime = time - $computeTime;
print "\$computeTime = $computeTime\n";

# Use python to visualize the results.
system "python plotResults.py $dataDir $checkError";

