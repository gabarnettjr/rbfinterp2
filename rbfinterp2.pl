# windows

# This is a perl script that can be called to interpolate (x, y, f) to 
# (xe, ye, fe_approx).  In order for this script to work, you need to have some
# text files saved in a common location, so they can be loaded in by the script.
# These four files should be saved in a coordinate directory:
# x.txt, y.txt, xe.txt, ye.txt (nodes and evaluation points)
# The following should be in a subfolder ($dataDir) of the coordinate directory:
# f.txt (function values at the nodes)
# In addition, if you wish to compare to some known exact values at the
# evaluation points, then you should also include this file in $dataDir:
# fe.txt (function values at the evaluation points)

# If you have all of this, then the script will produce the file fe_approx.txt,
# which will contain estimated values of the function at the evaluation points.

# Greg Barnett
# December 2022

################################################################################

use strict;
use warnings;

use lib ".";
use io;
use rbf2;

################################################################################

my ($dataDir, $checkError, $rbfPow, $deg, $nSubd, $mSubd);

$dataDir = "cdCoords\\smoothData";
$checkError = "y";
$rbfPow = 3;
$deg = 1;
$nSubd = "";
$mSubd = "";

if (scalar @ARGV) {
    if ((lc $ARGV[0]) =~ /^\-\-help$|^\-h$|^help$|^h$/) {
        print helpString(); exit;
    }
    $dataDir = shift;
    $checkError = "n";
}
if (scalar @ARGV) {
    $checkError = shift;
	$checkError = (lc $checkError);
	if ($checkError =~ /y/) {
		$checkError = "y";
	} elsif ($checkError =~ /n/) {
		$checkError = "n";
	} else {
		print STDERR "Invalid second input (\$checkError).  Should be \"y\" or \"n\".\n"; die;
	}
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
	print STDERR "Too many inputs.  Max number of inputs is 6.\n"; die;
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

################################################################################

sub helpString {
	my $s = "\nPlease use up to six inputs.\n";
    $s .= "(1) The first input is the folder containing input function data.\n";
    $s .= "(2) The second input is whether or not to calculate the error,\n";
    $s .= "    where \"y\" means calculate the error, and \"n\" means don\'t.\n";
    $s .= "(3) The third input is the rbf exponent (default 3).\n";
    $s .= "(4) The fourth input is the polynomial degree (default 1).\n";
    $s .= "(5) The fifth input is number of subdomains across (default auto calculate).\n";
    $s .= "(6) The sixth input is number of subdomains going down (default auto calculate).\n\n";
	return $s;
}

################################################################################

