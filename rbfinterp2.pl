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

# Process the input.

my ($dataDir, $checkError, $rbfPow, $deg, $nSubd, $mSubd);

$dataDir = "randomCoords\\smoothData";
$checkError = "y";
$rbfPow = 3;
$deg = 1;
$nSubd = -1;
$mSubd = -1;

if (scalar @ARGV) {
    if ((lc $ARGV[0]) =~ /^\-\-help$|^\-h$|^help$|^h$/) {
        print helpString(); exit;
    }
    $dataDir = shift;
	if (! -d $dataDir) {
		print STDERR "First input must be a data directory.\n"; die;
	} elsif (! -e "$dataDir\\f.txt") {
		print STDERR "Data directory must contain function values.\n"; die;
	}
    $checkError = "n";
}
if (scalar @ARGV) {
    $checkError = shift;
	$checkError = (lc $checkError);
	if (($checkError =~ /y/) && ($checkError !~ /n/)) {
		$checkError = "y";
	} elsif (($checkError =~ /n/) && ($checkError !~ /y/)) {
		$checkError = "n";
	} else {
		print STDERR "Invalid second input (\$checkError).  Should be \"y\" or \"n\".\n"; die;
	}
}
if (scalar @ARGV) { $rbfPow = shift; }
if (scalar @ARGV) { $deg = shift; }
if (scalar @ARGV) { $nSubd = shift; }
if (scalar @ARGV) { $mSubd = shift; }
if (scalar @ARGV) {
	print STDERR "Too many inputs.  Max number of inputs is 6.\n"; die;
}

################################################################################

# Interpolate using polyharmonic spline (PHS) radial basis functions (RBFs).

if (-e "$dataDir\\fe_approx.txt") {
	unlink "$dataDir\\fe_approx.txt";
}

# # USE PERL:
# my $x  = io::loadArray("$dataDir\\..\\x.txt");
# my $y  = io::loadArray("$dataDir\\..\\y.txt");
# my $f  = io::loadArray("$dataDir\\f.txt");
# my $xe = io::loadArray("$dataDir\\..\\xe.txt");
# my $ye = io::loadArray("$dataDir\\..\\ye.txt");
# my $computeTime = time;
# my $fe_approx = rbf2::interp($x, $y, $f, $xe, $ye, $rbfPow, $deg, $nSubd, $mSubd);
# $computeTime = time - $computeTime;
# print "\$computeTime = $computeTime\n";
# io::saveArray("$dataDir\\fe_approx.txt", $fe_approx);

# # USE JULIA:
# system "julia .\\julia\\rbfinterp2.jl $dataDir $rbfPow $deg $nSubd $mSubd";

# USE PYTHON:
system "python .\\python\\rbfinterp2.py $dataDir $rbfPow $deg $nSubd $mSubd";

if (! -e "$dataDir\\fe_approx.txt") {
	print STDERR "Please investigate error during fe_approx.txt creation.\n"; die;
}

################################################################################

# Use python to visualize the results.
system "python plotResults.py $dataDir $checkError";

################################################################################

sub helpString {
	my $s = "\n";
	$s .= "\"rbfinterp2.pl\": A script for interpolating (x,y,f) to (xe,ye,fe_approx).\n\n";
	$s .= "To run the script, you need a coordinates directory that contains these:\n";
	$s .= "x.txt, y.txt, xe.txt, ye.txt\n\n";
	$s .= "Inside the coordinates directory should be a function value subdirectory, containing these:\n";
	$s .= "f.txt  (required)\n";
	$s .= "fe.txt (optional, but needed for error calculation)\n\n";
	$s .= "This script accepts up to 6 command-line inputs.\n";
    $s .= "(1) The first input is the path to the folder that contains your function values (default: .\\randomCoords\\smoothData).\n";
    $s .= "(2) The second input is whether or not to calculate the error, y or n            (default: n).\n";
    $s .= "(3) The third input is the rbf exponent, an odd integer                          (default: 3).\n";
    $s .= "(4) The fourth input is the polynomial degree, an integer from 0 up to 4         (default: 1).\n";
    $s .= "(5) The fifth input is number of subdomains across, a positive integer           (default: auto calculate).\n";
    $s .= "(6) The sixth input is number of subdomains going down, a positive integer       (default: auto calculate).\n\n";
	return $s;
}

################################################################################

