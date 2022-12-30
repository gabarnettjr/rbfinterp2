# windows

use strict;
use warnings;

use Math::Trig;

use lib ".";
use io;
use linalg;

my $dataDir;
if ((scalar @ARGV) == 1) {
    $dataDir = shift;
} else {
    $dataDir = "cdCoords\\smoothData";
}

################################################################################

my $x = io::loadArray("$dataDir\\..\\x.txt");
my $y = io::loadArray("$dataDir\\..\\y.txt");

my $xe = io::loadArray("$dataDir\\..\\xe.txt");
my $ye = io::loadArray("$dataDir\\..\\ye.txt");

################################################################################

# sub eff {
    # my ($x, $y) = @_;
    # my $z = linalg::zeros(scalar @{$x});
    # for (my $i = 0; $i < (scalar @{$x}); $i++) {
        # @{$z}[$i] = cos(2*pi * @{$x}[$i] / 25000) * sin(2*pi * @{$y}[$i] / 25000);
    # }
    # return $z;
# }

sub eff {
    my ($x, $y) = @_;
    my $z = linalg::zeros(scalar @{$x});
    for (my $i = 0; $i < (scalar @{$x}); $i++) {
        @{$z}[$i] = cos(2*pi * @{$x}[$i] / 25000) * cos(2*pi * @{$y}[$i] / 25000)
                  + sin(2*pi * @{$x}[$i] / 25000) * sin(2*pi * @{$y}[$i] / 25000);
    }
    return $z;
}

# sub eff {
    # my ($x, $y) = @_;
    # my $z = linalg::zeros(scalar @{$x});
    # for (my $i = 0; $i < (scalar @{$x}); $i++) {
        # @{$z}[$i] = cos(2*pi * @{$x}[$i] / 50000) * cos(2*pi * @{$y}[$i] / 50000)
                  # + sin(2*pi * @{$x}[$i] / 50000) * sin(2*pi * @{$y}[$i] / 50000);
    # }
    # return $z;
# }

# def eff(x,y):
    # return np.cos(2*np.pi*x/50000) * np.cos(2*np.pi*y/50000) \
    # + np.sin(2*np.pi*x/50000) * np.sin(2*np.pi*y/50000)

################################################################################

io::saveArray("$dataDir\\f.txt", eff($x, $y));

io::saveArray("$dataDir\\fe.txt", eff($xe, $ye));

