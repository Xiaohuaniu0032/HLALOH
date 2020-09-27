use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw/$Bin/;
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);
use Data::Dumper;

my ($tbam,$nbam,$bed,$tbaf,$nbaf,$tlogR,$nlogR,)