use strict;
use warnings;
use File::Basename;

my ($resdir) = @ARGV;

my @files = glob "$resdir/*/*.HLALOH.sh";

for my $f (@files){
    my $d = dirname($f);
    my $name = (split /\./, basename($f))[0];
    my $a = "$d/hla_a_BAF.txt";
    my $b = "$d/hla_b_BAF.txt";
    my $c = "$d/hla_c_BAF.txt";
    my $cnv = "$d/$name\.CopyNumber.xls";

    if (-e $a and -e $b and -e $c and -e $cnv){
        print "determine_HLALOH_by_BAF for $name\n";
        `perl /home/fulongfei/workdir/git_repo/HLALOH/tools/determine_HLALOH_by_BAF.pl $d $name`;
    }
}
