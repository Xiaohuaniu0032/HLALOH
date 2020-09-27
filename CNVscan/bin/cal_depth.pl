use strict;
use warnings;
use Getopt::Long;


my ($bam,$name,$bed,$sambamba_bin,$outdir);

GetOptions(
    "bam:s" => \$bam,
    "n:s" => \$name,
    "bed|b:s" => \$bed,
    "sbb:s" => \$sambamba_bin,
    "outdir|od:s" => \$outdir,
    ) or die "unknown args\n";

if (not defined $sambamba_bin){
    $sambamba_bin = "/share/software/sambamba-0.6.6/bin/sambamba";
}

if (not defined $bam || not defined $name || not defined $bed || not defined $outdir){
    die "Usage: perl $0 -bam *.bam -n name -bed *.bed -outdir OUTDIR\n";
}

my $depth_tmp = "$outdir/$name\.depth.tmp";
my $cmd = "$sambamba_bin depth region -L $bed $bam -t 12 >$depth_tmp";
system($cmd) == 0 or die "sambamba depth failed\n";

my $depth_file = "$outdir/$name\.targetcoverage.cnn";
open O, ">$depth_file" or die;
print O "chromosome\tstart\tend\tgene\tdepth\tlog2\n";


open IN, "$depth_tmp" or die;
while (<IN>){
    chomp;
    next if /^\#/;
    my @arr = split /\t/;
    my $depth = int($arr[-2]);
    my $log2_depth;
    if ($depth == 0){
        $log2_depth = -20; # obey cnvkit rule
    }else{
        $log2_depth = sprintf "%.5f", &log2($depth);
    }

    print O "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$depth\t$log2_depth\n";
}
close IN;
close O;


sub log2{
    my $n = shift;
    return (log($n)/log(2));
}

