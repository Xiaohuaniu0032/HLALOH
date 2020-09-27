use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($bed,$depth_cnn,$bedtools_bin,$fasta,$outdir);
GetOptions(
    "bed:s" => \$bed,
    "depth:s" => \$depth_cnn,
    "bedtools_bin:s" => \$bedtools_bin,
    "fa:s" => \$fasta,
    "outdir:s" => \$outdir,
    ) or die "unknown args\n";


if (not defined $bedtools_bin){
    $bedtools_bin = "/share/public/software/bedtools2/bin/bedtools";
}

if (not defined $fasta){
    $fasta = "/share/public/software/Onc_Soft/database/hg19/hg19.rm_CRLF2_P2RY8.fasta";
}

# cal gc
my $gc_file = "$outdir/gc.xls";
my $cmd = "$bedtools_bin nuc -fi $fasta -bed $bed >$gc_file";
system($cmd) == 0 or die "bedtools gc fail\n";

# add gc to depth file
my %target2gc;
open IN, "$gc_file" or die;
my $h = <IN>;
close IN;

my $gc_col;
my @h = split /\t/, $h;
my $flag = 0;
for (@h){
    $flag += 1;
    if (/pct_gc/){
        $gc_col = $flag;
    }
}
undef $flag;

print "gc col is $gc_col\n";

open IN, "$gc_file" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $target = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
    $target2gc{$target} = $arr[($gc_col-1)];
}
close IN;

my $name = basename $depth_cnn;
my $depth_gc = "$outdir/$name\.with.gc.xls";
open O, ">$depth_gc" or die;
print O "chromosome\tstart\tend\tgene\tdepth\tlog2\tgc\n";

open IN, "$depth_cnn" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $target = "$arr[0]\t$arr[1]\t$arr[2]";
	my $gc = int($target2gc{$target} * 100); # 0,1,2,...,100%
    print O "$_\t$gc\n";
}
close IN;
close O;



