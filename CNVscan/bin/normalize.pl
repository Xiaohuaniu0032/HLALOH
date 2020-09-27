use strict;
use warnings;
use Getopt::Long;
use List::Util qw/uniq/;
use File::Basename;

my ($depthF,$outdir);

GetOptions(
    "depth|d:s" => \$depthF,
    "outdir|od:s" => \$outdir,
    ) or die "unknown args\n";

my $name = (split /\./, basename $depthF)[0];
print "Normalizing $name\n";
my $ratioF = "$outdir/$name\.norm.xls";
open O, ">$ratioF" or die;

open IN, "$depthF" or die;
my $h = <IN>;
chomp $h;
close IN;
print O "$h\tnorm_depth\n";

my $flag = 0;
my $aref_norm_depth = &lib_size_norm($depthF);
open IN, "$depthF" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $t = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
    my $norm = $aref_norm_depth->[$flag];
    print O "$_\t$norm\n";
    $flag += 1;
}
close IN;
close O;
undef $flag;


sub lib_size_norm{
    # by median cov of all exon, the cov of chrX/Y has been corrected to diploid
    my ($file) = @_;
    my @norm_depth;
    my @cov;
    open IN, "$file" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        push @cov, $arr[-1];
    }
    close IN;

    my $med_cov = &get_median(\@cov);
    print "median cov is $med_cov\n";

    open IN, "$file" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        my $norm_depth = sprintf "%.3f", $arr[-1]/$med_cov;
        push @norm_depth, $norm_depth;
    }
    close IN;

    return(\@norm_depth);
}

sub get_median{
    my ($aref) = @_;
    my @arr = sort {$a <=> $b} @{$aref}; # 由大到小排序
    my $num = scalar(@arr);
    my $median;
    if ($num % 2 == 0){
        # 偶数
        my $idx_1 = $num/2 - 1;
        my $idx_2 = $idx_1 + 1;
        my $val_1 = $arr[$idx_1];
        my $val_2 = $arr[$idx_2];
        $median = sprintf "%.3f", ($val_1+$val_2)/2;
    }else{
        my $idx = ($num + 1) / 2 - 1;
        $median = $arr[$idx];
    }

    return($median);
}
