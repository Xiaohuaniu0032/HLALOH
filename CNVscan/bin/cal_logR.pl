use strict;
use warnings;
use File::Basename;
use POSIX qw(ceil);
use Data::Dumper;

my ($depth,$ref,$outfile) = @ARGV;

my $name = (split /\./, basename $depth)[0];

my %ref;
open IN, "$ref" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $t = shift @arr;
    my $med = &median_arr(\@arr);
    $ref{$t} = $med;
}
close IN;

#print(Dumper(\%ref));

open O, ">$outfile" or die;
print O "sample\tchr\tstart\tend\tgene\traw_depth\tgc\tgc_corrected_depth\tnorm_depth\tref_median_depth\tlog2Ratio\n";
open IN, "$depth" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $t = "$arr[0]\.$arr[1]\-$arr[2]\.$arr[3]"; # chr.start-end.gene
    my $ref_median = $ref{$t};

    my $log2R;
    if ($arr[-1] == 0 || $ref_median == 0){
        $log2R = "NA"
    }else{
        $log2R = &log2($arr[-1]/$ref_median);
    }
    print O "$name\t$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[-3]\t$arr[-2]\t$arr[-1]\t$ref_median\t$log2R\n";
}
close IN;
close O;





sub log2{
    my $n = shift;
    return log($n)/log(2);
}


sub median_arr{
    my ($val_aref) = @_;
    my $median;
    my @val = sort {$a <=> $b} @{$val_aref};
    if (scalar(@val) % 2 == 0){
        my $v1 = $val[scalar(@val)/2-1];
        my $v2 = $val[scalar(@val)/2];
        $median = sprintf "%.2f", ($v1+$v2)/2;
    }else{
        $median = sprintf "%.2f", $val[ceil(scalar(@val)/2)];
    }

    return($median);
}