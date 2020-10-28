use strict;
use warnings;
use File::Basename;

my ($logR_BAF,$eas_snp_bed,$outfile) = @ARGV;

my $name = (split /\./, basename $logR_BAF)[0];

open O, ">$outfile" or die;
print O "\tchrs\tpos\t$name\n";

my %snp_pos;
open IN, "$eas_snp_bed" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/, $_;
    $snp_pos{$arr[0]}{$arr[2]} = $arr[-1]; # chr=>{pos=>rs}
}
close IN;

my $flag = 0;
open IN, "$logR_BAF" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    if (exists $snp_pos{$arr[1]}{$arr[2]}){
        $flag += 1;
        my $snp = "SNP".$flag;
        print O "$snp\t$arr[1]\t$arr[2]\t$arr[3]\n";
    }
}
close IN;
close O;

