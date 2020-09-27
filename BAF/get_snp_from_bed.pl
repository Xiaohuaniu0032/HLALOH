use strict;
use warnings;

my ($bed,$hg19_EAS,$outfile) = @ARGV;

# /data1/software/annovar/humandb/hg19_EAS.sites.2014_10.txt
my %snp;
open IN, "$hg19_EAS" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $ref_len = length($arr[2]);
    my $alt_len = length($arr[3]);
    next if ($ref_len != 1 or $alt_len != 1); # skip useless snp pos
    $snp{$arr[0]}{$arr[1]} = $_;
}
close IN;

open O, ">$outfile" or die;

open BED, "$bed" or die;
while (<BED>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }

    my $start = $arr[1] + 1;
    my $end = $arr[2];

    for my $i ($start..$end){
        if (exists $snp{$chr}{$i}){
            print O "$snp{$chr}{$i}\n";
        }
    }
}
close BED;
close O;
