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

open O, ">$outfile" or die; # hg19_EAS.sites.2014_10.txt is 1-based, we need 0-based bed file, which will be used to calculate specific SNP's BAF

open BED, "$bed" or die;
while (<BED>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }

    my $start = $arr[1] + 1; # BED file
    my $end = $arr[2];

    for my $i ($start..$end){
        if (exists $snp{$chr}{$i}){
            # get all snp that included by the bed file
            my $val = $snp{$chr}{$i};
            my @val = split /\t/, $val;
            my $pos = $val[1]; # 1-based
            my $new_pos_start = $pos - 1;
            my $new_pos_end = $pos;
            print O "$val[0]\t$new_pos_start\t$new_pos_end\t$val[2]\t$val[3]\t$val[4]\t$val[5]\n";
        }
    }
}
close BED;
close O;
