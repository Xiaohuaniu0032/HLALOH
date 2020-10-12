use strict;
use warnings;

my ($ascat_BAF,$logR,$outfile) = @ARGV;

my %logR; # for each single snp pos
open IN, "$logR" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	# 0-based
	my $start = $arr[2] + 1;
	my $end = $arr[3];
	for my $i ($start..$end){
		$logR{$arr[1]}{$i} = $arr[-1]; # chr=>{pos=>logR}
	}
}
close IN;

open O, ">$outfile";
print O "\tchrs\tpos\tS1\n";

open IN, "$ascat_BAF" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $logR;
	if (exists $logR{$arr[0]}{$arr[1]}){
		$logR = $logR{$arr[0]}{$arr[1]};
	}else{
		# obey BAF SNP info (if one SNP in BAF do not have effective LogR, then set its logR = 0, which means cn is 2)
		$logR = "0"; # set 0 means the copy number is 2
		print "$arr[0].$arr[1].$arr[2] in ascat.BAF file does not have the logR info, please pay atention\n";
	}
	print O "$arr[0]\t$arr[1]\t$arr[2]\t$logR\n";
}
close IN;
close O;
		
