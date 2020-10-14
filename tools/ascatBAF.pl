use strict;
use warnings;

# for translate file format
# input: pileup2vaf.v2.py's outfile
# output: ASCAT's BAF file format
#

my ($vaf,$snp_pos,$outfile) = @ARGV;

my %vaf;
open IN, "$vaf" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$vaf{$arr[0]}{$arr[1]} = $arr[-1];
}
close IN;

open O, ">$outfile" or die;
print O "\tchrs\tpos\tS1\n";

my $flag = 0;
open IN, "$snp_pos" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	if (exists $vaf{$arr[0]}{$arr[1]}){
		$flag += 1;
		my $snp = "SNP".$flag;
		my $vaf = $vaf{$arr[0]}{$arr[1]};
		print O "$snp\t$arr[0]\t$arr[1]\t$vaf\n";
	}
}
close IN;
close O;

