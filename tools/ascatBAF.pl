use strict;
use warnings;

# translate BAF into ascat's format


my ($vaf,$outfile) = @ARGV;


open O, ">$outfile" or die;
print O "\tchrs\tpos\tS1\n";

my $flag = 0;
open IN, "$vaf" or die;
<IN>;
while (<IN>){
	chomp;
	$flag += 1;
	my @arr = split /\t/;
	my $chr = $arr[0];
	if ($chr =~ /^chr/){
		$chr =~ s/^chr//;
	}
	my $snp = "SNP".$flag;
	print O "$snp\t$arr[0]\t$arr[1]\t$arr[-1]\n";
}
close IN;
close O;

