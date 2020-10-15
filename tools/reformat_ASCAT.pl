use strict;
use warnings;

my ($in,$out) = @ARGV;

my ($purity,$ploidy);
open IN, "$in" or die;
my $v1 = <IN>;
chomp $v1;
$purity = sprintf "%.2f", (split /\s+/, $v1)[1];
my $v2 = <IN>;
chomp $v2;
$ploidy = sprintf "%.2f", (split /\s+/, $v2)[1];

open O, ">$out" or die;
print O "Ploidy\ttumorPurity\ttumorPloidy\n";
print O "example_tumor_sorted\t$ploidy\t$purity\t$ploidy\n";
close O;
