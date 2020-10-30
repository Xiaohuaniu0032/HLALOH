use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($het_pos_file,$allele_1_pileup,$allele_2_pileup,$outfile);
GetOptions(
    "het:s" => \$het_pos_file,
    "a1:s" => \$allele_1_pileup,
    "a2:s" => \$allele_2_pileup,
    "of:s" => \$outfile,
    ) or die "unknown args\n";


my %hla_pileup;
open IN, "$allele_1_pileup" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $hla_pileup{$arr[0]}{$arr[1]} = "$arr[2]\t$arr[3]"; # allele=>{pos=>{"base\tcount"}}
}
close IN;

open IN, "$allele_2_pileup" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $hla_pileup{$arr[0]}{$arr[1]} = "$arr[2]\t$arr[3]";
}
close IN;

open O, ">$outfile" or die;
print O "Rank\tpos1\tpos2\tbase1\tbase2\tnum1\tnum2\tdepth\tbaf1\tbaf2\n";

my $flag = 0;
open IN, "$het_pos_file" or die;
my $h = <IN>;
my @h = split /\t/, $h;
my $a1 = $h[0];
my $a2 = $h[1];

while (<IN>){
    chomp;
    my @arr = split /\t/;
    if ($arr[2] == 1 and $arr[3] == 1){
        # het pos
        if (exists $hla_pileup{$a1}{$arr[0]} and exists $hla_pileup{$a2}{$arr[1]}){
            my @v1 = split /\t/, $hla_pileup{$a1}{$arr[0]};
            my @v2 = split /\t/, $hla_pileup{$a2}{$arr[1]};

            my $n1 = $v1[1];
            my $n2 = $v2[1];

            if ($n1 < 5 || $n2 < 5){
                # skip low depth positions
                next;
            }

            my $base1 = $v1[0];
            my $base2 = $v2[0];

            my ($baf1,$baf2);
            my $all_n = $n1 + $n2;
            if ($all_n == 0){
                next;
                # skip this pos
            }else{
                $baf1 = sprintf "%.2f", $n1/($n1+$n2);
                $baf2 = sprintf "%.2f", $n2/($n1+$n2);
            }

            $flag += 1;

            print O "$flag\t$arr[0]\t$arr[1]\t$base1\t$base2\t$n1\t$n2\t$all_n\t$baf1\t$baf2\n";
        }
    }
}
close IN;
