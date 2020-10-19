use strict;
use warnings;
use Getopt::Long;

my ($predict_file,$hla_allele,$cn_cutoff,$p_value_cutoff,$outfile);
GetOptions(
    "res:s" => \$predict_file,        # NEED
    "hla:s" => \$hla_allele,          # NEED (hla.result.new)
    "cn:f" => \$cn_cutoff,            # 0.5
    "p:f" => \$p_value_cutoff,        # 0.01
    "o:s" => \$outfile,               # NEED
    ) or die;

# default value
if (not defined $cn_cutoff){
    $cn_cutoff = 0.5;
}
if (not defined $p_value_cutoff){
    $p_value_cutoff = 0.01;
}

open O, ">$outfile" or die;
print O "hla_allele\tcopyNumber\tLOH\n";

my @alleles;
open IN, "$hla_allele" or die;
while (<IN>){
    chomp;
    push @alleles, $_;
}
close IN;

my %hla_predict;
open IN, "$predict_file" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $hla_predict{$arr[1]} = "$arr[22]"; # HLA_type1copyNum_withBAFBin
    $hla_predict{$arr[2]} = "$arr[28]"; # HLA_type2copyNum_withBAFBin
    my $abc = (split /\_/, $arr[1])[1];
    $hla_predict{$abc} = $arr[-6]; # PVal_unique
}
close IN;

for my $allele (@alleles){
    #print "$allele\n";
    my $abc = (split /\_/, $allele)[1];
    if (exists $hla_predict{$allele}){
        my $cn = sprintf "%.2f", $hla_predict{$allele};
        my $pval = $hla_predict{$abc};
        if ($cn <= $cn_cutoff and $pval <= $p_value_cutoff){
            # a loh
            print O "$allele\t$cn\tYES\n";
        }else{
            print O "$allele\t$cn\tNO\n";
        }
    }else{
        # no exists in result file
        print O "$allele\tNA\tNO\n";
    }
}

close O;