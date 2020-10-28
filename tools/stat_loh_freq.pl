use strict;
use warnings;
use File::Basename;

my ($dir,$of) = @ARGV;

my @files = glob "$dir/*/final_hla.xls";
my $n = scalar(@files);
print "find $n loh result files\n";


open O, ">$of" or die;
print O "sample\tHLA_A1_allele\tHLA_A1_cn\tHLA_A2_allele\tHLA_A2_cn\tHLA_B1_allele\tHLA_B1_cn\tHLA_B2_allele\tHLA_B2_cn\tHLA_C1_allele\tHLA_C1_cn\tHLA_C2_allele\tHLA_C2_cn\tif_LOH\tnumLOH\n";


for my $file (@files){
    my $dir = dirname $file;
    my @loh_sh = glob "$dir/*.HLALOH.sh";
    my $loh_sh = $loh_sh[0];
    my $name = (split /\./, basename $loh_sh)[0];

    my $line = (split /\s/, `wc -l $file`)[0];
    if ($line != 7){
        print "$file is empty, will be skipped\n";
        next;
    }

    open IN, "$file" or die;
    <IN>;
    my $hla_a1 = <IN>;my $hla_a2 = <IN>;
    my $hla_b1 = <IN>;my $hla_b2 = <IN>;
    my $hla_c1 = <IN>;my $hla_c2 = <IN>;
    close IN;

    my @val_a1 = split /\t/, $hla_a1;
    my @val_a2 = split /\t/, $hla_a2;

    my @val_b1 = split /\t/, $hla_b1;
    my @val_b2 = split /\t/, $hla_b2;

    my @val_c1 = split /\t/, $hla_c1;
    my @val_c2 = split /\t/, $hla_c2;

    my $hla_a = "$val_a1[0]\;$val_a2[0]";
    my $hla_b = "$val_b1[0]\;$val_b2[0]";
    my $hla_c = "$val_c1[0]\;$val_c2[0]";

    my ($hla_a_loh,$hla_b_loh,$hla_c_loh) = (0,0,0);

    # determine loh only by cn <= 0.5
    my $cn_cutoff = 0.5;

    if ($val_a1[1] eq "NA" || $val_a2[1] eq "NA"){
        $hla_a_loh = 0;
    }else{
        if ($val_a1[1] <= $cn_cutoff || $val_a2[1] <= $cn_cutoff){
            $hla_a_loh = 1;
        }else{
            $hla_a_loh = 0;
        }
    }

    if ($val_b1[1] eq "NA" || $val_b2[1] eq "NA"){
        $hla_b_loh = 0;
    }else{
        if ($val_b1[1] <= $cn_cutoff || $val_b2[1] <= $cn_cutoff){
            $hla_b_loh = 1;
        }else{
            $hla_b_loh = 0;
        }
    }

    if ($val_c1[1] eq "NA" || $val_c2[1] eq "NA"){
        $hla_c_loh = 0;
    }else{
        if ($val_c1[1] <= $cn_cutoff || $val_c2[1] <= $cn_cutoff){
            $hla_c_loh = 1;
        }else{
            $hla_c_loh = 0;
        }
    }

    my $if_loh;
    if ($hla_a_loh + $hla_b_loh + $hla_c_loh == 0){
        $if_loh = "NO";
    }else{
        $if_loh = "YES";
    }

    my $nLOH = $hla_a_loh + $hla_b_loh + $hla_c_loh;

    print O "$name\t$val_a1[0]\t$val_a1[1]\t$val_a2[0]\t$val_a2[1]\t$val_b1[0]\t$val_b1[1]\t$val_b2[0]\t$val_b2[1]\t$val_c1[0]\t$val_c1[1]\t$val_c2[0]\t$val_c2[1]\t$if_loh\t$nLOH\n";

}
close O;