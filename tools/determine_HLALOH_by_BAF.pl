use strict;
use warnings;
#use Getopt::Long;
use POSIX qw(ceil);

my ($hla_a_baf,$hla_b_baf,$hla_c_baf,$resdir,$name) = @ARGV;

#my $hla_a_baf = "$resdir/hla_a_BAF.txt"; # maybe empty
#my $hla_b_baf = "$resdir/hla_b_BAF.txt";
#my $hla_c_baf = "$resdir/hla_c_BAF.txt";

my $of = "$resdir/$name\.hlaloh.result.xls";


# default cutoff
my $loh_cutoff_lower = 0.4;
my $loh_cutoff_upper = 0.6;
my $loh_pct_cutoff = 0.6;



open O, ">$of" or die;
print O "sample\tA_het_pos_num\tA_loh_pos_num\tA_loh_pos_pct\tmedian_A_BAF\tHLA_A_CopyNumber\tif_A_loh\tB_het_pos_num\tB_loh_pos_num\tB_loh_pos_pct\tmedian_B_BAF\tHLA_B_CopyNumber\tif_B_loh\tC_het_pos_num\tC_loh_pos_num\tC_loh_pos_pct\tmedian_C_BAF\tHLA_C_CopyNumber\tif_C_loh\n";


my $cnFile = "$resdir/$name\.CopyNumber.xls";
if (!-e $cnFile){
    die "can not find $cnFile\n";
}


my $hla_a_info = &stat_loh_pos($hla_a_baf,$cnFile,"A");
my $hla_b_info = &stat_loh_pos($hla_b_baf,$cnFile,"B");
my $hla_c_info = &stat_loh_pos($hla_c_baf,$cnFile,"C");

my @a = split /\t/, $hla_a_info;
my @b = split /\t/, $hla_b_info;
my @c = split /\t/, $hla_c_info;

my $a_loh = &if_loh($hla_a_info);
my $b_loh = &if_loh($hla_b_info);
my $c_loh = &if_loh($hla_c_info);


print O "$name\t$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a_loh\t$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b_loh\t$c[0]\t$c[1]\t$c[2]\t$c[3]\t$c[4]\t$c_loh\n";

close O;

sub if_loh{
    my ($loh_info) = @_;
    my $loh;
    my @val = split /\t/, $loh_info; # my $val = "$het_pos_n\t$loh_pos_n\t$loh_pct\t$median_baf\t$hla_cn";
    #print($val[2]);
    #print "$loh_info\n";
    if ($val[2] > $loh_pct_cutoff){ # 0.4~0.6之外的点的百分比，默认60%
        $loh = "YES";
    }else{
        $loh = "NO";
    }

    return($loh);
}


sub stat_loh_pos{
    # 统计a/b/c每个allele的信息（杂合位点个数，loh位点个数，BAF中位值，是否为LOH）
    my ($baf_file,$cn_file,$gene) = @_;

    my ($het_pos_n,$loh_pos_n,$loh_pct,$median_baf,$hla_cn) = (0,0,0,0,0);
    
    my @baf;
    open IN, "$baf_file" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        $het_pos_n += 1;
        if ($arr[-1] > $loh_cutoff_upper || $arr[-1] < $loh_cutoff_lower){
            # BAF > 0.6 OR BAF < 0.4
            $loh_pos_n += 1;
        }
        push @baf, $arr[-1];
    }
    close IN;


    if ($het_pos_n != 0){
        $loh_pct = sprintf "%.2f", $loh_pos_n/$het_pos_n;
    }else{
        $loh_pct = 0;
    }

    if ($het_pos_n == 0){
        # empty file
        $median_baf = "NA";
    }else{
        $median_baf = &cal_median(\@baf);
    }

    
    my %gene_cn;
    #my @genes;
    open IN, "$cnFile" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        #$gene_cn{$arr[1]} = "$arr[2]\t$arr[3]\t$arr[-2]"; # gene=>chr/start/cn
        $gene_cn{$arr[1]} = $arr[-2];
        #push @genes,$arr[1];
    }
    close IN;

    if ($gene eq "A"){
        if (exists $gene_cn{"HLA-A"}){
            $hla_cn = $gene_cn{"HLA-A"};
        }else{
            $hla_cn = "NA";
        }
    }

    if ($gene eq "B"){
        if (exists $gene_cn{"HLA-B"}){
            $hla_cn = $gene_cn{"HLA-B"};
        }else{
            $hla_cn = "NA";
        }
    }

    if ($gene eq "C"){
        if (exists $gene_cn{"HLA-C"}){
            $hla_cn = $gene_cn{"HLA-C"};
        }else{
            $hla_cn = "NA";
        }
    }

    #print "$hla_cn\n";

    my $val = "$het_pos_n\t$loh_pos_n\t$loh_pct\t$median_baf\t$hla_cn";
    print("$val\n");
    return($val);
}



sub cal_median{
    my ($val_aref) = @_;
    my $median;
    my @val = sort {$a <=> $b} @{$val_aref};
    if (scalar(@val) % 2 == 0){
        my $v1 = $val[scalar(@val)/2-1];
        my $v2 = $val[scalar(@val)/2];
        $median = sprintf "%.2f", ($v1+$v2)/2;
    }else{
        # if only one value
        if (@val == 1){
            $median = $val[0];
        }else{
            $median = sprintf "%.2f", $val[ceil(scalar(@val)/2)];
        }
    }

    return($median);
}



