use strict;
use warnings;
#use Getopt::Long;
use POSIX qw(ceil);

my ($baf_dir,$normal_name,$tumor_name,$outdir) = @ARGV;

my $of = "$outdir/$tumor_name\.hlaloh.result.xls";

# default cutoff
my $loh_cutoff_lower = 0.4;
my $loh_cutoff_upper = 0.6;

my $loh_pct_cutoff = 0.65;
my $median_baf_cutoff = 0.6;

# how to determine the LOH by BAF distribution?
# median BAF & sd(BAF)

open O, ">$of" or die;
print O "Sample\tGene\tHetPosNum\tLikelyLoHPosNum\tLoHPosPercent(%)\tMedianBAF\tCopyNumber\tIfLoH\n";


# tumor copy number
my $cnFile = "$baf_dir/$tumor_name\.CopyNumber.xls";

if (!-e $cnFile){
    die "can not find $cnFile\n";
}


my $hla_a_baf_N = "$baf_dir/$normal_name\.hla_a_BAF.txt";
my $hla_a_baf_T = "$baf_dir/$tumor_name\.hla_a_BAF.txt";
my $hla_a_info = &stat_loh_info($hla_a_baf_N,$hla_a_baf_T,$cnFile,'HLA-A');


my $hla_b_baf_N = "$baf_dir/$normal_name\.hla_b_BAF.txt";
my $hla_b_baf_T = "$baf_dir/$tumor_name\.hla_b_BAF.txt";
my $hla_b_info = &stat_loh_info($hla_b_baf_N,$hla_b_baf_T,$cnFile,'HLA-B');

my $hla_c_baf_N = "$baf_dir/$normal_name\.hla_c_BAF.txt";
my $hla_c_baf_T = "$baf_dir/$tumor_name\.hla_c_BAF.txt";
my $hla_c_info = &stat_loh_info($hla_c_baf_N,$hla_c_baf_T,$cnFile,'HLA-C');

my $a_loh = &determine_LOH($hla_a_info);
my $b_loh = &determine_LOH($hla_b_info);
my $c_loh = &determine_LOH($hla_c_info);

print O "$tumor_name\t$hla_a_info\t$a_loh\n";
print O "$tumor_name\t$hla_b_info\t$b_loh\n";
print O "$tumor_name\t$hla_c_info\t$c_loh\n";

close O;




sub determine_LOH{
    my ($loh_info) = @_;
    my @val = split /\t/, $loh_info; # my $val = "$gene\t$het_pos_n\t$loh_pos_n\t$loh_pct\t$median_baf\t$hla_cn";

    # how to determine LOH by BAF distribution?
    # 1. loh pos pct >= 65% and
    # 2. median BAF >= 0.6

    # if het_pos_n <= 10, then can not determine LOH status

    my $loh_status;

    if ($val[1] <= 10){
        $loh_status = 'NA';
        print "can not determine LOH status for $val[0]: het pos is $val[1] (<=10)\n";

        return($loh_status);
    }

    if ($val[3] > $loh_pct_cutoff and $val[-2] >= $median_baf_cutoff){
        # 0.65 and 0.6
        $loh_status = 'YES';
    }else{
        $loh_status = 'NO';
    }

    return($loh_status);
}


sub stat_loh_info{
    # 统计A/B/C每个allele信息:杂合位点个数/可能的LOH位点个数/LOH位点百分比/BAF中位值/拷贝数
    my ($normal_baf,$tumor_baf,$tumor_cn_file,$gene) = @_;

    my %eff_het_pos_normal;
    open IN, "$normal_baf" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        next if ($arr[-3] <= 50); # skip low depth pos
        if ($arr[-1] >= $loh_cutoff_lower and $arr[-1] <= $loh_cutoff_upper){
            $eff_het_pos_normal{$arr[1]} = 1;
        }
    }
    close IN;

    # how many eff het pos in normal?
    my $eff_het_n = scalar(keys %eff_het_pos_normal);

    my ($het_pos_n,$loh_pos_n,$loh_pct,$median_baf,$hla_cn) = (0,0,0,0,0);

    my @baf;
    open IN, "$tumor_baf" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        if (exists $eff_het_pos_normal{$arr[1]}){
            print "$_\n";
            # this pos is ok in normal, and can be used to stat tumor baf info
            # stat baf info
            $het_pos_n += 1;
            
            if ($arr[-1] > 0.5){
                push @baf, $arr[-1];
            }else{
                push @baf, $arr[-2];
            }

            if ($arr[-1] <= $loh_cutoff_lower or $arr[-1] >= $loh_cutoff_upper){
                $loh_pos_n += 1;
            }
        }
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
    open IN, "$cnFile" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        $gene_cn{$arr[1]} = $arr[-2]; # gene=>cn
    }
    close IN;

    # gene is HLA-A or HLA-B or HLA-C
    if (exists $gene_cn{$gene}){
        $hla_cn = $gene_cn{$gene};
    }else{
        $hla_cn = "NA";
    }

    my $val = "$gene\t$het_pos_n\t$loh_pos_n\t$loh_pct\t$median_baf\t$hla_cn";
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



