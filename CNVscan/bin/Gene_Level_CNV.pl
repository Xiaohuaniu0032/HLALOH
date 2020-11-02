use strict;
use warnings;
use Getopt::Long;

my ($logR,$gain_cutoff,$loss_cutoff,$target_num_cutoff,$cnv_pct_cutoff,$outfile);

GetOptions(
    "in|i:s" => \$logR,                    # Need
    "gain:f" => \$gain_cutoff,             # Default: 2.7
    "loss:f" => \$loss_cutoff,             # Default: 1.3
    "tnum:i" => \$target_num_cutoff,       # Default: 3
    "pct:f" => \$cnv_pct_cutoff,           # Default: 70%
    "o:s" => \$outfile,                    # Need
    ) or die "unknown args\n";

# input file: *.log2Ratio.xls
# output file: *.Gene_Level_CNV.xls

# filter rules:
# 1. skip gene with <= 3 targets
# 2. skip target with low depth (<=30X)


# default value
if (not defined $gain_cutoff){
    $gain_cutoff = 2.7;
}

if (not defined $loss_cutoff){
    $loss_cutoff = 1.3;
}

if (not defined $target_num_cutoff){
    $target_num_cutoff = 3;
}

if (not defined $cnv_pct_cutoff){
    $cnv_pct_cutoff = 70;
}



open O, ">$outfile" or die;
print O "sample\tgene\tchr\tstart\tend\tcopynumber\tcnvtype\n";


my %geneStartEnd;
my %geneChr;
my %geneCN;
my @genes;
my @sampleName;
open IN, "$logR" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $cn;
    if ($arr[-2] != 0){
        $cn = sprintf "%.2f", $arr[-3]/$arr[-2]*2;
    }else{
        $cn = "NA";
    }
    
    if ($arr[5] <= 30 || $cn eq "NA"){
        push @{$geneCN{$arr[4]}{"fail"}}, $cn;
    }else{
        push @{$geneCN{$arr[4]}{"ok"}}, $cn;
    }

    $geneChr{$arr[4]} = $arr[1];
    push @{$geneStartEnd{$arr[4]}}, $arr[2];
    push @{$geneStartEnd{$arr[4]}}, $arr[3];
    push @genes, $arr[4];
    push @sampleName, $arr[0];
}
close IN;

my %geneRegion;
my %geneFlag;
my @uniqGenes;

for my $gene (@genes){
    if (!exists $geneFlag{$gene}){
        $geneFlag{$gene} = 1;
        my @reg = sort {$a <=> $b} @{$geneStartEnd{$gene}};
        $geneRegion{$gene}{"start"} = $reg[0];
        $geneRegion{$gene}{"end"} = $reg[-1];
        push @uniqGenes, $gene;
    }else{
        next;
    }
}

for my $gene (@uniqGenes){
    #print "$gene\n";
    my @ok_target_cn;
    if (exists $geneCN{$gene}{"ok"}){
        @ok_target_cn = @{$geneCN{$gene}{"ok"}};
    }else{
        @ok_target_cn = qw//;
    }

    my @fail_target_cn;
    if (exists $geneCN{$gene}{"fail"}){
        @fail_target_cn = @{$geneCN{$gene}{"fail"}};
    }else{
        @fail_target_cn = qw//;
    }

    my $target_num = scalar(@ok_target_cn) + scalar(@fail_target_cn);
    #print "$target_num\n";
    my $ok_pct = sprintf "%.2f", scalar(@ok_target_cn)/$target_num * 100;
    #print "$ok_pct\n";
    if ($target_num <= 3 || $ok_pct <= 70){
        next; # filter low quality gene or gene with less targets
    }

    my $cnvtype;
    my $gain_num = 0;
    my $loss_num = 0;
    for my $cn (@ok_target_cn){
        if ($cn >= $gain_cutoff){
            $gain_num += 1;
        }
        if ($cn <= $loss_cutoff){
            $loss_num += 1;
        }
    }

    my ($gain_pct,$loss_pct) = (0,0);
    $gain_pct = sprintf "%.2f", $gain_num/$target_num * 100;
    $loss_pct = sprintf "%.2f", $loss_num/$target_num * 100;


    my ($sum_cn,$mean_cn) = (0,0);

    if ($gain_pct >= $cnv_pct_cutoff){
        # re-check cn
        for my $cn (@ok_target_cn){
            $sum_cn += $cn;
        }
        
        $mean_cn = sprintf "%.2f", $sum_cn/scalar(@ok_target_cn);
        
        if ($mean_cn >= $gain_cutoff){
            $cnvtype = "gain";
        }else{
            $cnvtype = "normal";
        }
    }elsif ($loss_pct >= $cnv_pct_cutoff){
        for my $cn (@ok_target_cn){
            $sum_cn += $cn;
        }

        $mean_cn = sprintf "%.2f", $sum_cn/scalar(@ok_target_cn);

        if ($mean_cn <= $loss_cutoff){
            $cnvtype = "loss";
        }else{
            $cnvtype = "normal";
        }
    }else{
        for my $cn (@ok_target_cn){
            $sum_cn += $cn;
        }

        $mean_cn = sprintf "%.2f", $sum_cn/scalar(@ok_target_cn);

        $cnvtype = "normal";
    }

    my $start = $geneRegion{$gene}{"start"};
    my $end = $geneRegion{$gene}{"end"};

    print O "$sampleName[0]\t$gene\t$geneChr{$gene}\t$start\t$end\t$mean_cn\t$cnvtype\n";
}

close O;

