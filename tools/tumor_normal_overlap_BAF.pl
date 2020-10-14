use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($normal_vaf,$tumor_vaf,$outdir);
GetOptions(
    "nvaf:s" => \$normal_vaf,
    "tvaf:s" => \$tumor_vaf,
    "od:s" => \$outdir,
    ) or die "unknown args\n";

my $n_name = (split /\./, basename $normal_vaf)[0];
my $t_name = (split /\./, basename $tumor_vaf)[0];

my $out_n_vaf = "$outdir/$n_name\.normal.overlap.vaf";
my $out_t_vaf = "$outdir/$t_name\.tumor.overlap.vaf";

print "will creat $out_n_vaf file\n";
print "will creat $out_t_vaf file\n";


my %all_vaf;
open IN, "$normal_vaf" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }
    $all_vaf{"n"}{$chr}{$arr[1]} = 1;
}
close IN;

open IN, "$tumor_vaf" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }
    $all_vaf{"t"}{$chr}{$arr[1]} = 1;
}
close IN;







open O1, ">$out_n_vaf" or die;
my $h = "chr\tpos\tref_base\tref_n\talt_n\tall_n\talt_freq\n";
print O1 "$h";
open IN, "$normal_vaf" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }
    
    if (exists $all_vaf{"t"}{$chr}{$arr[1]}){
        # occur in tumor.vaf
        print O1 "$_\n";
    }
}
close IN;
close O1;





open O2, ">$out_t_vaf" or die;
print O2 "$h";
open IN, "$tumor_vaf" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $chr = $arr[0];
    if ($chr =~ /^chr/){
        $chr =~ s/^chr//;
    }

    if (exists $all_vaf{"n"}{$chr}{$arr[1]}){
        # occur in normal.vaf
        print O2 "$_\n";
    }
}
close IN;
close O2;
