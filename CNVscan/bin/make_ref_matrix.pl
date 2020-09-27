use strict;
use warnings;
use File::Basename;

my ($ref_dir,$outfile) = @ARGV;

my @norm = glob "$ref_dir/*.norm.xls";
my $ref_num = scalar(@norm);

print "find $ref_num ref files (*.norm.xls)\n";

my @samples;
for my $s (@norm){
    my $name = (split /\./, basename $s)[0];
    push @samples, $name;
}

open O, ">$outfile" or die;
print O "target";

# header
for my $s (@samples){
    print O "\t$s";
}
print O "\n";


# get target
my @targets;
my $this_sample = $norm[0];
open IN, "$this_sample" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $t = "$arr[0]\.$arr[1]\-$arr[2]\.$arr[3]"; # chr.start-end.gene
    push @targets, $t;
}
close IN;


#
my %ref_norm;
for my $s (@norm){
    print "process $s\n";
    open IN, "$s" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        my $t = "$arr[0]\.$arr[1]\-$arr[2]\.$arr[3]"; # chr.start-end.gene
        push @{$ref_norm{$t}}, $arr[-1];
    }
    close IN;
}

# output
for my $t (@targets){
    my @norm_val = @{$ref_norm{$t}};
    print O "$t";
    for my $v (@norm_val){
        print O "\t$v";
    }
    print O "\n";
}
close O;

