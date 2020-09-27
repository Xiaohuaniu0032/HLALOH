# CNVscan
detect somatic copy number from capture NGS data using Pool of Normal (PoN) method

## Usage
`python3 CNVscan.py -bam <*.bam> -bed <*.bed> -n <name> -fa <*.fasta> -m <ref|cnv> -od <outdir>`



```
$ python CNVscan.py -h
usage: detect somatic copy number from capture NGS data [-h] [-bam BAM]
                                                        [-bed BED] [-n NAME]
                                                        [-fa FASTA] [-m MODE]
                                                        [-ref REF]
                                                        [-od OUTDIR]

optional arguments:
  -h, --help  show this help message and exit
  -bam BAM    bam file
  -bed BED    bed file
  -n NAME     sample name
  -fa FASTA   fasta file
  -m MODE     analysis mode. can be [ref|cnv]
  -ref REF    control dir
  -od OUTDIR  out dir
```

### *Note:*
please change `config.ini` file when you run this pipeline

```
$ cat config.ini
[software]
perl = /home/fulongfei/miniconda3/bin/perl
python3 = /home/fulongfei/miniconda3/bin/python3
bedtools = /home/fulongfei/miniconda3/bin/bedtools
sambamba = /home/fulongfei/miniconda3/bin/sambamba
Rscript = /home/fulongfei/miniconda3/bin/Rscript


[db]
hg19 = /data1/database/b37/human_g1k_v37.fasta
```


## Method
1. calculate region coverage
2. gc correction
3. library size normalization
4. calculate logR for each region
5. calculate gene-level copy number


### parameter specification
`-m`: can be `ref` or `cnv`

`-ref`: just get the normalized depth file (*.norm.xls)

`-cnv`: calculate copy number using PoN method

### how to construct PoN control?
1. get a set of wbc sample's bam file (30~50 samples is ok)
2. use `-m ref` to generate a set of `*.norm.xls` files
3. copy all of `*.norm.xls` files into a directory
4. call cnv using `-m cnv` (please specify `-ref`)


## Testing

1. `sh /path/CNVscan/test/run.sh`
2. `cd /path/CNVscan/test/BB20062440_CL01167`
3. `sh ref.BB20062440_CL01167.sh >ref.log 2>&1 &`






