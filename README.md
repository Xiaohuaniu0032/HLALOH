# HLALOH
detect HLA LOH from NGS data

## Usage
`python3 predictHLALOH.py -nbam <normal.bam> -tbam <tumor.bam> -bed <bed file> -nname <normal name> -tname <tumor name> -fa <hg19 fasta> -ref <cnv control dir> [-p <889|338> | -snp <snp bed file>] -od <outdir>`


```
$ python predictHLALOH.py -h
usage: detect HLA LOH from capture NGS data using paired tumor/normal
       [-h] [-nbam NBAM] [-tbam TBAM] [-bed BED] [-nname NNAME] [-tname TNAME]
       [-fa FASTA] [-ref REF] [-p PANEL] [-snp SNP_BED] [-od OUTDIR]

optional arguments:
  -h, --help    show this help message and exit
  -nbam NBAM    normal bam file
  -tbam TBAM    tumor bam file
  -bed BED      bed file
  -nname NNAME  normal sample name
  -tname TNAME  tumor sample name
  -fa FASTA     fasta file
  -ref REF      cnv ref control dir
  -p PANEL      panel, can be <889|338>
  -snp SNP_BED  snp bed
  -od OUTDIR    output dir

```

### parameter specification
`-nbam`: normal bam

`-tbam`: tumor bam

`-bed`: bed file

`-nname`: normal sample name

`-tname`: tumor sample name

`fa`: hg19 fasta

`-ref`: cnv control dir

`-p`: panel, this argument specified the `snp bed file`

`-snp`: snp bed

`-od`: output dir


***NOTE:***
you can specifiy `-p` or `-snp`, `-p` will be considered first. `-p` can give you a snp bed file (this bed file is included in this git repo) by `-p`.

## Method
*we will follow this paper's method:* `Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell, 2017`

This paper published an open-source software called `lohhla`, you can see official doc `https://bitbucket.org/mcgranahanlab/lohhla/src/master/`

*the main steps of `lohhla` are:*

1. extract HLA reads
2. create HLA allele specific BAM files
3. determine coverage at mismatch positions between homologous HLA alleles
4. obtain HLA specific logR and BAF
5. determine HLA haplotype specific copy number

the `lohhla` software needs several files as input, which means that you need to prepare these files before you run `lohhla`):

* hla typing result (we use Optitype software)
* tumor purity and ploidy (we use ASCAT software)

the `ASCAT` software needs below files:

* logR (we use self-developed software: CNVscan)
* BAF (we use self-developed python script)


## Testing
1. `cd /path/HLALOH/test`
1. `sh run.sh`
2. `sh 20091701T.HLALOH.sh`

you can find three sub-dirs under `/path/HLALOH/test` dir:

* bamdir (contain normal and tumor bams' soft link)
* lohhla (lohhla software's results)
* 2020\_10\_14\_15\_30\_28 (created by Optitype)

and some files.

## Software Needed
1. lohhla
    * novoindex (included by this git repo, you do not need to install it)
    * Jellyfish (https://www.cbcb.umd.edu/software/jellyfish/, a fast k-mer counting tool for DNA, you can use conda to install it. `conda install -c bioconda jellyfish`. you need to add binary path into your $PATH)
    * bedtools (you need to add binary path into your $PATH)
    * samtools (you need to add binary path into your $PATH)
    * picard (SortSam.jar,FilterSamReads.jar. included by this git repo, you do not need to install it.)
    * R
2. ASCAT
    * R package: ASCAT
3. CNVscan
    * see `/path/HLALOH/CNVscan/config.ini`
4. OptiType
    * razers3 (you need to install it and modify `/path/HLALOH/OptiType-1.3.2/config.ini.example`)
    * glpk (you need to install glpk and add binary path into your $PATH)

## 
