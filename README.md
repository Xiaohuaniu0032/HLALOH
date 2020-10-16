# HLALOH
detect HLA LOH from NGS data

## Usage
`python3 predictHLALOH.py -nbam <normal.bam> -tbam <tumor.bam> -bamDir <DIR> -bed <*.bed> -nname <normal.name> -tname <tumor.name> -fa <*.fasta> -ref <cnv control dir> -gatkDir <DIR> -novoDir <DIR> -od <outdir>`

```
$ python predictHLALOH.py -h
usage: detect HLA LOH from capture NGS data using paired tumor/normal
       [-h] [-nbam NBAM] [-tbam TBAM] [-bamDir BAMDIR] [-bed BED]
       [-nname NNAME] [-tname TNAME] [-fa FASTA] [-ref REF] [-gatkDir GATKDIR]
       [-novoDir NOVODIR] [-od OUTDIR]

optional arguments:
  -h, --help        show this help message and exit
  -nbam NBAM        normal bam file
  -tbam TBAM        tumor bam file
  -bamDir BAMDIR    dir contains tumor and normal bam files
  -bed BED          bed file
  -nname NNAME      normal sample name
  -tname TNAME      tumor sample name
  -fa FASTA         fasta file
  -ref REF          cnv ref control dir
  -gatkDir GATKDIR  gatk dir
  -novoDir NOVODIR  novoalign bin dir
  -od OUTDIR        output dir
```

## Method
*we will follow this paper's method:* `Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell, 2017`

*you can see `METHOD DETAILS` for more details*

*the main steps are:*

1. extract HLA reads
2. create HLA allele specific BAM files
3. determine coverage at mismatch positions between homologous HLA alleles
4. obtain HLA specific logR and BAF
5. determine HLA haplotype specific copy number


**there are several main step for HLA LOH:**

1. HLA Class-I A/B/C alleles inference (we use Optitype software)
2. CNV calling (we use self-developed software: CNVscan)
3. BAF calculate (we use self-developed python script)
4. tumor purity and ploidy estimate (we use ASCAT software, which needs tumor and normal bam files as input files)
5. HLA LOH inference (we use `Cell,2017` paper's software: lohhla)


## Testing

## Software Needed
1. lohhla
    * novoindex
    * Jellyfish (https://www.cbcb.umd.edu/software/jellyfish/, a fast k-mer counting tool for DNA)
    `conda install -c bioconda jellyfish`
    `which jellyfish`
    * 
