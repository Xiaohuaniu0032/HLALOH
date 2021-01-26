# HLALOH
detect HLA LOH from NGS data

## Usage
`python3 predictHLALOH.v3.py -nbam <normal.bam> -tbam <tumor.bam> -nname <normal_name> -tname <tumor_name> -r <hg19.fasta> -cnvRefDir <Dir> -bed <*.BED> -py2 <python2> -py3 <python3> -sbb <sambamba> -rscript <Rscript> -bedtools <BEDTOOLS> -bwa <BWA> -samts <SAMTOOLS> -jre <JRE> -od <outputDir>`

```
$ python3 predictHLALOH.v3.py -h
usage: detect HLA LOH by BAF using paired tumor-normal data [-h] [-tbam TBAM] [-nbam NBAM] [-tname TNAME]
                                                            [-nname NNAME] [-r FASTA] [-cnvRefDir CNVREFDIR]
                                                            [-bed BED] [-py2 PY2] [-py3 PY3] [-sbb SBB]
                                                            [-rscript RSCRIPT] [-bedtools BEDTOOLS] [-bwa BWA]
                                                            [-samts SAMTOOLS] [-jre JRE] [-od OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -tbam TBAM            tumor bam file
  -nbam NBAM            normal bam file
  -tname TNAME          tumor sample name
  -nname NNAME          tumor sample name
  -r FASTA              fasta file
  -cnvRefDir CNVREFDIR  cnv ref control dir
  -bed BED              bed file
  -py2 PY2              python2 path
  -py3 PY3              python3 path
  -sbb SBB              sambamba path
  -rscript RSCRIPT      R path
  -bedtools BEDTOOLS    bedtools path
  -bwa BWA              bwa path
  -samts SAMTOOLS       samtools bin
  -jre JRE              java JRE
  -od OUTDIR            output dir

```
## Method
For WES data, there exists a published tool `lohhla` for HLA LOH detection, but for small and medium-size capture panel data, we find that the `lohhla` is not very appropriate because `lohhla` needs to estimate tumor purity and ploidy, which is hard to estimate from panel data. besides the tumor purity and ploidy, we also find that the results output by `lohhla` are a little confused.

We just use the BAF distribution to infer the HLA LOH status directly, although the method is simple and intuitional, but the result is effective (see below fig).

![tumor purity & BAF distribution](https://github.com/Xiaohuaniu0032/HLALOH/blob/master/qpure.plos.one.png)

For detail of `lohhla`, you can see `Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell, 2017` and official doc `https://bitbucket.org/mcgranahanlab/lohhla/src/master/`

*the main steps of our methods are:*

1. hla typing using normal bam
2. get six hla alleles from IMGT/HLA database to make patient-specific HLA reference
3. extract hla region reads from tumor bam
4. align tumor fastq to patient-specific HLA ref
5. get allele-specific bam
6. filter allele-specific bam (mismatch + ins + del need to <= 1)
7. for A/B/C allele, determine het positions
8. calculate BAF at those het positions (for normal and tumor sample)
9. calculate normal and tumor sample's copy number (copy number information is not used to infer LoH status at present, but the pipeline still need to call CNV at present, please note this)
10. infer LoH status by BAF distribution at het positions (first filter out abnormal het positions in normal sample, then infer LoH status using left het positions BAF)

## How to determine LoH by BAF distribution?
after we get filtered tumor BAF file, we use two information to infer LoH status:

1. likely loh het pos number and total het pos

2. median BAF

let us mark `likely loh het pos number` as `n`, and mark `total het pos` as `N`, then we can calculate loh pos percentage: `n/N`.

For a specific HLA allele, if `n/N` >= 0.65 and `median_baf` >= 0.6, we think this allele in tumor is in LoH status, otherwise not in LoH status.


In theory, if one allele is lost in tumor and the tumor's ploidy is 2 and purity is 100%, then for a het position in normal, the BAF of this position in tumor will be 0 or 100%. if the purity is not 100%, then the BAF will be centered around 0.5 (for example, if tumor ploidy is 2 and purity is 70%, then then BAF will be 0.77/0.23)

when the tumor purity is 30% (assume ploidy is 2), then BAF will be 0.59/0.41.

if the BAF is < 0.4 or > 0.6, then we think this het pos is a likely LoH position (or you can say that this het pos is likely covered by a LoH region)

so, the cutoff we used here is based on the assumption that we can detect the LoH as long as the tumor purity >= 30%. for the tumor with purity < 30%, it is hard to detect the LoH with hgih confidence.



## Testing
1. `cd /path/HLALOH/test` you will see two directoires: `Negative_Sample` and `Positive_Sample`
2. `Negative_Sample` is a HLA-A LoH negative sample (validate by STR + GeneMapper), `Positive_Sample` is a HLA-A LoH positive sample (validate by STR + GeneMapper)
3. cd `Positive_Sample` and `sh run.sh` (please modify `run.sh`)
4. `nohup sh 201023033T.HLALOH.sh &`

The final result file is `*.hlaloh.result.xls`, the content in this file is:

```
Sample	Gene	HetPosNum	LikelyLoHPosNum	LoHPosPercent(%)	MedianBAF	CopyNumber	IfLoH
201023033T	HLA-A	17	17	1.00	0.73	2.76	YES
201023033T	HLA-B	12	12	1.00	0.72	2.80	YES
201023033T	HLA-C	21	19	0.90	0.63	2.81	YES
```

1th col (Sample): sample name

2th col (Gene): gene

3th col (HetPosNum): the number of het positions in tumor BAF file

4th col (LikelyLoHPosNum): the number of likely LoH positions

5th col (LoHPosPercent(%)): LikelyLoHPosNum/HetPosNum

6th col (MedianBAF): median BAF (for BAF, its value is mirror symmetry at 0.5, so if BAF < 0.5, we will use 1-BAF to calculate median BAF)

7th (CopyNumber): copy number (not used to infer LoH, and we find the copy number for HLA region is not accurate)

8th (IfLoH): YES/NO

*Note:* if an allele has <= 10 het pos, we will not infer the LoH, so the IfLoH col will be `NA`, please note this.



## Software Needed
1. OptiType (used for HLA typing, see `https://github.com/FRED-2/OptiType` for install)
2. CNVscan (used for calculate CNV, see `https://github.com/Xiaohuaniu0032/CNVscan` for usage)
3. python2
4. python3
5. sambamba
6. Rscript
7. Bedtools
8. BWA
9. Samtools
10. Java JRE


## Questions
###1) can be used for WES data?
yes

###2) can be used for panel data?
yes, make sure your panel capures HLA region

###3) how to validate HLA LoH?
you can see `lohhla` Cell paper. it use polymorphic STR that close to HLA-A/B/C to detect allele imbalance.
