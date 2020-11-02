import sys
import os
import configparser
import argparse
import glob

def parse_args():
    AP = argparse.ArgumentParser("detect HLA LOH by BAF using paired tumor-normal data")
    AP.add_argument('-tbam',help='tumor bam file',dest='tbam')
    AP.add_argument('-nbam',help='normal bam file',dest='nbam')
    AP.add_argument('-tname',help='tumor sample name',dest='tname')
    AP.add_argument('-nname',help='tumor sample name',dest='nname')
    AP.add_argument('-fa',help='fasta file',dest='fasta',default='/data1/database/b37/human_g1k_v37.fasta')
    #AP.add_argument('-p',help='panel, can be <889|338>',dest='panel')
    AP.add_argument('-snpBED',help='snp bed file',dest='snpBED',default='/data1/workdir/fulongfei/git_repo/HLALOH/BAF/889gene.snp.bed')
    AP.add_argument('-py2',help='python2 path',dest='py2',default='/home/fulongfei/miniconda3/envs/py27/bin/python2')
    AP.add_argument('-py3',help='python3 path',dest='py3',default='/home/fulongfei/miniconda3/bin/python3')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()


#py3 = '/home/fulongfei/miniconda3/bin/python3'
#py2 = '/home/fulongfei/miniconda3/envs/py27/bin/python2'

def main():
    '''
    main steps:
        1. hla typing
        2. get each allele's fasta sequence
        3. get het pos by pairwise alignment
        4. extract HLA reads from tumor bam
        5. align HLA reads using hla ref fasta
        6. make pileup for each allele and cal BAF
        7. plot HLA BAF fig
    '''

    args = parse_args()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]

    runsh = "%s/%s.HLALOH.sh" % (args.outdir,args.tname)
    of = open(runsh,'w')

    # cal snp vaf for tumor to visually estimate tumor's purity and ploidy
    of.write("###cal tumor snp vaf"+'\n')
    
    '''
    if args.panel:
        # if panel exists
        if args.panel == '889':
            snp_bed = "%s/BAF/889gene.snp.bed" % (bin_dir)
        elif args.panel == '338':
            snp_bed = "%s/BAF/338gene.snp.bed" % (bin_dir)
        else:
            snp_bed = "NA"
    else:
        # if panel not exists
        snp_bed = args.snp_bed
    '''

    snp_bed = args.snpBED
    snp_vaf = "%s/tumor.snp.vaf" % (args.outdir)
    tumor_pileup = "%s/tumor.mpileup" % (args.outdir)
    cmd = "samtools mpileup -d 8000 -f %s -l %s %s >%s" % (args.fasta,snp_bed,args.tbam,tumor_pileup)
    of.write(cmd+'\n')
    cmd = "%s %s/BAF/pileup2vaf.py %s %s" % (args.py3,bin_dir,tumor_pileup,snp_vaf)
    of.write(cmd+'\n')

    # plot fig
    cmd = "Rscript %s/BAF/plot_vaf_by_chr.r %s %s %s" % (bin_dir,snp_vaf,args.tname,args.outdir)
    of.write(cmd+'\n\n')





    # extract HLA region read1/read2 from normal bam to HLA typing
    of.write("###extract HLA reads"+'\n')
    chrname = chrNaming(args.fasta)
    extract_HLA_reads(args.nbam,args.nname,args.outdir,of,chrname)

    # HLA typing using OptiType
    fq1 = "%s/%s.chr6region.1.fastq" % (args.outdir,args.nname)
    fq2 = "%s/%s.chr6region.2.fastq" % (args.outdir,args.nname)

    cmd = "%s %s/OptiType-1.3.2/OptiTypePipeline.py -i %s %s --dna -v -o %s -c %s/OptiType-1.3.2/config.ini.example" % (args.py2,
        bin_dir,
        fq1,
        fq2,
        args.outdir,
        bin_dir
        )

    of.write('\n'+"###HLA Typing"+'\n')
    of.write(cmd+'\n')

    # transformat HLA result
    cmd = "cp %s/*/*_result.tsv %s/hla.result.raw" % (args.outdir,args.outdir)
    of.write(cmd+'\n')

    #hla_res = glob.glob("%s/*/*_result.tsv" % args.outdir)[0]
    #print(hla_res)
    raw_hla_res = "%s/hla.result.raw" % (args.outdir)
    new_hla_res = "%s/hla.result.new" % (args.outdir)
    cmd = "%s %s/tools/translate_HLA_format.py -i %s -o %s" % (args.py3,bin_dir,raw_hla_res,new_hla_res)

    of.write(cmd+'\n')


    # get HLA alleles' fasta sequence
    hla_fa = "%s/patient.hlaFasta.fa" % (args.outdir)
    cmd = "%s %s/tools/get_hla_fasta.py -i %s -o %s" % (args.py3,bin_dir,raw_hla_res,hla_fa)

    of.write('\n'+"###get hla fasta"+'\n')
    of.write(cmd+'\n')


    # extract HLA reads from tumor
    of.write('\n'+"###extract HLA reads from tumor BAM file"+'\n')
    extract_HLA_reads(args.tbam,args.tname,args.outdir,of,chrname)

    fq1_tumor = "%s/%s.chr6region.1.fastq" % (args.outdir,args.tname)
    fq2_tumor = "%s/%s.chr6region.2.fastq" % (args.outdir,args.tname)

    # aln fq using hla fasta as ref
    of.write('\n'+'###aln fq to hla ref fasta using novoalign'+'\n')
    hla_sam = "%s/%s.chr6region.patient.reference.hlas.sam" % (args.outdir,args.tname)
    novoIndexFile = "%s/patient.hlaFasta.nix" % (args.outdir)
    novoindexBinary = "%s/novoDir/novocraft/novoindex" % (bin_dir)
    novoalignBinary = "%s/novoDir/novocraft/novoalign" % (bin_dir)
    cmd = "%s %s %s" % (novoindexBinary,novoIndexFile,hla_fa)
    of.write(cmd+'\n')

    cmd = "%s -d %s -f %s %s -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW >%s" % (
                                                                                novoalignBinary,
                                                                                novoIndexFile,
                                                                                fq1_tumor,
                                                                                fq2_tumor,
                                                                                hla_sam
                                                                                        )

    of.write(cmd+'\n')

    # convert sam->bam
    hla_bam = "%s/%s.chr6region.patient.reference.hlas.bam" % (args.outdir,args.tname)
    cmd = "samtools view -b -o %s %s" % (hla_bam,hla_sam)
    of.write(cmd+'\n')

    # sort
    hla_sort_bam = "%s/%s.chr6region.patient.reference.hlas.csort.bam" % (args.outdir,args.tname)
    cmd = "samtools sort -@ 8 -o %s %s" % (hla_sort_bam,hla_bam)
    of.write(cmd+'\n')

    # rmdup
    hla_sort_rmdup_bam = "%s/%s.chr6region.patient.reference.hlas.csort.noduplicates.bam" % (args.outdir,args.tname)
    cmd = "samtools rmdup %s %s" % (hla_sort_bam,hla_sort_rmdup_bam)
    of.write(cmd+'\n')

    # only take reads that are in proper pair
    proper_align_bam = "%s/%s.chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam" % (args.outdir,args.tname)
    cmd = "samtools view -f 2 -b -o %s %s" % (proper_align_bam,hla_sort_rmdup_bam)
    of.write(cmd+'\n')

    # let's index the aligned bam
    cmd = "samtools index %s" % (proper_align_bam)
    of.write(cmd+'\n')

    of.write('\n'+'###baf2loh main script'+'\n')
    cmd = "%s %s/BAF2LOH.py -indir %s -tname %s" % (args.py3,bin_dir,args.outdir,args.tname)
    of.write(cmd+'\n')

    of.close()



def chrNaming(ref):
    '''
    check chr naming
    UCSC's chr naing is "chrN"
    '''
    infile = open(ref,'r')
    firstLine = infile.readline()
    if '>chr' in firstLine:
        chrNaming = 'with_prefix'
    else:
        chrNaming = 'no_prefix'

    return(chrNaming)

def extract_HLA_reads(bam,name,outdir,runshFH,chrNaming):
    cmd = "samtools view -H %s >%s/%s.hla.sam" % (
                                            bam,
                                            outdir,
                                            name
                                            )

    runshFH.write(cmd+'\n')

    if chrNaming == 'no_prefix':
        region1 = '6:29909037-29913661'
        region2 = '6:31321649-31324964'
        region3 = '6:31236526-31239869'
    else:
        region1 = 'chr6:29909037-29913661'
        region2 = 'chr6:31321649-31324964'
        region3 = 'chr6:31236526-31239869'

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (
                                            bam,
                                            region1,
                                            outdir,
                                            name
                                            )

    
    runshFH.write(cmd+'\n')
    #os.system(cmd)

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (
                                                bam,
                                                region2,
                                                outdir,
                                                name
                                                )

    runshFH.write(cmd+'\n')
    #os.system(cmd)

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (
                                                bam,
                                                region3,
                                                outdir,
                                                name
                                                )

    runshFH.write(cmd+'\n')

    
    ### ignore chr6 contig

    #1.chr6_apd_hap1
    #2.chr6_cox_hap2
    #3.chr6_dbb_hap3
    #4.chr6_mann_hap4
    #5.chr6_mcf_hap5
    #6.chr6_qbl_hap6
    #7.chr6_ssto_hap7
    
    
    # turn into fastq
    fq1 = "%s/%s.chr6region.1.fastq" % (outdir,name)
    fq2 = "%s/%s.chr6region.2.fastq" % (outdir,name)

    hla_sam = "%s/%s.hla.sam" % (outdir,name)

    # convert to bam
    hla_bam = "%s/%s.hla.bam" % (outdir,name)
    cmd = "samtools view -b -o %s %s" % (hla_bam,hla_sam)
    runshFH.write(cmd+'\n')

    # sort by name
    bam_sort_by_name = "%s/%s.sort_by_name.bam" % (outdir,name)
    cmd = "samtools sort -n %s -o %s" % (hla_bam,bam_sort_by_name)
    runshFH.write(cmd+'\n')

    # remove temp sam/bam

    # bedtools bamtofastq
    cmd = "bedtools bamtofastq -i %s -fq %s -fq2 %s" % (bam_sort_by_name,fq1,fq2)
    runshFH.write(cmd+'\n')
    

    # Illegal Mate State
    # https://www.biostars.org/p/59521/

    #cmd = "%s -jar %s/SamToFastq.jar I=%s F=%s F2=%s VALIDATION_STRINGENCY=SILENT" % (java, gatk_dir, hla_sam, fq1, fq2)


    
if __name__ == "__main__":
    main()
