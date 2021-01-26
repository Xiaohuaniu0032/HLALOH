import sys
import os
import configparser
import argparse
import glob

import tools.chr_naming
import tools.extract_HLA_reads

def parse_args():
    AP = argparse.ArgumentParser("detect HLA LOH by BAF using paired tumor-normal data")
    AP.add_argument('-tbam',help='tumor bam file',dest='tbam')
    AP.add_argument('-nbam',help='normal bam file',dest='nbam')
    AP.add_argument('-tname',help='tumor sample name',dest='tname')
    AP.add_argument('-nname',help='tumor sample name',dest='nname')
    AP.add_argument('-r',help='fasta file',dest='fasta',default='/data1/database/b37/human_g1k_v37.fasta')
    AP.add_argument('-cnvRefDir',help='cnv ref control dir',dest='cnvRefDir',default='/data1/workdir/fulongfei/git_repo/HLALOH/DB/cnv_ref/889/DB')
    AP.add_argument('-bed',help='bed file',dest='bed',default='/home/wangce/workdir/database/humandb/panel/889genes_20191225.bed')
    #AP.add_argument('-snpBED',help='snp bed file',dest='snpBED',default='/data1/workdir/fulongfei/git_repo/HLALOH/BAF/889gene.snp.bed')
    AP.add_argument('-py2',help='python2 path',dest='py2',default='/home/fulongfei/miniconda3/envs/py27/bin/python2')
    AP.add_argument('-py3',help='python3 path',dest='py3',default='/home/fulongfei/miniconda3/bin/python3')
    AP.add_argument('-sbb',help='sambamba path',dest='sbb',default='/home/fulongfei/miniconda3/bin/sambamba')
    AP.add_argument('-rscript',help='R path',dest='rscript',default='/home/fulongfei/miniconda3/bin/Rscript')
    AP.add_argument('-bedtools',help='bedtools path',dest='bedtools',default='/home/fulongfei/miniconda3/bin/bedtools')
    AP.add_argument('-bwa',help='bwa path',dest='bwa',default='/usr/bin/bwa')
    AP.add_argument('-samts',help='samtools bin',dest='samtools',default='/home/fulongfei/miniconda3/bin/samtools')
    AP.add_argument('-jre',help='java JRE',dest='jre',default='/usr/bin/java')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()


#py3 = '/home/fulongfei/miniconda3/bin/python3'
#py2 = '/home/fulongfei/miniconda3/envs/py27/bin/python2'

def main():
    args = parse_args()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]

    runsh = "%s/%s.HLALOH.sh" % (args.outdir,args.tname)
    of = open(runsh,'w')

    # cal snp vaf for tumor to visually estimate tumor's purity and ploidy
    #of.write("###cal tumor snp vaf"+'\n')
    #tumor_pileup = "%s/tumor.mpileup" % (args.outdir)
    #cmd = "%s mpileup -x -Q 0 -d 8000 -f %s -l %s %s >%s" % (args.samtools,args.fasta,args.snpBED,args.tbam,tumor_pileup) # -x will output right baseQ, -Q will not skip any low baseQ
    #of.write(cmd+'\n')

    # pileup to vaf
    #snp_vaf = "%s/tumor.snp.vaf" % (args.outdir)
    #cmd = "%s %s/BAF/pileup2vaf.py %s %s" % (args.py3,bin_dir,tumor_pileup,snp_vaf)
    #of.write(cmd+'\n')

    # plot fig
    #cmd = "%s %s/BAF/plot_vaf_by_chr.r %s %s %s" % (args.rscript,bin_dir,snp_vaf,args.tname,args.outdir)
    #of.write(cmd+'\n\n')


    ############################# MAIN STEP ###############################

    # extract HLA region read1/read2 from normal bam to do HLA typing
    of.write("###extract HLA reads"+'\n')
    chrname = tools.chr_naming.chrNaming(args.fasta)
    tools.extract_HLA_reads.extract_HLA_reads(args.nbam,args.nname,args.outdir,of,chrname)
    

    # HLA typing using OptiType
    fq1 = "%s/%s.chr6region.1.fastq" % (args.outdir,args.nname)
    fq2 = "%s/%s.chr6region.2.fastq" % (args.outdir,args.nname)

    cmd = "%s %s/OptiType-1.3.2/OptiTypePipeline.py -i %s %s --dna -v -o %s -c %s/OptiType-1.3.2/config.ini" % (
        args.py2,
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

    raw_hla_res = "%s/hla.result.raw" % (args.outdir)
    new_hla_res = "%s/hla.result.new" % (args.outdir)
    cmd = "%s %s/tools/translate_HLA_format.py -i %s -o %s" % (args.py3,bin_dir,raw_hla_res,new_hla_res)
    of.write(cmd+'\n')

    # get HLA allele fa seq
    hla_fa = "%s/patient.hlaFasta.fa" % (args.outdir)
    hla_fa_all = "%s/patient.hlaFasta.all.4digit.alleles.fa" % (args.outdir) # not used for the align info is so complex
    cmd = "%s %s/tools/get_hla_fasta.py -i %s -o %s -o2 %s" % (args.py3,bin_dir,raw_hla_res,hla_fa,hla_fa_all)
    of.write('\n'+"###get hla fasta"+'\n')
    of.write(cmd+'\n')



    ############### extract HLA reads from normal
    of.write('\n'+"###extract HLA reads from normal BAM file"+'\n')
    #tools.extract_HLA_reads.extract_HLA_reads(args.nbam,args.nname,args.outdir,of,chrname)
    
    # aln fq and get proper_aln reads
    fq1_normal = "%s/%s.chr6region.1.fastq" % (args.outdir,args.nname)
    fq2_normal = "%s/%s.chr6region.2.fastq" % (args.outdir,args.nname)

    cmd = "%s %s/tools/get_proper_aln_bam.py %s %s %s %s %s %s %s %s" % (args.py3,bin_dir,fq1_normal,fq2_normal,args.nname,hla_fa,args.bwa,args.samtools,args.sbb,args.outdir)
    of.write(cmd+'\n')

    # aln two allele and get het pos and cal BAF
    of.write('\n'+'###baf2loh main script for normal'+'\n')
    normal_proper_aln_bam = "%s/%s.chr6region.patient.reference.hlas.sort.markdup.proper_aln.bam" % (args.outdir,args.nname)
    cmd = "%s %s/BAF2LOH.py -indir %s -name %s -bam %s -jre %s" % (args.py3,bin_dir,args.outdir,args.nname,normal_proper_aln_bam,args.jre)
    of.write(cmd+'\n')


    ############### extract HLA reads from tumor
    of.write('\n'+"###extract HLA reads from tumor BAM file"+'\n')
    tools.extract_HLA_reads.extract_HLA_reads(args.tbam,args.tname,args.outdir,of,chrname)

    # aln fq and get proper_aln reads
    fq1_tumor = "%s/%s.chr6region.1.fastq" % (args.outdir,args.tname)
    fq2_tumor = "%s/%s.chr6region.2.fastq" % (args.outdir,args.tname)

    cmd = "%s %s/tools/get_proper_aln_bam.py %s %s %s %s %s %s %s %s" % (args.py3,bin_dir,fq1_tumor,fq2_tumor,args.tname,hla_fa,args.bwa,args.samtools,args.sbb,args.outdir)
    of.write(cmd+'\n')

    # aln two allele and get het pos and cal BAF
    of.write('\n'+'###baf2loh main script for tumor'+'\n')
    tumor_proper_aln_bam = "%s/%s.chr6region.patient.reference.hlas.sort.markdup.proper_aln.bam" % (args.outdir,args.tname)
    cmd = "%s %s/BAF2LOH.py -indir %s -name %s -bam %s -jre %s" % (args.py3,bin_dir,args.outdir,args.tname,tumor_proper_aln_bam,args.jre)
    of.write(cmd+'\n\n\n\n')

    
    # copy number analysis

    # for normal
    cmd = "%s %s/tools/cnv.py %s %s %s %s %s %s %s %s %s %s %s" % (
        args.py3,
        bin_dir,
        args.nbam,
        args.nname,
        args.bed,
        args.fasta,
        args.cnvRefDir,
        args.sbb,
        args.bedtools,
        args.rscript,
        bin_dir,
        args.py3,
        args.outdir
        )
    of.write(cmd+'\n')

    # for tumor
    cmd = "%s %s/tools/cnv.py %s %s %s %s %s %s %s %s %s %s %s" % (
        args.py3,
        bin_dir,
        args.tbam,
        args.tname,
        args.bed,
        args.fasta,
        args.cnvRefDir,
        args.sbb,
        args.bedtools,
        args.rscript,
        bin_dir,
        args.py3,
        args.outdir
        )
    of.write(cmd+'\n')



    # BAF to LOH main method
    cmd = "perl %s/tools/determine_HLALOH_by_BAF.pl %s %s %s %s" % (bin_dir,args.outdir,args.nname,args.tname,args.outdir)
    of.write(cmd+'\n')

    of.close()

    
if __name__ == "__main__":
    main()
