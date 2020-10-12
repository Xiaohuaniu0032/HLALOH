import sys
import os
import configparser
import argparse
import glob

def parse_args():
    AP = argparse.ArgumentParser("detect HLA LOH from capture NGS data using paired tumor/normal")
    AP.add_argument('-nbam',help='normal bam file',dest='nbam')
    AP.add_argument('-tbam',help='tumor bam file',dest='tbam')
    AP.add_argument('-bamDir',help='dir contains tumor and normal bam files',dest='bamdir')
    AP.add_argument('-bed',help='bed file',dest='bed')
    AP.add_argument('-nname',help='normal sample name',dest='nname')
    AP.add_argument('-tname',help='tumor sample name',dest='tname')
    AP.add_argument('-fa',help='fasta file',dest='fasta',default='/data1/database/b37/human_g1k_v37.fasta')
    AP.add_argument('-ref',help='cnv ref control dir',dest='ref')
    AP.add_argument('-gatkDir',help='gatk dir',dest='gatkDir')
    AP.add_argument('-novoDir',help='novoalign bin dir',dest='novoDir')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()


def main():
    args = parse_args()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]

    config_file = bin_dir + '/config.ini'
    config = configparser.ConfigParser()
    config.read(config_file)

    # software
    perl = config['software']['perl']
    py3  = config['software']['python3']
    py2 = config['software']['python2']
    samtools = config['software']['samtools']
    java = config['software']['java']
    gatk_dir = config['software']['gatk_dir']
    sambamba = config['software']['sambamba']
    bedtools = config['software']['bedtools']
    rscript = config['software']['rscript']


    runsh = "%s/%s.HLALOH.sh" % (args.outdir,args.tname)
    of = open(runsh,'w')

    # extract HLA region read1/read2
    extract_HLA_reads(args.nbam,
                        args.nname,
                        samtools,
                        bedtools,
                        java,
                        gatk_dir,
                        args.outdir,
                        of
                        )

    # HLA typing using OptiType
    fq1 = "%s/%s.chr6region.1.fastq" % (args.outdir,args.nname)
    fq2 = "%s/%s.chr6region.2.fastq" % (args.outdir,args.nname)

    cmd = "%s %s/OptiType-1.3.2/OptiTypePipeline.py -i %s %s --dna -v -o %s -c %s/OptiType-1.3.2/config.ini.example" % (py2,
        bin_dir,
        fq1,
        fq2,
        args.outdir,
        bin_dir
        )

    of.write(cmd+'\n')


    # transformat HLA result
    cmd = "cp %s/*/*_result.tsv %s/hla.result.raw" % (args.outdir,args.outdir)
    of.write(cmd+'\n')

    #hla_res = glob.glob("%s/*/*_result.tsv" % args.outdir)[0]
    #print(hla_res)
    raw_hla_res = "%s/hla.result.raw" % (args.outdir)
    new_hla_res = "%s/hla.result.new" % (args.outdir)
    cmd = "%s %s/tools/translate_HLA_format.py -i %s -o %s" % (py3,bin_dir,raw_hla_res,new_hla_res)

    of.write(cmd+'\n')


    # get HLA alleles' fasta sequence



    # BAF
    nvaf = "%s/%s.normal.vaf" % (args.outdir,args.nname)
    tvaf = "%s/%s.tumor.vaf" % (args.outdir,args.tname)

    cmd = "%s %s/BAF/pileup2vaf.v2.py -bam %s -bed %s -outfile %s" % (py3,
                                                                        bin_dir,
                                                                        args.nbam,
                                                                        args.bed,
                                                                        nvaf
                                                                        )
    of.write(cmd+'\n')

    cmd = "%s %s/BAF/pileup2vaf.v2.py -bam %s -bed %s -outfile %s" % (py3,
                                                                        bin_dir,
                                                                        args.tbam,
                                                                        args.bed,
                                                                        tvaf
                                                                        )
    of.write(cmd+'\n')


    # logR
    bams = [args.nbam,args.tbam]
    names = [args.nname,args.tname]
    idx = 0
    for i in bams:
        name = names[idx]
        # s1. calculate depth
        cmd = "%s %s/CNVscan/bin/cal_depth.pl -bam %s -n %s -bed %s -sbb %s -outdir %s" % (perl,
                                                                                        bin_dir,
                                                                                        i,
                                                                                        name,
                                                                                        args.bed,
                                                                                        sambamba,
                                                                                        args.outdir
                                                                                        )
        of.write(cmd+'\n')

        
        # s2. add gc
        depth_file = "%s/%s.targetcoverage.cnn" % (args.outdir,name)
        cmd = "%s %s/CNVscan/bin/add_gc.pl -bed %s -depth %s -bedtools_bin %s -fa %s -outdir %s" % (perl,
                                                                                                bin_dir,
                                                                                                args.bed,
                                                                                                depth_file,
                                                                                                bedtools,
                                                                                                args.fasta,
                                                                                                args.outdir
                                                                                                )

        of.write(cmd+'\n')

        # s3. infer sex
        cmd = "%s %s/CNVscan/bin/infer_sex.py -cov %s -o %s/sex.txt" % (py3,bin_dir,depth_file,args.outdir)
        of.write(cmd+'\n')

        # s4. gc correct
        add_gc_depth = "%s/%s.targetcoverage.cnn.with.gc.xls" % (args.outdir,name)
        sex_file = "%s/sex.txt" % (args.outdir)
        gc_correct_depth = "%s/%s.targetcoverage.cnn.with.gc.xls.gc.corrected.xls" % (args.outdir,name)
        cmd = "%s %s/CNVscan/bin/gc_correct.r %s %s %s" % (rscript,bin_dir,add_gc_depth,sex_file,gc_correct_depth)
        of.write(cmd+'\n')

        # s5. normalize
        cmd = "%s %s/CNVscan/bin/normalize.pl -d %s -od %s" % (perl,bin_dir,gc_correct_depth,args.outdir)
        of.write(cmd+'\n')

        # s6. make ref matrix
        ref_mat = "%s/ref.matrix.txt" % (args.outdir)
        cmd = "%s %s/CNVscan/bin/make_ref_matrix.pl %s %s" % (perl,bin_dir,args.ref,ref_mat)
        of.write(cmd+'\n')

        # s7. cal logR
        norm_file = "%s/%s.norm.xls" % (args.outdir,name)
        logR = "%s/%s.logR.xls" % (args.outdir,name)
        cmd = "%s %s/CNVscan/bin/cal_logR.pl %s %s %s" % (perl,bin_dir,norm_file,ref_mat,logR)
        of.write(cmd+'\n')

        idx += 1




    # estimate tumor purity/ploidy using ASCAT
    # get ascat.logR file
    tumor_ascatlogR = "%s/%s.ascat.logR.xls" % (args.outdir,args.tname)
    normal_ascatlogR = "%s/%s.ascat.logR.xls" % (args.outdir,args.nname)

    tumor_logR = "%s/%s.logR.xls" % (args.outdir,args.tname)
    normal_logR = "%s/%s.logR.xls" % (args.outdir,args.nname)

    cmd = "%s %s/ASCAT/ascatLogR.pl %s %s %s" % (perl,bin_dir,tvaf,tumor_logR,tumor_ascatlogR)
    of.write(cmd+'\n')

    cmd = "%s %s/ASCAT/ascatLogR.pl %s %s %s" % (perl,bin_dir,nvaf,normal_logR,normal_ascatlogR)
    of.write(cmd+'\n')

    # ASCAT
    ascatPurityPloidyFile = "%s/PurityPloidyEst.txt" % (args.outdir)
    cmd = "%s %s/ASCAT/ascat.r %s %s %s %s %s" % (rscript,bin_dir,tumor_ascatlogR,tvaf,normal_ascatlogR,nvaf,ascatPurityPloidyFile)
    of.write(cmd+'\n')

    # HLA LOH main script
    hla_alleles = "%s/%s.hla_alleles" % (args.outdir,args.nname)
    hla_fa = "%s/patient.hla.fa" % (args.outdir)

    cmd = '%s %s/lohhla/LOHHLAscript.R --patientId %s --outputDir %s --normalBAMfile %s --BAMDir %s --hlaPath %s --HLAfastaLoc %s --CopyNumLoc %s --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir %s --novoDir %s' % (rscript,
                                        bin_dir,
                                        args.tname,
                                        args.outdir,
                                        args.nbam,
                                        args.bamdir,
                                        hla_alleles,
                                        hla_fa,
                                        ascatPurityPloidyFile,
                                        args.gatkDir,
                                        args.novoDir
                                        )
    of.write(cmd+'\n')
    of.close()


def chr_naming(bam,samtools,outdir):
    cmd = "%s view %s | head -n 1 >%s/first.line.sam" % (samtools,bam,outdir)
    os.system(cmd)
    samfile = "%s/first.line.sam" % (outdir)
    sam = open(samfile,'r')
    line = sam.readline()
    chrom = line.split('\t')[2]

    if 'chr' in chrom:
        chr_naming = 'with_prefix'
    else:
        chr_naming = 'no_prefix'


    os.remove(samfile)

    return chr_naming


def extract_HLA_reads(bam,name,samtools,bedtools,java,gatk_dir,outdir,runshFH):
    cmd = "%s view -H %s >%s/%s.hla.sam" % (samtools,
                                            bam,
                                            outdir,
                                            name
                                            )

    runshFH.write(cmd+'\n')

    #os.system(cmd)

    #chr_name = chr_naming(bam,samtools,outdir)

    chr_name = 'no_prefix'

    if chr_name == 'no_prefix':
        region1 = '6:29909037-29913661'
        region2 = '6:31321649-31324964'
        region3 = '6:31236526-31239869'
    else:
        region1 = 'chr6:29909037-29913661'
        region2 = 'chr6:31321649-31324964'
        region3 = 'chr6:31236526-31239869'

    cmd = "%s view %s %s >>%s/%s.hla.sam" % (samtools,
                                            bam,
                                            region1,
                                            outdir,
                                            name
                                            )

    
    runshFH.write(cmd+'\n')
    #os.system(cmd)

    cmd = "%s view %s %s >>%s/%s.hla.sam" % (samtools,
                                                bam,
                                                region2,
                                                outdir,
                                                name
                                                )

    runshFH.write(cmd+'\n')
    #os.system(cmd)

    cmd = "%s view %s %s >>%s/%s.hla.sam" % (samtools,
                                                bam,
                                                region3,
                                                outdir,
                                                name
                                                )

    runshFH.write(cmd+'\n')
    #os.system(cmd)

    '''
    ignore chr6 contig

    1.chr6_apd_hap1
    2.chr6_cox_hap2
    3.chr6_dbb_hap3
    4.chr6_mann_hap4
    5.chr6_mcf_hap5
    6.chr6_qbl_hap6
    7.chr6_ssto_hap7
    '''
    
    # turn into fastq
    fq1 = "%s/%s.chr6region.1.fastq" % (outdir,name)
    fq2 = "%s/%s.chr6region.2.fastq" % (outdir,name)

    hla_sam = "%s/%s.hla.sam" % (outdir,name)

    # convert to bam
    hla_bam = "%s/%s.hla.bam" % (outdir,name)
    cmd = "%s view -b -o %s %s" % (samtools,hla_bam,hla_sam)
    runshFH.write(cmd+'\n')

    # sort by name
    bam_sort_by_name = "%s/%s.sort_by_name.bam" % (outdir,name)
    cmd = "%s sort -n %s -o %s" % (samtools,hla_bam,bam_sort_by_name)
    runshFH.write(cmd+'\n')

    # remove temp sam/bam

    # bedtools bamtofastq
    cmd = "%s bamtofastq -i %s -fq %s -fq2 %s" % (bedtools,bam_sort_by_name,fq1,fq2)
    runshFH.write(cmd+'\n')
    

    # Illegal Mate State
    # https://www.biostars.org/p/59521/

    #cmd = "%s -jar %s/SamToFastq.jar I=%s F=%s F2=%s VALIDATION_STRINGENCY=SILENT" % (java, gatk_dir, hla_sam, fq1, fq2)


    
if __name__ == "__main__":
    main()




