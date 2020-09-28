import sys
import os
import configparser


def parse_args():
    AP = argparse.ArgumentParser("detect HLA LOH from capture NGS data using paired tumor/normal")
    AP.add_argument('-nbam',help='bam file',dest='bam')
    AP.add_argument('-bed',help='bed file',dest='bed')
    AP.add_argument('-n',help='sample name',dest='name')
    AP.add_argument('-fa',help='fasta file',dest='fasta',default='/data1/database/b37/human_g1k_v37.fasta')
    AP.add_argument('-ref',help='cnv ref control dir',dest='ref')
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

    java = config['software']['java']
    gatk_dir = config['software']['gatk_dir']


    runsh = args.outdir + '/%s.HLALOH.sh' % (args.name)
    of = open(runsh,'w')

    # extract HLA region read1/read2
    extract_HLA_reads(args.bam,
                        args.name,
                        args.samtools,
                        java,
                        gatk_dir,
                        outdir,
                        of
                        )

    # HLA typing using OptiType
    fq1 = "%s/%s.chr6region.1.fastq"
    fq2 = "%s/%s.chr6region.2.fastq"

    cmd = "%s %s/OptiType-1.3.2/OptiTypePipeline.py -i %s %s --dna -v -o %s -c %s/OptiType-1.3.2/config.ini.example" % (py3,
        bin_dir,
        fq1,
        fq2,
        args.outdir,
        bin_dir
        )

    of.write(cmd+'\n')

    # transformat HLA result
    hla_result_file = "%s/%s." % ()


    # get HLA two alleles's fasta sequence

    # BAF
    cmd = 

    # logR
    cmd = 

    # estimate tumor purity/ploidy using ASCAT
    cmd = 

    # HLA main script
    cmd = '%s %s/lohhla/LOHHLAscript.R --patientId %s --outputDir %s --normalBAMfile %s --BAMDir %s --hlaPath %s --HLAfastaLoc %s --CopyNumLoc %s --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir %s --novoDir %s' % ()
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



    return $chr_naming


def extract_HLA_reads(bam,name,samtools,java,gatk_dir,outdir,runshFH):
    cmd = "%s view -H %s >%s/%s.hla.sam" % (samtools,
                                            bam,
                                            outdir,
                                            name
                                            )

    runshFH.write(cmd+'\n')

    #os.system(cmd)

    chr_name = chr_naming(bam,
                        samtools,
                        outdir)

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

    cmd = "%s -jar %s/SamToFastq.jar I=%s F=%s F2=%s VALIDATION_STRINGENCY=SILENT" % (java,
                                                                                        gatk_dir,
                                                                                        fq1,
                                                                                        fq2
                                                                                        )

    runshFH.write(cmd+'\n')
    #os.system(cmd)


def translate_HLA_format(raw_res,new_res):
    

def get_HLA_fasta(hla_a1,hla_a2,hla_b1,hla_b2,hla_c1,hla_c2,hla_fasta,outfile):
    




