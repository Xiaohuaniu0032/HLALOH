import sys
import configparser
import argparse
import os
#import pandas as pd

def parse_args():
    AP = argparse.ArgumentParser("detect somatic copy number from capture NGS data")
    AP.add_argument('-bam',help='bam file',dest='bam')
    AP.add_argument('-bed',help='bed file',dest='bed')
    AP.add_argument('-n',help='sample name',dest='name')
    AP.add_argument('-fa',help='fasta file',dest='fasta',default='/data1/database/b37/human_g1k_v37.fasta')
    AP.add_argument('-m',help='analysis mode. can be [ref|cnv]',default='ref',dest='mode')
    AP.add_argument('-ref',help='control dir',dest='ref')
    AP.add_argument('-od',help='out dir',dest='outdir')

    return AP.parse_args()
    



def main():
    args = parse_args()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    config_file = bin_dir + '/config.ini'
    #print(config_file)
    config = configparser.ConfigParser()
    config.read(config_file)

    # software
    perl = config['software']['perl']
    py3  = config['software']['python3']
    bedtools = config['software']['bedtools']
    sambamba = config['software']['sambamba']
    rscript  = config['software']['Rscript']
    hg19 = config['db']['hg19']

    bin_dir = os.path.split(os.path.realpath(__file__))[0]

    runsh = args.outdir + '/ref.%s.sh' % (args.name)
    print("analysis mode is: ref...")
    f = open(runsh,'w')

    # cal depth
    cmd = "%s %s/bin/cal_depth.pl -bam %s -n %s -bed %s -sbb %s -outdir %s" % (perl,bin_dir,args.bam,args.name,args.bed,sambamba,args.outdir)
    f.write(cmd+'\n')

    # add gc
    depth_file = args.outdir + "/%s.targetcoverage.cnn" % (args.name)
    cmd = "%s %s/bin/add_gc.pl -bed %s -depth %s -bedtools_bin %s -fa %s -outdir %s" % (perl,bin_dir,args.bed,depth_file,bedtools,args.fasta,args.outdir)
    f.write(cmd+'\n')

    # infer sex
    sex_of = "%s/sex.txt" % (args.outdir)
    cmd = "%s %s/bin/infer_sex.py -cov %s -o %s" % (py3,bin_dir,depth_file,sex_of)
    f.write(cmd+'\n')

    # gc correct
    cmd = "%s %s/bin/gc_correct.r %s/%s.targetcoverage.cnn.with.gc.xls %s %s/%s.targetcoverage.cnn.with.gc.xls.gc.corrected.xls" % (rscript,bin_dir,args.outdir,args.name,sex_of,args.outdir,args.name)
    f.write(cmd+'\n')

    # lib normalize
    cmd = "%s %s/bin/normalize.pl -d %s/%s.targetcoverage.cnn.with.gc.xls.gc.corrected.xls -od %s" % (perl,bin_dir,args.outdir,args.name,args.outdir)
    f.write(cmd+'\n')

    if args.mode == 'ref':
        pass
    else:
        # cal logR
        cmd = "%s %s/bin/cal_logR.pl"
        f.write(cmd)

    f.close()

if __name__ == "__main__":
    main()

