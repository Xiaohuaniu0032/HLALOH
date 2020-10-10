import os
import sys
import argparse
import pysam
from collections import defaultdict

def parse_args():
    AP = argparse.ArgumentParser("Convert SAM/BAM file into fastq files")
    AP.add_argument('-f',help='sam/bam file',dest='infile')
    AP.add_argument('-t',help='file type <sam or bam>',dest='type',default='sam')
    AP.add_argument('-n',help='sample name',dest='name')
#    AP.add_argument('-fq1',help='fq1 file',dest='fastq1')
#    AP.add_argument('-fq2',help='fq2 file',dest='fastq2')
    AP.add_argument('-od',help='outdir',dest='outdir')

    return AP.parse_args()

def main():
    args = parse_args()
    fq1 = "%s/%s.R1.fastq" % (args.outdir, args.name)
    fq2 = "%s/%s.R2.fastq" % (args.outdir, args.name)

    fq1_fh = open(fq1,'w')
    fq2_fh = open(fq2,'w')

    proper_paired = get_paired_names(args.infile)
    #print(proper_paired)

    samfile = pysam.AlignmentFile(file,'r')
    



def file_type_check(file):
    '''
    auto detect input file type <SAM or BAM> by its suffix (*.[bam|sam])
    '''

    next

def get_paired_names(file):
    samfile = pysam.AlignmentFile(file,'r')
    readNameDic = defaultdict(int)
    
    for read in samfile.fetch():
        if read.is_qcfail or read.is_secondary or read.is_supplementary:
            continue
        read_name = read.query_name
        #print(read_name)
        readNameDic[read_name] += 1

    proper_paired = []
    for x,y in readNameDic.items():
        if y == 2:
            proper_paired.append(x)

    samfile.close()
    #print(proper_paired)

    # re-check if the two reads are r1 and r2
    samfile = pysam.AlignmentFile(file,'r') # re-open
    final_name = []
    r1r2 = defaultdict(list)

    for read in samfile.fetch():
        if read.is_qcfail or read.is_secondary or read.is_supplementary:
            continue

        read_name = read.query_name
        #print(read_name)

        if read_name in proper_paired:
            #print(read_name)
            #print(read_name)
            # this read has two records
            if read.is_read1:
                #print("this read is R1")
                r1r2[read_name].append('r1')
            if read.is_read2:
                #print("this read is R2")
                r1r2[read_name].append('r2')


    #print(r1r2)

    for x,y in r1r2.items():
        if 'r1' in y and 'r2' in y:
            final_name.append(x)

    #print(final_name)
    samfile.close()

    return(final_name)







if __name__ == "__main__":
    main()
