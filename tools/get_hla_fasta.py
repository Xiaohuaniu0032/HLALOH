import os
import sys
import argparse
import re
from collections import defaultdict

def parse_args():
    AP = argparse.ArgumentParser("get hla A/B/C allels' fasta sequence")
    AP.add_argument('-i',help='OptiType result file',dest='hla_res')
    AP.add_argument('-o',help='outfile',dest='outfile')
    AP.add_argument('-o2',help='output all alleles belong to this 4-digit',dest='outfile2')

    return AP.parse_args()

def main():
    args = parse_args()

    hla2fa = {}
    allele = {}

    this_dir = os.path.split(os.path.realpath(__file__))[0]
    upper_dir = os.path.dirname(this_dir)
    #hla_fa = "%s/OptiType-1.3.2/data/hla_reference_dna.fasta" % (upper_dir)
    #hla_fa = "%s/DB/hla_fa_hg19/hla_abc_gen.fasta" % (upper_dir)
    hla_fa = "%s/DB/hla_fa_hg19/hla_abc_cds.fasta" % (upper_dir)
    print("IMGT-HLA database is: %s" % (hla_fa))

    hla = open(hla_fa,'r')
    for line in hla:
        if '>HLA' in line:
            hla_allele = line.strip().split(' ')[1]
            hla2fa[hla_allele] = [] # store each line
        else:
            hla2fa[hla_allele].append(line.strip())
    hla.close()


    # cat into one seq
    hla2fa_new = {}
    for x in hla2fa:
        seqs = hla2fa[x] # a list
        new_seq = ""
        for s in seqs:
            new_seq  = new_seq + s

        hla2fa_new[x] = new_seq

    #print(hla2fa_new)


    # get each A/B/C alleles' fasta
    infile = open(args.hla_res,'r')
    infile.readline() # skip header
    val = infile.readline().split('\t')
    infile.close()

    hla_a1 = val[1] # A*29:01
    hla_a2 = val[2] # A*30:01

    hla_b1 = val[3]
    hla_b2 = val[4]

    hla_c1 = val[5]
    hla_c2 = val[6]


    outfile = open(args.outfile,'w')
    outfile2 = open(args.outfile2,'w')

    final_allele = defaultdict(list)

    for x in hla2fa_new:
        # 循环每个hla allele
        seq = hla2fa_new[x]
        a = ">%s" % (x)

        if hla_a1 in x:
            final_allele[hla_a1].append(seq) # hla_a1是4位,x可能为4-8位

            # write for outfile2
            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')

        if hla_a2 in x:
            final_allele[hla_a2].append(seq)

            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')

        if hla_b1 in x:
            final_allele[hla_b1].append(seq)

            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')

        if hla_b2 in x:
            final_allele[hla_b2].append(seq)

            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')

        if hla_c1 in x:
            final_allele[hla_c1].append(seq)

            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')

        if hla_c2 in x:
            final_allele[hla_c2].append(seq)

            outfile2.write(a+'\n')
            outfile2.write(seq+'\n')
            
        
    # write file
    hla_a1_12 = val[1].split('*')[1].split(':')[0] # 1-2位
    hla_a1_34 = val[1].split('*')[1].split(':')[1] # 3-4位
    a = ">hla_a_%s_%s" % (hla_a1_12,hla_a1_34)
    seq = final_allele[hla_a1][0] # 取第一个allele（可能为4-8位）作为该4位allele的ref
    outfile.write(a+'\n')
    outfile.write(seq+'\n')




    hla_a2_12 = val[2].split('*')[1].split(':')[0]
    hla_a2_34 = val[2].split('*')[1].split(':')[1]
    a = ">hla_a_%s_%s" % (hla_a2_12,hla_a2_34)
    seq = final_allele[hla_a2][0]
    outfile.write(a+'\n')
    outfile.write(seq+'\n')


    hla_b1_12 = val[3].split('*')[1].split(':')[0]
    hla_b1_34 = val[3].split('*')[1].split(':')[1]
    a = ">hla_b_%s_%s" % (hla_b1_12,hla_b1_34)
    seq = final_allele[hla_b1][0]
    outfile.write(a+'\n')
    outfile.write(seq+'\n')


    hla_b2_12 = val[4].split('*')[1].split(':')[0]
    hla_b2_34 = val[4].split('*')[1].split(':')[1]
    a = ">hla_b_%s_%s" % (hla_b2_12,hla_b2_34)
    seq = final_allele[hla_b2][0]
    outfile.write(a+'\n')
    outfile.write(seq+'\n')


    hla_c1_12 = val[5].split('*')[1].split(':')[0]
    hla_c1_34 = val[5].split('*')[1].split(':')[1]
    a = ">hla_c_%s_%s" % (hla_c1_12,hla_c1_34)
    seq = final_allele[hla_c1][0]
    outfile.write(a+'\n')
    outfile.write(seq+'\n')


    hla_c2_12 = val[6].split('*')[1].split(':')[0]
    hla_c2_34 = val[6].split('*')[1].split(':')[1]
    a = ">hla_c_%s_%s" % (hla_c2_12,hla_c2_34)
    seq = final_allele[hla_c2][0]
    outfile.write(a+'\n')
    outfile.write(seq+'\n')


    outfile.close()



if __name__ == "__main__":
    main()






