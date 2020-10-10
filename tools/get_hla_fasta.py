import os
import sys
import argparse
import re

def parse_args():
    AP = argparse.ArgumentParser("get hla A/B/C allels' fasta sequence")
    AP.add_argument('-i',help='OptiType result file',dest='hla_res')
    AP.add_argument('-fa',help='HLA all alleles fasta database',dest='fasta')
    AP.add_argument('-o',help='outfile',dest='outfile')

    return AP.parse_args()

def main():
    args = parse_args()

    hla2fa = {}
    allele = {}
    hla_fa = open(args.fasta,'r')
    for line in hla_fa:
        #hla = line.strip().split(' ')[1].split('-')[1]
        #abc = hla.split('*')[0]
        #first = hla.split('*')[1].split(':')[0]
        #second = hla.split('*')[1].split(':')[1]


        if re.findall("^>HLA",line):
            hla_allele = line.strip().split(' ')[1]
            if hla_allele not in allele:
                hla2fa[hla_allele] = []
                allele[hla_allele] = 1
        else:
            hla2fa[hla_allele].append(line.strip())

    #print(hla2fa)
    hla_fa.close()



    # cat into one seq
    hla2fa_new = {}
    for x in hla2fa:
        seqs = hla2fa[x]
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

    hla_a1 = val[1]
    hla_a2 = val[2]

    hla_b1 = val[3]
    hla_b2 = val[4]

    hla_c1 = val[5]
    hla_c2 = val[6]


    outfile = open(args.outfile,'w')

    abc_allele_flag = {}

    for x in hla2fa_new:
        #print(x)
        seq = hla2fa_new[x]
        #print(hla_a1)
        if hla_a1 in x:
            print(x)
            # get new allele format
            #print(val[1])
            print(val[1].split('*')[1].split(':'))
            hla_a1_12 = val[1].split('*')[1].split(':')[0]
            #print(hla_a1_12)
            hla_a1_34 = val[1].split('*')[1].split(':')[1]

            #val = "hla_a_%s_%s" % (hla_a1_12,hla_a1_34)
            #print(hla_a1_12)
            #print(hla_a1_34)
            #val = str(hla_a1_12)
            #print(val)
            #val = "hla_a_" + str(hla_a1_12) + '_' + str(hla_a1_34)
            #outfile.write(val+'\n')
            #outfile.write(seq+'\n')
            
            #abc_allele_flag["hla_a1"] = 1 # mark


        '''
        if hla_a2 in x:
            hla_a2_12 = val[2].split('*')[1].split(':')[0]
            hla_a2_34 = val[2].split('*')[1].split(':')[1]

            val = ">hla_a_%s_%s" % (hla_a2_12,hla_a2_34)
            outfile.write(val+'\n')
            outfile.write(seq+'\n')
            
            abc_allele_flag["hla_a2"] = 1

        if hla_b1 in x:
            hla_b1_12 = val[3].split('*')[1].split(':')[0]
            hla_b1_34 = val[3].split('*')[1].split(':')[1]

            val = ">hla_b_%s_%s" % (hla_b1_12,hla_b1_34)
            outfile.write(val+'\n')
            outfile.write(seq+'\n')

            abc_allele_flag["hla_b1"] = 1

        if hla_b2 in x:
            hla_b2_12 = val[4].split('*')[1].split(':')[0]
            hla_b2_34 = val[4].split('*')[1].split(':')[1]

            val = "hla_b_%s_%s" % (hla_b2_12,hla_b2_34)
            outfile.write(val+'\n')
            outfile.write(seq+'\n')

            abc_allele_flag["hla_b2"] = 1

        if hla_c1 in x:
            hla_c1_12 = val[5].split('*')[1].split(':')[0]
            hla_c1_34 = val[5].split('*')[1].split(':')[1]

            val = "hla_c_%s_%s" % (hla_c1_12,hla_c1_34)
            outfile.write(val+'\n')
            outfile.write(seq+'\n')

            abc_allele_flag["hla_c1"] = 1

        if hla_c2 in x:
            hla_c2_12 = val[6].split('*')[1].split(':')[0]
            hla_c2_34 = val[6].split('*')[1].split(':')[1]

            val = "hla_c_%s_%s" % (hla_c2_12,hla_c2_34)
            outfile.write(val+'\n')
            outfile.write(seq+'\n')

            abc_allele_flag["hla_c2"] = 1
        '''
    outfile.close()



if __name__ == "__main__":
    main()






