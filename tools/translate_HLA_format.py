import os
import sys
import argparse

def parse_args():
    AP = argparse.ArgumentParser("translate OptiType HLA result into lohhla input format")
    AP.add_argument('-i',help='infile',dest='infile')
    AP.add_argument('-o',help='outfile',dest='outfile')
    return AP.parse_args()


def main():
    args = parse_args()

    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')

    infile.readline() # skip header
    val = infile.readline().split('\t')
    infile.close()

    #print(val[1])

    
    hla_a1_12 = val[1].split('*')[1].split(':')[0]
    hla_a1_34 = val[1].split('*')[1].split(':')[1]

    hla_a2_12 = val[2].split('*')[1].split(':')[0]
    hla_a2_34 = val[2].split('*')[1].split(':')[1]

    
    hla_b1_12 = val[3].split('*')[1].split(':')[0]
    hla_b1_34 = val[3].split('*')[1].split(':')[1]

    hla_b2_12 = val[4].split('*')[1].split(':')[0]
    hla_b2_34 = val[4].split('*')[1].split(':')[1]

    hla_c1_12 = val[5].split('*')[1].split(':')[0]
    hla_c1_34 = val[5].split('*')[1].split(':')[1]

    hla_c2_12 = val[6].split('*')[1].split(':')[0]
    hla_c2_34 = val[6].split('*')[1].split(':')[1]

    val = "hla_a_%s_%s" % (hla_a1_12,hla_a1_34)
    outfile.write(val+'\n')

    val = "hla_a_%s_%s" % (hla_a2_12,hla_a2_34)
    outfile.write(val+'\n')

    val = "hla_b_%s_%s" % (hla_b1_12,hla_b1_34)
    outfile.write(val+'\n')

    val = "hla_b_%s_%s" % (hla_b2_12,hla_b2_34)
    outfile.write(val+'\n')

    val = "hla_c_%s_%s" % (hla_c1_12,hla_c1_34)
    outfile.write(val+'\n')

    val = "hla_c_%s_%s" % (hla_c2_12,hla_c2_34)
    outfile.write(val+'\n')
    

    outfile.close()


if __name__ == "__main__":
    main()


