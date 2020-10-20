import sys
import os

infile, outfile = sys.argv[1:]

class pileupConfig(object):
    def __init__(self):
        self.base_set = ['A','T','C','G','N','a','t','c','g','n']
        self.alt_plus_ref_base_set = ['A','T','C','G','N']
        self.alt_rev_ref_base_set = ['a','t','c','g','n']
        self.base_start = ['.',',','^']


class vcf_format(object):
    def __init__(self):
        self.fileformat = 'VCFv4.3'
        self.source = os.path.abs(sys.argv[0])
        self.reference = 'NA'
        self.phasing = 'NA'

        self.INFO = []
        self.FILTER = []
        self.FORMAT = []

    def INFO(self,ID,Num,Type,Des):
        val = "##<ID=%s,Number=%s,Type=%s,Description=\"%s\"> % (ID,Num,Type,Des))"
        self.INFO.append(val)

    def FILTER(self,ID,Des):
        val = "##<ID=%s,Description=\"%s\"> % (ID,Des)"
        self.FILTER.append(val)

    def FORMAT(self,ID,Num,Type,Des):
        val = "##<ID=%s,Number=%s,Type=%s,Description=\"%s\"> % (ID,Num,Type,Des))"
        self.FORMAT.append(val)


def indel_len(pos,seq,lenSeq):# pos is 1-based
    indel_chr = seq[pos-1] # +/-
    indel_base_num = 0
    nnn = [0,1,2,3,4,5,6,7,8,9]
    left_string = seq[pos:]
    flag = 0
    NUM = []
    while (flag <= len(left_string)-1):
        flag += 1
        val = left_string[flag-1]
        if val in nnn:
            NUM.append(string(val))
            continue
        else:
            break

    indel = indel_chr + "".join(NUM)

    return(len(NUM)) # +3ATC, return 3

def open_file(file):
    of = open(file,'w')
    return(of)

def write_vcf_meta_Info(outfile):
    vcf_format()
    vcf_format.INFO("DP",1,"Integer","Depth of this sample")
    vcf_format.INFO("AF","A","Float","Alt allele frequency")
    vcf_format.FILTER("q10","Base quality <= 10")
    vcf_format.FILTER("lowDP","Depth <= 30")
    vcf_format.FORMAT("GT",1,"String","Genotype")
    vcf_format.FORMAT("GQ",1,"Integer","Genotype Quality")
    vcf_format.FORMAT("DP",1,"Integer","Read Depth")
    vcf_format.FORMAT("HQ",2,"Integer","Haplotype Quality")

    outfile.write(vcf_format.fileformat+'\n')
    outfile.write(vcf_format.source+'\n')
    outfile.write(vcf_format.reference+'\n')
    outfile.write(vcf_format.phasing+'\n')

    for i in vcf_format.INFO():
        outfile.write(i+'\n')

    for i in vcf_format.FILTER():
        outfile.write(i+'\n')

    for i in vcf_format.FORMAT():
        outfile.write(i+'\n')

of = open_file(outfile)
#write_vcf_meta_Info(of)

#header = ['chr','pos','ref_base','depth_samtools','depth_myself','ref_num','alt_num','vaf']
header = ['chr','pos','ref_base','ref_n','alt_n','all_n','alt_freq']
of.write('\t'.join(header)+'\n')

with open(infile,'r') as f:
    for line in f:
        arr = line.strip().split('\t')
        #print(arr)
        depth = int(arr[3])
        if depth == 0:
            continue
        seq = arr[4]
        baseQ = arr[5]

        plus_ref, rev_ref, plus_alt, rev_alt = (0,0,0,0)
        ins_num, del_num = 0,0

        flag = 0
        while (flag <= len(seq)-1):
            # for each base
            flag += 1
            c = seq[flag-1]
            # print(c)
            if c == '.':
                plus_ref += 1
            if c == ',':
                rev_ref += 1
            if c == '+':
                ins_num += 1
                n = indel_len(flag,seq,len(seq)) # skip len
                skip_len = n + len(str(n))
                flag += skip_len
                continue
            if c == '-':
                del_num += 1
                n = indel_len(flag,seq,len(seq))
                skip_len = n + len(str(n))
                flag += skip_len
                continue
            if c == '^':
                flag += 1
                continue
            if c in pileupConfig().alt_plus_ref_base_set:
                plus_alt += 1
            if c in pileupConfig().alt_rev_ref_base_set:
                rev_alt += 1

        ref_num = plus_ref + rev_ref
        alt_num = plus_alt + rev_alt
        depth = ref_num + alt_num

        if depth == 0:
            vaf = 0 # set as 0
        else:
            vaf = round(float(alt_num)/depth,3)

        val = arr[0] + '\t' + arr[1] + '\t' + arr[2] + '\t' + str(ref_num) + '\t' + str(alt_num) + '\t' + str(depth) + '\t' + str(vaf)
        of.write(val+'\n')
of.close()
