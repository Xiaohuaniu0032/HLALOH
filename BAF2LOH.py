import sys
import os
import configparser
import argparse
import glob

def parse_args():
    AP = argparse.ArgumentParser("Convert BAF into LOH")
    AP.add_argument('-indir',help='result dir',dest='indir')
    AP.add_argument('-name',help='sample name',dest='name')
    AP.add_argument('-bam',help='bam',dest='bam')
    #AP.add_argument('-r',help='hla ref',dest='ref')
    
    return AP.parse_args()

def main():
    args = parse_args()
    hla_typing_file = "%s/hla.result.new" % (args.indir)
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    
    # get all alleles
    new_hla_res = "%s/hla.result.new" % (args.indir)
    hla_alleles = []
    hla_allele = open(new_hla_res,'r')
    for line in hla_allele:
        hla_alleles.append(line.strip())
    hla_allele.close()

    proper_align_bam = args.bam
    
    # for each allele, get its bam file
    for a in hla_alleles:
        temp_bam = "%s/%s.temp.%s.bam" % (args.indir,args.name,a)
        cmd = "samtools view -b -o %s %s %s" % (temp_bam,proper_align_bam,a)
        os.system(cmd)
        #of.write(cmd+'\n')

        # sort
        sort_bam = "%s/%s.type.%s.bam" % (args.indir,args.name,a)
        cmd = "samtools sort %s -o %s" % (temp_bam,sort_bam)
        os.system(cmd)
        #of.write(cmd+'\n')
        
        # index
        cmd = "samtools index %s" % (sort_bam)
        os.system(cmd)
        #of.write(cmd+'\n')
        
        # pileup file
        pileupFile = "%s/%s.%s.pileup" % (args.indir,args.name,a)
        hla_fa = "%s/patient.hlaFasta.fa" % (args.indir)
        cmd = "samtools mpileup -x -Q 0 -f %s %s >%s" % (sort_bam,pileupFile)
        os.system(cmd)
        #of.write(cmd+'\n')
        
        # remove temp files
        if os.path.exists(temp_bam):
            os.remove(temp_bam)

    # get het pos
    cmd = "Rscript %s/tools/hla_het_snp.R %s" % (bin_dir,args.indir)
    os.system(cmd)
    #of.write(cmd+'\n')



    # for homo alleles, you will not get hla_*_aln.txt file(s)

    # cal BAF
    hla_a1 = hla_alleles[0]
    hla_a2 = hla_alleles[1]

    hla_b1 = hla_alleles[2]
    hla_b2 = hla_alleles[3]

    hla_c1 = hla_alleles[4]
    hla_c2 = hla_alleles[5]

    hla_a1_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_a1)
    hla_a2_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_a2)

    hla_b1_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_b1)
    hla_b2_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_b2)

    hla_c1_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_c1)
    hla_c2_pileup = "%s/%s.%s.pileup" % (args.indir,args.name,hla_c2)

    hla_a_het_file = "%s/hla_a_aln.txt" % (args.indir)
    hla_b_het_file = "%s/hla_b_aln.txt" % (args.indir)
    hla_c_het_file = "%s/hla_c_aln.txt" % (args.indir)

    hla_a_baf = "%s/%s.hla_a_BAF.txt" % (args.indir,args.name)
    hla_b_baf = "%s/%s.hla_b_BAF.txt" % (args.indir,args.name)
    hla_c_baf = "%s/%s.hla_c_BAF.txt" % (args.indir,args.name)

    # for empty file's header
    h = 'Rank' + '\t' + 'pos1' + '\t' + 'pos2' + '\t' + 'base1' + '\t' + 'base2' + '\t' + 'num1' + '\t' + 'num2' + '\t' + 'depth' + '\t' + 'baf1' + '\t' + 'baf2'

    if os.path.exists(hla_a_het_file):
        cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_a_het_file,hla_a1_pileup,hla_a2_pileup,hla_a_baf)
        #print("CMD is: %s" % (cmd))
        os.system(cmd)
    else:
        # make an empty file
        of = open(hla_a_baf,'w')
        of.write(h+'\n')
        of.close()

    if os.path.exists(hla_b_het_file):
        cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_b_het_file,hla_b1_pileup,hla_b2_pileup,hla_b_baf)
        #print(cmd)
        os.system(cmd)
    else:
        of = open(hla_b_baf,'w')
        of.write(h+'\n')
        of.close()

    if os.path.exists(hla_c_het_file):
        cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_c_het_file,hla_c1_pileup,hla_c2_pileup,hla_c_baf)
        #print(cmd)
        os.system(cmd)
    else:
        of = open(hla_c_baf,'w')
        of.write(h+'\n')
        of.close()

    # plot baf fig
    cmd = "Rscript %s/tools/hla_baf.r %s %s %s %s %s" % (bin_dir,hla_a_baf,hla_b_baf,hla_c_baf,args.name,args.indir)
    os.system(cmd)
    #of.write(cmd+'\n')



if __name__ == "__main__":
    main()
