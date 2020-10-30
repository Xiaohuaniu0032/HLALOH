import sys
import os
import configparser
import argparse
import glob

def parse_args():
    AP = argparse.ArgumentParser("Convert BAF into LOH")
    AP.add_argument('-indir',help='result dir',dest='indir')
    AP.add_argument('tname',help='tumor name',dest='tname')
    
    return AP.parse_args()

def main():
    hla_typing_file = "%s/hla.result.new" % (args.indir)
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    #runsh = "%s/%s.BAF2LOH.sh" % (args.indir,args.tname)
    #of = open(runsh,'w')

    # get each allele's bam file
    hla_alleles = []
    hla_allele = open(new_hla_res,'r')
    for line in hla_allele:
        hla_alleles.append(strip(line))
    hla_allele.close()

    proper_align_bam = "%s/%s.chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam" % (args.indir,args.tname)

    for a in hla_alleles:
        # for each allele, get its bam file
        temp_bam = "%s/%s.temp.%s.bam" % (args.indir,args.tname,a)
        cmd = "samtools view -b -o %s %s %s" % (temp_bam,proper_align_bam,a)
        os.system(cmd)
        #of.write(cmd+'\n')

        # sort
        sort_bam = "%s/%s.type.%s.bam" % (args.indir,args.tname,a)
        cmd = "samtools sort %s -o %s" % (temp_bam,sort_bam)
        os.system(cmd)
        #of.write(cmd+'\n')
        
        # index
        cmd = "samtools index %s" % (sort_bam)
        os.system(cmd)
        #of.write(cmd+'\n')
        
        # pileup file
        pileupFile = "%s/%s.%s.pileup" % (args.indir,args.tname,a)
        cmd = "samtools mpileup -f %s %s >%s" % (hla_fa,sort_bam,pileupFile)
        os.system(cmd)
        #of.write(cmd+'\n')

    # get het pos
    cmd = "Rscript %s/tools/hla_het_snp.R %s" % (bin_dir,args.indir)
    os.system(cmd)
    #of.write(cmd+'\n')

    # cal BAF
    hla_a1 = hla_alleles[0]
    hla_a2 = hla_alleles[1]

    hla_b1 = hla_alleles[2]
    hla_b2 = hla_alleles[3]

    hla_c1 = hla_alleles[4]
    hla_c2 = hla_alleles[5]

    hla_a1_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_a1)
    hla_a2_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_a2)

    hla_b1_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_b1)
    hla_b2_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_b2)

    hla_c1_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_c1)
    hla_c2_pileup = "%s/%s.%s.pileup" % (args.indir,args.tname,hla_c2)

    hla_a_het_file = "%s/hla_a_aln.txt" % (args.indir)
    hla_b_het_file = "%s/hla_b_aln.txt" % (args.indir)
    hla_c_het_file = "%s/hla_c_aln.txt" % (args.indir)

    hla_a_baf = "%s/hla_a_BAF.txt" % (args.indir)
    hla_b_baf = "%s/hla_b_BAF.txt" % (args.indir)
    hla_c_baf = "%s/hla_c_BAF.txt" % (args.indir)

    cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_a_het_file,hla_a1_pileup,hla_a2_pileup,hla_a_baf)
    os.system(cmd)
    #of.write(cmd+'\n')

    cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_b_het_file,hla_b1_pileup,hla_b2_pileup,hla_b_baf)
    os.system(cmd)
    #of.write(cmd+'\n')

    cmd = "perl %s/tools/cal_het_pos_baf_using_align_info.pl -het %s -a1 %s -a2 %s -of %s" % (bin_dir,hla_c_het_file,hla_c1_pileup,hla_c2_pileup,hla_c_baf)
    os.system(cmd)
    #of.write(cmd+'\n')

    # plot baf fig
    cmd = "Rscript %s/tools/hla_baf.r %s %s %s %s %s" % (bin_dir,hla_a_baf,hla_b_baf,hla_c_baf,args.tname,args.outdir)
    os.system(cmd)
    #of.write(cmd+'\n')



if __name__ == "__main__":
    main()
