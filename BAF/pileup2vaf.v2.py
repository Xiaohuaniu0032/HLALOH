import pysam
import sys
import argparse

class Options:
    def __init__(self):
        self.parser =  argparse.ArgumentParser("Calculate VAF in Samtools mpileup Way")
        self.parser.add_argument('-bam',help="bam file",dest="bam")
        self.parser.add_argument('-bed',help="bed file",dest="bed")
        self.parser.add_argument('-fa',help="fasta file",dest="fasta",default="/data1/database/b37/human_g1k_v37.fasta")
        self.parser.add_argument('-outfile',help="output file",dest="outfile")

        args = self.parser.parse_args()
        
        if args.bam:
            self.bam = args.bam
        else:
            self.parser.error("BAM file not supplied")

        if args.bed:
            self.bed = args.bed
        else:
            self.parser.error("BED file not supplied")

        if args.fasta:
            self.fasta = args.fasta
        else:
            self.parser.error("Fasta file not supplied")

        if args.outfile:
            self.outfile = args.outfile
        else:
            self.parser.error("Output file not supplied")


def main():
    options = Options()
    pysam_sam = pysam.AlignmentFile(options.bam,"rb")
    pysam_fa = pysam.FastaFile(options.fasta)

    # get target
    target = []
    with open(options.bed,"r") as b_fh:
        for line in b_fh:
            arr = line.rstrip().split('\t')
            t = arr[0] + ':' + str(int(arr[1]) + 1) + '-' + arr[2] # chr:start-end, 1-based
            target.append(t)

    print(target)

    # output file
    of = open(options.outfile,"w")
    h = ['chr','pos','ref_base','ref_n','alt_n','all_n','alt_freq']
    hh = '\t'.join(h)
    of.write(hh)
    of.write('\n')

    for t in target:
        for pileupcol in pysam_sam.pileup(region=t,reference=options.fasta,truncate=True,stepper="samtools",fastafile=pysam_fa,max_depth=5000,min_base_quality=10,min_mapping_quality=10):
            # max depth = 5000X, mapQ = 10, baseQ = 10
            querybase = pileupcol.get_query_sequences(mark_matches=False,mark_ends=False,add_indels=False)
            col = pileupcol.reference_pos + 1
            t_chr = t.split(":")[0]
            reg = t_chr + ':' + str(col) + '-' + str(col) # chr1:2-2
            ref_base = pysam_fa.fetch(region=reg).upper() # ref base

            ref_n, alt_n = 0,0

            for base in querybase:
                if base.upper() == ref_base:
                    ref_n += 1
                else:
                    alt_n += 1

            all_n = ref_n + alt_n

            if all_n > 0:
                alt_freq = round(alt_n/all_n,3)
            else:
                alt_freq = "NA"

            if alt_freq == "NA":
                pass

            val = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t_chr,col,ref_base,ref_n,alt_n,all_n,alt_freq)
            of.write(val)


if __name__ == "__main__":
    main()







