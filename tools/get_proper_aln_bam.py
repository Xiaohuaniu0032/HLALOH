import sys
import os

fq1,fq2,name,fa,bwa,samtools_bin,sambamba_bin,outdir = sys.argv[1:]

# index ref
cmd = "%s index %s" % (bwa,fa)
os.system(cmd)

# aln
hla_sam = "%s/%s.chr6region.patient.reference.hlas.sam" % (outdir,name) # normal/tumor
cmd = "%s mem -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:illumina\" %s %s %s >%s" % (bwa,name,name,fa,fq1,fq2,hla_sam)
os.system(cmd)

# sam to bam
hla_bam = "%s/%s.chr6region.patient.reference.hlas.bam" % (outdir,name)
cmd = "%s view -b -o %s %s" % (samtools_bin,hla_bam,hla_sam)
os.system(cmd)

# sort by coord
hla_sort_bam = "%s/%s.chr6region.patient.reference.hlas.sort.bam" % (outdir,name)
cmd = "%s sort -o %s %s" % (samtools_bin,hla_sort_bam,hla_bam)
os.system(cmd)

# markdup
markdup_bam = "%s/%s.chr6region.patient.reference.hlas.sort.markdup.bam" % (outdir,name)
cmd = "%s markdup -t 6 --tmpdir %s %s %s" % (sambamba_bin,outdir,hla_sort_bam,markdup_bam)
os.system(cmd)

# get proper aln reads
proper_aln_bam = "%s/%s.chr6region.patient.reference.hlas.sort.markdup.proper_aln.bam" % (outdir,name)
'''

0x2   PROPER_PAIR
samtools flags DUP,SECONDARY => 0x500 1280    SECONDARY,DUP

'''
#cmd = "%s view -f 0x2 -F 0x500 -b -o %s %s" % (samtools_bin,proper_aln_bam,markdup_bam)
cmd = "%s view -F 0x500 -b -o %s %s" % (samtools_bin,proper_aln_bam,markdup_bam)
os.system(cmd)


# index
cmd = "%s index %s" % (samtools_bin,proper_aln_bam)
os.system(cmd)
 
