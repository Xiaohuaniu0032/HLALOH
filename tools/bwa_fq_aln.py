# aln fq1 & fq2 using bwa
# convert sam to bam and sort bam by coord using samtools
# markdup using sambamba's markdup

# the origin propose is to check why normal sample's HLA BAF is not centered at 0.5
# 2021-1-18

import os
import sys
import datetime

fq1,fq2,name,fa_ref,bwa_bin,samtools_bin,sambamba_bin,outdir = sys.argv[1:]

print("start bwa alignment for %s on %s" % (name,datetime.datetime.now()))

# aln
sam = "%s/%s.sam" % (outdir,name)
cmd = "%s mem -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:illumina\" %s %s %s >%s" % (bwa,name,name,fa_ref,fq1,fq2,sam)
os.system(cmd)

# sam to bam
bam = "%s/%s.bam" % (outdir,name)
cmd = "%s view -b -o %s %s" % (samtools_bin,bam,sam)
os.system(cmd)

# sort by coord
sort_bam = "%s/%s.sort.bam" % (outdir,name)
cmd = "%s sort -o %s %s" (samtools_bin,sort_bam,bam)
os.system(cmd)

# markdup
sort_markdup_bam = "%s/%s.sort.markdup.bam" % (outdir,name)
cmd = "%s markdup -t 6 --tmpdir %s %s %s" % (sambamba_bin,outdir,sort_bam,sort_markdup_bam)
os.system(cmd)

# index
cmd = "%s index %s" (samtools_bin,sort_markdup_bam)
os.system(cmd)

print("end bwa alignment for %s on %s" % (name,datetime.datetime.now()))
