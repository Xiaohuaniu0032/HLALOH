# use Panel of Normal (PoN) method to analysis normal and tumor CNV
import sys
import os
import datetime

bam,name,bed,fa_ref,cnv_ref_dir,sbb_bin,bedtools_bin,rscript,bin_dir,py3_bin,outdir = sys.argv[1:]

print("start cnv analysis for %s on %s" % (name,datetime.datetime.now()))
# cal depth
cmd = "perl %s/CNVscan/bin/cal_depth.pl -bam %s -n %s -bed %s -sbb %s -outdir %s" % (bin_dir,bam,name,bed,sbb_bin,outdir)
os.system(cmd)

# add gc
depth_f = "%s/%s.targetcoverage.cnn" % (outdir,name)
cmd = "perl %s/CNVscan/bin/add_gc.pl -bed %s -depth %s -bedtools_bin %s -fa %s -outdir %s" % (bin_dir,bed,depth_f,bedtools_bin,fa_ref,outdir)
os.system(cmd)

# infer sex
cmd = "%s %s/CNVscan/bin/infer_sex.py -cov %s -o %s/%s.sex.txt" % (py3_bin,bin_dir,depth_f,outdir,name)
os.system(cmd)

# gc correct
add_gc_depth = "%s/%s.targetcoverage.cnn.with.gc.xls" % (outdir,name)
sex_file = "%s/%s.sex.txt" % (outdir,name)
gc_correct_depth = "%s/%s.targetcoverage.cnn.with.gc.xls.gc.corrected.xls" % (outdir,name)

cmd = "%s %s/CNVscan/bin/gc_correct.r %s %s %s" % (rscript,bin_dir,add_gc_depth,sex_file,gc_correct_depth)
os.system(cmd)

# library norm
cmd = "perl %s/CNVscan/bin/normalize.pl -d %s -od %s" % (bin_dir,gc_correct_depth,outdir)
os.system(cmd)

# make ref matrix
ref_mat = "%s/ref.matrix.txt" % (outdir)
cmd = "perl %s/CNVscan/bin/make_ref_matrix.pl %s %s" % (bin_dir,cnv_ref_dir,ref_mat)
os.system(cmd)

# cal logR
norm_file = "%s/%s.norm.xls" % (outdir,name)
logR = "%s/%s.logR.xls" % (outdir,name)
cmd = "perl %s/CNVscan/bin/cal_logR.pl %s %s %s" % (bin_dir,norm_file,ref_mat,logR)
os.system(cmd)

# calculate copy number
cnFile = "%s/%s.CopyNumber.xls" % (outdir,name)
cmd = "perl %s/CNVscan/bin/Gene_Level_CNV.pl -in %s -o %s" % (bin_dir,logR,cnFile)
os.system(cmd)


print("end cnv analysis for %s on %s" % (name,datetime.datetime.now()))
