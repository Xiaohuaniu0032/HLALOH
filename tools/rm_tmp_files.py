import os
import sys

outdir,nname,tname = sys.argv[1:]

samlist = [nname,tname]


def rm_file(infile):
    if os.path.exists(infile):
        print("find %s file, will be removed." % (infile))
        os.remove(infile)
    else:
        print("not find %s file." % (infile))


for s in samlist:
    print("remove temp files...")
    # rm extract_HLA_reads() tmp files
    hla_sam = "%s/%s.hla.sam" % (outdir,s)
    rm_file(hla_sam)
    
    hla_bam = "%s/%s.hla.bam" % (outdir,s)
    rm_file(hla_bam)

    sort_bam = "%s/%s.sort_by_name.bam" % (outdir,s)
    rm_file(hla_bam)

    # rm get_proper_aln_bam.py tmp files
    sam = "%s/%s.chr6region.patient.reference.hlas.sam" % (outdir,s)
    rm_file(sam)

    bam = "%s/%s.chr6region.patient.reference.hlas.bam" % (outdir,s)
    rm_file(bam)

    sort_bam = "%s/%s.chr6region.patient.reference.hlas.sort.bam" % (outdir,s)
    rm_file(sort_bam)
