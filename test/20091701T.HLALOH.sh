###extract HLA reads
/home/fulongfei/miniconda3/bin/samtools view -H /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam >/home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.sam
/home/fulongfei/miniconda3/bin/samtools view /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam 6:29909037-29913661 >>/home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.sam
/home/fulongfei/miniconda3/bin/samtools view /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam 6:31321649-31324964 >>/home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.sam
/home/fulongfei/miniconda3/bin/samtools view /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam 6:31236526-31239869 >>/home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.sam
/home/fulongfei/miniconda3/bin/samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.bam /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.sam
/home/fulongfei/miniconda3/bin/samtools sort -n /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla.bam -o /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.sort_by_name.bam
/home/fulongfei/miniconda3/bin/bedtools bamtofastq -i /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.sort_by_name.bam -fq /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.chr6region.1.fastq -fq2 /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.chr6region.2.fastq

###HLA Typing
/home/fulongfei/miniconda3/envs/py27/bin/python2 /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/OptiTypePipeline.py -i /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.chr6region.2.fastq --dna -v -o /home/fulongfei/workdir/git_repo/HLALOH/test -c /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/config.ini.example
cp /home/fulongfei/workdir/git_repo/HLALOH/test/*/*_result.tsv /home/fulongfei/workdir/git_repo/HLALOH/test/hla.result.raw
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/translate_HLA_format.py -i /home/fulongfei/workdir/git_repo/HLALOH/test/hla.result.raw -o /home/fulongfei/workdir/git_repo/HLALOH/test/hla.result.new

###cal BAF for normal
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF/pileup2vaf.v2.py -bam /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -outfile /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.normal.vaf

###cal BAF for tumor
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF/pileup2vaf.v2.py -bam /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701T.sort.bam -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -outfile /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.tumor.vaf

### call CNV for normal
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/cal_depth.pl -bam /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam -n 20091701N -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -sbb /home/fulongfei/miniconda3/bin/sambamba -outdir /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/add_gc.pl -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -depth /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.targetcoverage.cnn -bedtools_bin /home/fulongfei/miniconda3/bin/bedtools -fa /data1/database/b37/human_g1k_v37.fasta -outdir /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/infer_sex.py -cov /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.targetcoverage.cnn -o /home/fulongfei/workdir/git_repo/HLALOH/test/sex.txt
/home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/gc_correct.r /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.targetcoverage.cnn.with.gc.xls /home/fulongfei/workdir/git_repo/HLALOH/test/sex.txt /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.targetcoverage.cnn.with.gc.xls.gc.corrected.xls
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/normalize.pl -d /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.targetcoverage.cnn.with.gc.xls.gc.corrected.xls -od /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/make_ref_matrix.pl /home/fulongfei/workdir/git_repo/HLALOH/DB/cnv_ref/889 /home/fulongfei/workdir/git_repo/HLALOH/test/ref.matrix.txt
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/cal_logR.pl /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.norm.xls /home/fulongfei/workdir/git_repo/HLALOH/test/ref.matrix.txt /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.logR.xls

### call CNV for tumor
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/cal_depth.pl -bam /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701T.sort.bam -n 20091701T -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -sbb /home/fulongfei/miniconda3/bin/sambamba -outdir /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/add_gc.pl -bed /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed -depth /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.targetcoverage.cnn -bedtools_bin /home/fulongfei/miniconda3/bin/bedtools -fa /data1/database/b37/human_g1k_v37.fasta -outdir /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/infer_sex.py -cov /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.targetcoverage.cnn -o /home/fulongfei/workdir/git_repo/HLALOH/test/sex.txt
/home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/gc_correct.r /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.targetcoverage.cnn.with.gc.xls /home/fulongfei/workdir/git_repo/HLALOH/test/sex.txt /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.targetcoverage.cnn.with.gc.xls.gc.corrected.xls
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/normalize.pl -d /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.targetcoverage.cnn.with.gc.xls.gc.corrected.xls -od /home/fulongfei/workdir/git_repo/HLALOH/test
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/make_ref_matrix.pl /home/fulongfei/workdir/git_repo/HLALOH/DB/cnv_ref/889 /home/fulongfei/workdir/git_repo/HLALOH/test/ref.matrix.txt
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/CNVscan/bin/cal_logR.pl /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.norm.xls /home/fulongfei/workdir/git_repo/HLALOH/test/ref.matrix.txt /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.logR.xls

###make ascat logR for tumor
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/ASCAT/ascatLogR.pl /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.tumor.vaf /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.logR.xls /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.ascat.logR.xls

###make ascat logR for normal
/home/fulongfei/miniconda3/bin/perl /data1/workdir/fulongfei/git_repo/HLALOH/ASCAT/ascatLogR.pl /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.normal.vaf /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.logR.xls /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.ascat.logR.xls

###estimate tumor purity & ploidy By ASCAT
/home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH/ASCAT/ascat.r /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.ascat.logR.xls /home/fulongfei/workdir/git_repo/HLALOH/test/20091701T.tumor.vaf /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.ascat.logR.xls /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.normal.vaf /home/fulongfei/workdir/git_repo/HLALOH/test/PurityPloidyEst.txt

###detect HLA LOH by lohhla software
/home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH/lohhla/LOHHLAscript.R --patientId 20091701T --outputDir /home/fulongfei/workdir/git_repo/HLALOH/test --normalBAMfile /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam --BAMDir /data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln --hlaPath /home/fulongfei/workdir/git_repo/HLALOH/test/20091701N.hla_alleles --HLAfastaLoc /home/fulongfei/workdir/git_repo/HLALOH/test/patient.hla.fa --CopyNumLoc /home/fulongfei/workdir/git_repo/HLALOH/test/PurityPloidyEst.txt --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir /home/fulongfei/workdir/git_repo/HLALOH/gatkDir/picard-tools-1.34 --novoDir /home/fulongfei/workdir/git_repo/HLALOH/novoDir/novocraft
