###cal tumor snp vaf
###extract HLA reads
samtools view -H /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam >/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam 6:29909000-29914000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam 6:31321000-31326000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam 6:31236000-31241000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.sam
samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.bam /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.sam
samtools sort -n /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.hla.bam -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.sort_by_name.bam
bedtools bamtofastq -i /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.sort_by_name.bam -fq /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.1.fastq -fq2 /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.2.fastq

###HLA Typing
/home/fulongfei/miniconda3/envs/py27/bin/python2 /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/OptiTypePipeline.py -i /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.2.fastq --dna -v -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample -c /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/config.ini
cp /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/*/*_result.tsv /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/hla.result.raw
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/translate_HLA_format.py -i /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/hla.result.raw -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/hla.result.new

###get hla fasta
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/get_hla_fasta.py -i /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/hla.result.raw -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/patient.hlaFasta.fa -o2 /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/patient.hlaFasta.all.4digit.alleles.fa

###extract HLA reads from normal BAM file
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/get_proper_aln_bam.py /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.2.fastq 201031005N /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/patient.hlaFasta.fa /usr/bin/bwa /home/fulongfei/miniconda3/bin/samtools /home/fulongfei/miniconda3/bin/sambamba /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample

###baf2loh main script for normal
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF2LOH.py -indir /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample -name 201031005N -bam /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005N.chr6region.patient.reference.hlas.sort.markdup.proper_aln.bam -jre /usr/bin/java

###extract HLA reads from tumor BAM file
samtools view -H /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam >/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam 6:29909000-29914000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam 6:31321000-31326000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.sam
samtools view /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam 6:31236000-31241000 >>/home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.sam
samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.bam /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.sam
samtools sort -n /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.hla.bam -o /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.sort_by_name.bam
bedtools bamtofastq -i /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.sort_by_name.bam -fq /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.chr6region.1.fastq -fq2 /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.chr6region.2.fastq
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/get_proper_aln_bam.py /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.chr6region.2.fastq 201031005T /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/patient.hlaFasta.fa /usr/bin/bwa /home/fulongfei/miniconda3/bin/samtools /home/fulongfei/miniconda3/bin/sambamba /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample

###baf2loh main script for tumor
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF2LOH.py -indir /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample -name 201031005T -bam /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample/201031005T.chr6region.patient.reference.hlas.sort.markdup.proper_aln.bam -jre /usr/bin/java



/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/cnv.py /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam 201031005N /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed /data1/database/b37/human_g1k_v37.fasta /data1/workdir/fulongfei/git_repo/HLALOH/DB/cnv_ref/889/DB /home/fulongfei/miniconda3/bin/sambamba /home/fulongfei/miniconda3/bin/bedtools /home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH /home/fulongfei/miniconda3/bin/python3 /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/cnv.py /data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam 201031005T /home/wangce/workdir/database/humandb/panel/889genes_20191225.bed /data1/database/b37/human_g1k_v37.fasta /data1/workdir/fulongfei/git_repo/HLALOH/DB/cnv_ref/889/DB /home/fulongfei/miniconda3/bin/sambamba /home/fulongfei/miniconda3/bin/bedtools /home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/git_repo/HLALOH /home/fulongfei/miniconda3/bin/python3 /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample
perl /data1/workdir/fulongfei/git_repo/HLALOH/tools/determine_HLALOH_by_BAF.pl /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample 201031005N 201031005T /home/fulongfei/workdir/git_repo/HLALOH/test/Negative_Sample
