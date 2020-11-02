###cal tumor snp vaf
samtools mpileup -d 8000 -f /data1/database/b37/human_g1k_v37.fasta -l /data1/workdir/fulongfei/git_repo/HLALOH/BAF/889gene.snp.bed /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam >/home/fulongfei/workdir/git_repo/HLALOH/test2/tumor.mpileup
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF/pileup2vaf.py /home/fulongfei/workdir/git_repo/HLALOH/test2/tumor.mpileup /home/fulongfei/workdir/git_repo/HLALOH/test2/tumor.snp.vaf
Rscript /data1/workdir/fulongfei/git_repo/HLALOH/BAF/plot_vaf_by_chr.r /home/fulongfei/workdir/git_repo/HLALOH/test2/tumor.snp.vaf 20091701T /home/fulongfei/workdir/git_repo/HLALOH/test2

###extract HLA reads
samtools view -H /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam >/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam 6:29909037-29913661 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam 6:31321649-31324964 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam 6:31236526-31239869 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.sam
samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.sam
samtools sort -n /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.hla.bam -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.sort_by_name.bam
bedtools bamtofastq -i /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.sort_by_name.bam -fq /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.chr6region.1.fastq -fq2 /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.chr6region.2.fastq

###HLA Typing
/home/fulongfei/miniconda3/envs/py27/bin/python2 /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/OptiTypePipeline.py -i /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701N.chr6region.2.fastq --dna -v -o /home/fulongfei/workdir/git_repo/HLALOH/test2 -c /data1/workdir/fulongfei/git_repo/HLALOH/OptiType-1.3.2/config.ini.example
cp /home/fulongfei/workdir/git_repo/HLALOH/test2/*/*_result.tsv /home/fulongfei/workdir/git_repo/HLALOH/test2/hla.result.raw
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/translate_HLA_format.py -i /home/fulongfei/workdir/git_repo/HLALOH/test2/hla.result.raw -o /home/fulongfei/workdir/git_repo/HLALOH/test2/hla.result.new

###get hla fasta
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/tools/get_hla_fasta.py -i /home/fulongfei/workdir/git_repo/HLALOH/test2/hla.result.raw -o /home/fulongfei/workdir/git_repo/HLALOH/test2/patient.hlaFasta.fa

###extract HLA reads from tumor BAM file
samtools view -H /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam >/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam 6:29909037-29913661 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam 6:31321649-31324964 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.sam
samtools view /home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam 6:31236526-31239869 >>/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.sam
samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.sam
samtools sort -n /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.hla.bam -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.sort_by_name.bam
bedtools bamtofastq -i /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.sort_by_name.bam -fq /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.1.fastq -fq2 /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.2.fastq

###aln fq to hla ref fasta using novoalign
/data1/workdir/fulongfei/git_repo/HLALOH/novoDir/novocraft/novoindex /home/fulongfei/workdir/git_repo/HLALOH/test2/patient.hlaFasta.nix /home/fulongfei/workdir/git_repo/HLALOH/test2/patient.hlaFasta.fa
/data1/workdir/fulongfei/git_repo/HLALOH/novoDir/novocraft/novoalign -d /home/fulongfei/workdir/git_repo/HLALOH/test2/patient.hlaFasta.nix -f /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.1.fastq /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.2.fastq -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW >/home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.sam
samtools view -b -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.sam
samtools sort -@ 8 -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csort.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.bam
samtools rmdup /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csort.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csort.noduplicates.bam
samtools view -f 2 -b -o /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csort.noduplicates.bam
samtools index /home/fulongfei/workdir/git_repo/HLALOH/test2/20091701T.chr6region.patient.reference.hlas.csorted.noduplicates.filtered.bam

###baf2loh main script
/home/fulongfei/miniconda3/bin/python3 /data1/workdir/fulongfei/git_repo/HLALOH/BAF2LOH.py -indir /home/fulongfei/workdir/git_repo/HLALOH/test2 -tname 20091701T
