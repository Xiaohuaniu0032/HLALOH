perl get_snp_from_bed.pl /data1/workdir/fulongfei/genetic/GeExCNV_db/338/bed/338genes_20191225.bed /data1/software/annovar/humandb/hg19_EAS.sites.2014_10.txt $PWD/338gene.snp.txt

python pileup2vaf.v2.py -bam /data2/Projects/genetic_tumor_338/200613_A00838_0217_AH5TCFDSXY/02_aln/BB20061128_CL01167.rmdup.bam -bed /data1/workdir/fulongfei/genetic/GeExCNV_db/338/bed/338genes_20191225.bed -outfile test.vaf
