if [ ! -d $PWD/BB20062440_CL01167 ];then
	mkdir $PWD/BB20062440_CL01167
fi

python ../CNVscan.py -bam /data2/Projects/genetic_tumor_338/200627_A00262_0472_AHCWFWDSXY/02_aln/BB20062440_CL01167.rmdup.bam -bed /data1/workdir/fulongfei/genetic/GeExCNV_db/338/bed/338genes_20191225.bed -n BB20062440_CL01167 -m ref -od $PWD/BB20062440_CL01167
