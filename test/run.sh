nbam='/home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam'
tbam='/home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam'
bamdir='/home/fulongfei/workdir/hla_loh/program_test/downsample_bam'
bed='/home/wangce/workdir/database/humandb/panel/889genes_20191225.bed'
gatkdir='/home/fulongfei/workdir/git_repo/HLALOH/gatkDir/picard-tools-1.119'
ref='/home/fulongfei/workdir/git_repo/HLALOH/DB/cnv_ref/889/DB'

/home/fulongfei/miniconda3/bin/python3 ../predictHLALOH.py -nbam $nbam -tbam $tbam -bed $bed -nname 20091701N -tname 20091701T -ref $ref -p 889 -od $PWD
