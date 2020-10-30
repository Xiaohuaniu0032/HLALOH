nbam='/home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701N.sort.bam'
tbam='/home/fulongfei/workdir/hla_loh/program_test/downsample_bam/20091701T.sort.bam'

/home/fulongfei/miniconda3/bin/python3 ../predictHLALOH.v2.py -tbam $tbam -nbam $nbam -tname 20091701T -nname 20091701N -od $PWD -p 889
