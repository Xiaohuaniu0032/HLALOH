nbam='/data3/Projects/panel889/201027_A00869_0321_AHKJM2DSXY/02_aln/201023033N.rmdup.recal.bam'
tbam='/data3/Projects/panel889/201027_A00869_0321_AHKJM2DSXY/02_aln/201023033T.rmdup.bam'
nname='201023033N'
tname='201023033T'


/home/fulongfei/miniconda3/bin/python3 ../../predictHLALOH.v3.py -tbam $tbam -nbam $nbam -nname $nname -tname $tname -od $PWD
