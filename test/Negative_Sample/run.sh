nbam='/data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005N.rmdup.recal.bam'
tbam='/data3/Projects/panel889/201102_A00403_0449_AHKY5WDSXY/02_aln/201031005T.rmdup.bam'
nname='201031005N'
tname='201031005T'

/home/fulongfei/miniconda3/bin/python3 ../../predictHLALOH.v3.py -tbam $tbam -nbam $nbam -nname $nname -tname $tname -od $PWD
