nbam='/data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701N.sort.bam'
tbam='/data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln/20091701T.sort.bam'
bamdir='/data3/Projects/panel889/200921_A00869_0301_AHKJC3DSXY/02_aln'
bed='/home/wangce/workdir/database/humandb/panel/889genes_20191225.bed'
gatkdir='/home/fulongfei/workdir/git_repo/HLALOH/gatkDir/picard-tools-1.34'
ref='/home/fulongfei/workdir/git_repo/HLALOH/CNVscan/ref'

/home/fulongfei/miniconda3/bin/python3 ../predictHLALOH.py -nbam $nbam -tbam $tbam -bamDir $bamdir -bed $bed -nname 20091701N -tname 20091701T -ref $ref -gatkDir $gatkdir -od $PWD -novoDir /home/fulongfei/workdir/git_repo/HLALOH/novoDir/novocraft
