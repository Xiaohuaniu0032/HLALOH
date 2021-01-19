runsh=$1
Num=$2
mem=$3


cat $runsh | while read line;
do
	while true;do
		# check runing jobs and waiting jobs
		# nJobs=$(qstat|grep $JobIDer|wc -l)
		nJobs=$(qstat|wc -l)
		# echo $num
		if [ $nJobs -lt $Num ];then
			d=$(dirname $line)
			cd $d
			bn=$(basename $line)
			echo $bn
			qsub -l "nodes=1:ppn=12,mem=$mem" -q all $line
			break
		else
			sleep 30
		fi
	done
done
