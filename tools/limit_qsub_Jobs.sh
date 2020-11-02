runsh=$1
#JobIDer=$2
Num=15


cat $runsh | while read line;
do
	while true;do
		# check runing jobs and waiting jobs
		#nJobs=$(qstat|grep $JobIDer|wc -l)
		nJobs=$(qstat|wc -l)
		#echo $num
		if [ $nJobs -lt $Num ];then
			d=$(dirname $line)
			cd $d
			qsub -l "mem=5g" -q all $line
			break
		else
			sleep 120
		fi
	done
done
