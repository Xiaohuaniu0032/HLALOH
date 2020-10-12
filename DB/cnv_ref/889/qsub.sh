resdir=$PWD
ls $resdir/*/*sh | while read line;
do
	d=$(dirname $line)
	cd $d
	qsub -l "mem=2g" -q all $line
done
