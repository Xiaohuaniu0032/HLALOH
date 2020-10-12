resdir=$PWD
bed='/home/wangce/workdir/database/humandb/panel/889genes_20191225.bed'

cat bam.list | while read line;
do
	name=$(echo $(basename $line) | cut -d '.' -f 1)
	if [ ! -d $resdir/$name ];
	then
		`mkdir $resdir/$name`
	fi

	python /home/fulongfei/workdir/git_repo/HLALOH/CNVscan/CNVscan.py -bam $line -bed $bed -n $name -m ref -od $resdir/$name
done
