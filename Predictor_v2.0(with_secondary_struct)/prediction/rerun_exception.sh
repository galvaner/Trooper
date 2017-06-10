#!/bin/sh
start_dir=`pwd` 
for F in `find -L ./*/*/*/*/core -type f`
do
	#echo $F
	dir=`echo "$F" | cut -d "/" -f1-5`
	cd "$dir"
	rm ./start_script.sh.e*
	rm ./start_script.sh.o*
	rm ./core
	qsub -l select=1:ncpus=4:mem=16gb -l walltime=23:59:00 start_script.sh
	cd "$start_dir"
done
