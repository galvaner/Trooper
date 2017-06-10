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
	qsub start_script.sh -l nodes=1:ppn=6,mem=16gb,walltime=23h
	cd "$start_dir"
done
