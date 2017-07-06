#!/bin/bash
source config.py
# reruns all FARFAR predictionwhich ended with core dump (not exception because it could be walltime exceeded)
start_dir=`pwd` 
echo "++Start rerun exceptions++"
for F in `find -L ../prediction/*/*/*/*/core -type f`
do
	#echo "Rerun in folder: ${F}"
	dir=`echo "$F" | cut -d "/" -f1-5`
	cd "$dir"
	currentDir=`pwd`
	echo "	Currently rerunning: ${currentDir}"
	rm ./start_script.sh.e*
	rm ./start_script.sh.o*
	rm ./core
	qsub -l select=1:ncpus=4:mem=16gb -l walltime=${configWallTime} start_script.sh
	cd "$start_dir"
done
