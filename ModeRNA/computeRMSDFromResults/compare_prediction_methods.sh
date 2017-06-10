#!/bin/bash
for D in `find ./prediction/*/*/ -maxdepth 0 -mindepth 0 -type d`
do
	line=""
	pred_name=`echo $D | cut -d / -f4`
	#check if the file was predicted by modeRNA and process result
	if [ -a ${D}RMSD_from_pymol.txt ];
	then
		line=`cat ${D}RMSD_from_pymol.txt | grep 'Executive: RMS =' | cut -d '=' -f2 | cut -d '(' -f1`
		line="${pred_name} ${line} MEDZERA"
	else
		line="${pred_name} NOT_PREDICTED MEZDERA"
	fi
	#check if the file was predicted by US and process result
	path="/auto/praha1/galvaner/tst_comparasion/${pred_name}/results/RMSD_from_pymol.txt"
	if [ -a ${path} ];
	then
		tmp=`cat ${path}| grep 'Executive: RMS =' | cut -d '=' -f2 | cut -d '(' -f1`
		line="${line}${tmp}"
		echo ${line} >> results.txt	
	else
		line="${line} NOT_PREDICTED"
		echo ${line} >> results.txt	
	fi
	
done