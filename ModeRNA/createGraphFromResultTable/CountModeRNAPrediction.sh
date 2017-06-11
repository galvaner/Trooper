#!/bin/bash

# compute how many files were predicted/unpredicted by ModeRNA 

predicted='0'
unpredicted='0'
total='0'
for D in `find ./prediction/*/*/ -maxdepth 0 -mindepth 0 -type d`
do
	line=""
	pred_name=`echo $D | cut -d / -f4`
# check if the file was predicted by modeRNA and process result
	total=$((total+1))
	if [ -a ${D}RMSD_from_pymol.txt ];
	then
		predicted=$((predicted+1))
		line=`cat ${D}RMSD_from_pymol.txt | grep 'Executive: RMS =' | cut -d '=' -f2 | cut -d '(' -f1`
		line="${pred_name} ${line}"
		echo ${line}
	else
		unpredicted=$((unpredicted+1))
		line="${pred_name} NOT_PREDICTED MEZDERA"
	fi
	#echo ${line}
done
echo ${predicted}
echo ${unpredicted}
echo ${total}