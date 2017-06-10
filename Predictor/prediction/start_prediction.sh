#!/bin/sh
echo "++++++FIRST_PART++++++++"
pwd
statements=$(<prepared_statements.txt)
while read -r line; do
	dir=$(echo $line | cut -f 5 -d ' ' | cut -f 2 -d '/' | cut -f 1 -d '.')
	cd $dir
	line="../../../../../../rosetta/rosetta_bin_linux_2015.38.58158_bundle/tools/rna_tools/bin/$line"
	$line
	cd ..
done <<< "$statements"
echo "++++++SECOND_PART++++++++"
pwd
wt=1
for D in `find ./?* -type d`
do
	wt=$((wt+10))
	cd "$D"
	path=$(pwd)
	echo 	"#!/bin/bash
	 cd $path
	 sleep $wt
	 source README_FARFAR" > start_script.sh
	echo "$D.pdb:"
	qsub start_script.sh -l nodes=1:ppn=6,mem=16gb,walltime=24h
	cd ".."
done

