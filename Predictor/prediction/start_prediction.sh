#!/bin/sh
# take statements created by python pipeline and run them (statements are for python tool which prepares statements for actual FARFAR prediction)
# be aware of path to correct folder structure regarding rna_tools

echo "++Execute prepared statements (rostetta_rna_tools)++"
#pwd
statements=$(<prepared_statements.txt)
while read -r line; do
	echo "	Currently executed: $line"
	dir=$(echo $line | cut -f 5 -d ' ' | cut -f 2 -d '/' | cut -f 1 -d '.')
	cd $dir
	line="../../../../../../rosetta/rosetta_bin_linux_2015.38.58158_bundle/tools/rna_tools/bin/$line" # relative path
	$line
	cd ..
done <<< "$statements"

# this part is specific for metacentrum - scripts which will run tasks which will run later (planed by metacentrum scheduller)
# there is added a sleep statement to the script which will be executed - reason is that there were collisions in some source file of rosetta 
# (probably both starting prediction tried to access the same file in the same moment which crashed - those sleeps resolved the problem) 

currentDirectory=$(pwd)
echo "++Creating tasks for metacentrum (${currentDirectory})++"
#pwd
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
	qsub -l select=1:ncpus=4:mem=16gb -l walltime=23:59:00 start_script.sh
	cd ".."
done

