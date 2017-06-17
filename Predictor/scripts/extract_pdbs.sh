#!/bin/bash
source config.py
# install tools
export ROSETTA="${configPathToRosetta}"
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH
source $ROSETTA/tools/rna_tools/INSTALL
# extract predicted pdbs from *.out files into result folder
startFolder=`pwd`
for D in `find ../prediction/*/*/results/predicted_parts`
do
	cd "${D}"
	for F in `find ../../rosetta/*/*.out`
	do
		${configPathToRNATools}extract_lowscore_decoys.py $F 1
	done
	cd "${startFolder}"
done
