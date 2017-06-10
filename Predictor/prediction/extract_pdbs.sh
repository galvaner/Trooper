#!/bin/bash
# install tools
export ROSETTA='/storage/praha1/home/galvaner/rosetta/rosetta_bin_linux_2015.38.58158_bundle/'
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH
source $ROSETTA/tools/rna_tools/INSTALL
# extract predicted pdbs from *.out files into result folder
for D in `find ./*/*/results/predicted_parts`
do
	cd "$D"
	for F in `find ../../rosetta/*/*.out`
	do
		../../../../../../rosetta/rosetta_bin_linux_2015.38.58158_bundle/tools/rna_tools/bin/extract_lowscore_decoys.py $F 1
	done
	cd ../../../..
done
