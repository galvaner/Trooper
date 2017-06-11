#!/bin/bash

# be aware of paths structure, also script expects target structure to be present as target.pdb
# compute RMDS using PyMol

module add pymol-1.7.6-gcc
module add pymol-1.8.2.1-gcc

for D in `find ./prediction/*/*/ -maxdepth 0 -mindepth 0 -type d`
do
	if [ -a ${D}model.pdb ];
	then
		mv ${D}model.pdb ${D}modelModeRNA.pdb
		pml_script="	load ${D}target.pdb \n 
				load ${D}modelModeRNA.pdb \n 
				align target, modelModeRNA, cycles=0 " 
		echo -e $pml_script > temp.pml
		pymol -qc temp.pml > ${D}RMSD_from_pymol.txt
	else
		echo 	"In ${D} the model file is missing."
	fi
done