#!/bin/bash

#get the first line and then pick the chainID into variable (it is on 22. position)
#pattern="../prediction/*/*/results/*/*.pdb"
#files=( $pattern )
#read -r line<${files[0]}
#chainID=${line:21:1}

#call python script and give him chainID as parameter

module add python-2.7.10-gcc python-2.7.10-intel
module add python27-modules-gcc
#module add emboss-6.5.7

for D in `find ../prediction/*/*/results/ -maxdepth 0 -mindepth 0 -type d`
do
	chainID="A"    #`echo $D | cut -d_ -f4`  rosetta dava vzdy chain A?????? 
	python conect_predicted_pdbs.py -i "$chainID" -d "$D"
done

#get native pdb structure into results files

for D in `find ../prediction/*/*/results/ -maxdepth 0 -mindepth 0 -type d`
do
	chainID=`echo $D | cut -d_ -f4`
	native_struct=`echo $D | cut -d_ -f3 | cut -d/ -f2 | tr '[:upper:]' '[:lower:]'`
	native_pdb=../pdbs/$native_struct.pdb
	grep "^ATOM.................$chainID" $native_pdb > $D/$native_struct.pdb
done

#compute RMDS using pymol

module add pymol-1.7.6-gcc
module add pymol-1.8.2.1-gcc

for D in `find ../prediction/*/*/results/ -maxdepth 0 -mindepth 0 -type d`
do
	native_struct=`echo $D | cut -d_ -f3 | cut -d/ -f2 | tr '[:upper:]' '[:lower:]'`
	pml_script="	load $D$native_struct.pdb \n 
			load ${D}predicted_structure.pdb \n 
			align $native_struct, predicted_structure, cycles=0 " 
	echo -e $pml_script > temp.pml
	pymol -qc temp.pml > ${D}RMSD_from_pymol.txt
done