#!/bin/bash

# script used to compare MY prediction and ModeRNA prediction results (first the script which will run Pymol and get RMSD has to be used)
# be aware of paths - best to use same structure of folders
# output can be used as input for PlotGraphModeRNAvsMyPred.py script

MY_PREDICTION="predictionparalel"
#rm ./result_table.txt
#printf "%9s %9s %6s %8s %8s %11s %12s\n" "TARGET" "TEMPLATE" "LEN" "SIMIL" "GAP" "MY_RMSD" "MODERNA_RMSD" >> ./result_table.txt
for F in `find ./${MY_PREDICTION}/prediction/*/*/results/RMSD_from_pymol.txt`
do
	TARGET=`echo ${F} | cut -d_ -f3,4 | cut -d/ -f2`
	TEMPLATE=`echo ${F} | cut -d_ -f5,6 | cut -d/ -f1`
	MY_RMSD=`tail ${F} -n1 | cut -d= -f2 | cut -d\( -f1`
	GENERAL_PATH=`echo ${F} | cut -d/ -f1-5`		# expects $F to be: ./predictionparalel/prediction/75-90s15-30g_51_100l/486D_C_1I9V_A/results/RMSD_from_pymol.txt
	ALN_FILE=${GENERAL_PATH}/files/emboss.aln
	SIMILARITY=`grep Similarity ${ALN_FILE}| cut -d\( -f2 | cut -d\) -f1`
	GAP=`grep Gaps ${ALN_FILE}| cut -d\( -f2 | cut -d\) -f1`
	LENGTH=`grep Length ${ALN_FILE}| cut -d\: -f2`
	MODERNA_PATH=./modeRNA/prediction/`echo ${F} | cut -d/ -f4,5`/RMSD_from_pymol.txt
	MODERNA_RMSD=`tail ${MODERNA_PATH} -n1 | cut -d= -f2 | cut -d\( -f1`
	re='[0-9]+[.][0-9]+'
	if ! [[ ${MY_RMSD} =~ $re ]] ; then
   		MY_RMSD="N/A"
	fi
	if ! [[ ${MODERNA_RMSD} =~ $re ]] ; then
   		MODERNA_RMSD="N/A"
	fi

	#echo "${TARGET} ${TEMPLATE} ${LENGTH} ${SIMILARITY} ${GAP} ${MY_RMSD} ${MODERNA_RMSD}" #>> ./result_table.txt
	printf "%9s %9s %6s %8s %8s %11s %12s\n" ${TARGET} ${TEMPLATE} ${LENGTH} ${SIMILARITY} ${GAP} ${MY_RMSD} ${MODERNA_RMSD} >> ./result_table.txt
done