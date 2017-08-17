#!/bin/bash

# start prediction of all prepared target-template pairs

source config.py

# install tools
export ROSETTA=${configPathToRosetta}
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH
source $ROSETTA/tools/rna_tools/INSTALL

# create directories in prediction folder and copy all the files (including script /Predictor/prediction/start_prediction.sh)
cd ../prediction

# clear predictions which resulted in error
for D in `find ./?*/?*/files -maxdepth 0 -mindepth 0 -type d`
do
	parentDir=${D::-6}
	if [ ! -f ${D}/prepared_statements.txt ]; then
		# echo ${parentDir}
    		rm -rf ${parentDir}
	fi
done


# copy prepared parts of .pdb files to newly created rosetta folder

for filename in ./*/*/files/*/*.pdb; do 
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2
	mkdir -p rosetta
	cd rosetta
	dirname=$(echo $filename | cut -f 6 -d '/' | cut -f 1 -d '.')
	mkdir -p ./$dirname
	cp ../../../$filename ./$dirname
	cp ../../../../scripts/start_prediction.sh .
	chmod +x start_prediction.sh
	#rm $filename
	cd ../../..
done

# copy prepared parts of preparedStatements.txt files to rosetta folder

for filename in ./*/*/files/*.txt; do
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2/rosetta
	cp ../../../$filename .
	cd ../../..
done

# copy target fasta file to rosetta folder

for filename in ./*/*/files/*.fasta; do
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2/rosetta
	cp ../../../$filename .
	cd ../../..
done

for filename in ./*/*/files/*.secstr; do
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2/rosetta
	cp ../../../$filename .
	cd ../../..
done

# run start_prediction.sh for each target-template pair

for D in `find ./?*/?*/rosetta -maxdepth 0 -mindepth 0 -type d`
do
	cd "$D"
	#pwd
	#result=${PWD##*/}
	#echo $result
	echo "++Run start_prediction.sh++"
	./start_prediction.sh
	cd ../../..
done