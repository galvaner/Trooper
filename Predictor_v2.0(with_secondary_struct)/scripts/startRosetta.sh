#!/bin/sh
# install tools
export ROSETTA='/storage/praha1/home/galvaner/rosetta/rosetta_bin_linux_2015.38.58158_bundle/'
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH
source $ROSETTA/tools/rna_tools/INSTALL

# create directories in prediction folder and copy all the files (including script /main_folder/prediction/start_prediction.sh)
cd ../prediction
for filename in ./*/*/files/*/*.pdb; do #*
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2
	mkdir -p rosetta
	cd rosetta
	dirname=$(echo $filename | cut -f 6 -d '/' | cut -f 1 -d '.')
	mkdir -p ./$dirname
	cp ../../../$filename ./$dirname
	cp ../../../start_prediction.sh .
	chmod +x start_prediction.sh
	#rm $filename
	cd ../../..
done

for filename in ./*/*/files/*.txt; do
	predictionDirname1=$(echo $filename | cut -f 2 -d '/')
	predictionDirname2=$(echo $filename | cut -f 3 -d '/')
	cd $predictionDirname1/$predictionDirname2/rosetta
	cp ../../../$filename .
	cd ../../..
done

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

for D in `find ./?*/?*/rosetta -maxdepth 0 -mindepth 0 -type d`
do
	cd "$D"
	pwd
	#result=${PWD##*/}
	#echo $result
	echo "++++++++start_prediction.sh+++++++++"
	./start_prediction.sh
	echo "--------start_prediction.sh---------"
	cd ../../..
done