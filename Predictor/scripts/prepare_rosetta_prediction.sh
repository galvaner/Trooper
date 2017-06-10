#!/bin/sh
#import modules
echo "Importing modules..."
module add python-2.7.10-gcc python-2.7.10-intel
module add python27-modules-gcc
module add emboss-6.5.7
echo "Modules imported!"

# start prepare rosetta prediction for files in pairs
echo "Starting prediction loop..."
for pairs in ../pairs/*; do
	filename=$(echo $pairs | cut -d/ -f3 | cut -d. -f1)
	echo "Starting python script for $filename" 
	python connect_scripts.py -f $pairs
	rm -rf ../prediction/$filename
	mkdir ../prediction/$filename
	mv ../prediction/????_?_????_? ../prediction/$filename
done

