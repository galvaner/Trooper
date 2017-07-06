#!/bin/bash
for file in ./pdbs/*pdb; 
do 
	./x3dna-dssr.exe -input="$file" -output="json.temp" --json; cat "json.temp" | python x3dna-dssr.py; 
done