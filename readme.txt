Prediction functionality is in Predictor folder. 
It needs pdb of template, fastas of template and target and file with target-template pairs as input.
Pdbs and fastas are included (downloaded from protein data bank in 1/2016).
For this to work you hawe to change all absolute paths in scripts (paths for rosetta python tools), 
have placed Rosetta instalation as described in rosetta folder,
have Python v2.7, BioPython, emboss-6.5.7 (for prediction with secondary structure also viennaRNA-2.0.7) available form command line.
Also all *.sh script prepared are made to be used in Metacentrum environment - this means it creates tasks which are comuted later (FARFAR prediction).

Then you firstly run prepare_rosetta_prediction.sh, secondly run startFARFAR.sh which is currently set to run 24 hours after each task starts.

When FARFAR is done run concat_pdbs.sh which will concat partly predicted pdbs, and it will compare them to downloaded pdbs and save the RMSD of this comparasion. 
