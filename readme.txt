Prediction functionality is in Predictor folder. 
It needs pdb of template, fastas of template and target and file with target-template pairs as input.
Pdbs and fastas are included (downloaded from protein data bank in 1/2016).
For this to work you have to change set paths in config.py (paths for rosetta python tools), 
have placed Rosetta instalation as described in rosetta folder,
have Python v2.7, BioPython, emboss-6.5.7, PyMol (for prediction with secondary structure also viennaRNA-2.0.7) available form command line.
All *.sh script are made to be used in Metacentrum environment (https://wiki.metacentrum.cz) - this means it creates tasks which are comuted later (FARFAR prediction) - if you want to run it on different unix environment you have to modify module imports (Python, emboss, PyMol) and 'qsub' command which runs task in metacentrum.

Be aware that by not runnig the tasks parallel you will lost big advantage of this tool and the prediction will take longer time.

More info how to run the prediction can be found in Predictor/RunPrediction.txt
