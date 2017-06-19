WHAT IS THIS PROJECT ABOUT

TOOLS USED
 0. UNIX-based operating system
 1. Rosetta (https://www.rosettacommons.org/)
 2. Emboss (http://emboss.sourceforge.net/what/)
 3. VienaRNA (https://www.tbi.univie.ac.at/RNA/), only if you want to use version witch firstly predicts secondary structure.
 4. PyMOL, only if you want to compare predicted and native structure at the end of prediction.
 5. Python v2.7
 6. BioPython

INSTALATION
 1. Download RNA_PREDICTOR/Predictor folder
 2. Download Rosetta at (https://www.rosettacommons.org/software/license-and-download). It is free for academic users, I recommend to       download alredy built version for UNIX - then it is only copy and run. Be prepared that aproximetly 18 GiB diskspace is needed.         Place built Rosetta software to "RNA_PREDICTOR/rosetta/rosetta_bin_linux_2015.38.58158_bundle" (important thig is that "rosetta"       folder has to be in the same folder as "Predictor" folder and "rosetta/rosetta_bin_linux_2015.38.58158_bundle" is the path where       compiled rosetta software is located.) 
 3. Easiest way to run our tool is to use VO Metacentrum (https://metavo.metacentrum.cz/en/index.html) which is also free for academic       users. In this case only Rosetta download is required - other tools are alredy available in Metacentrum.
 4. If you will use your own UNIX-based system, you will need to modify some bash scripts. 
    First issue is to replace module import statements (imports Python, Emboss, BioPython, Pymol, VienaRNA to Metacentrum). 
    Secondly "qsub" command and necessary code around (This commant puts generated task for FARFAR prediction to queue. Task will be         executed later and many tasks can be executed parallel in Metacentrum. Prediction can generate many tasks, and each task can run up     to one day or 100 generated samples according default configuration. Time depends on the length of predicted gap and number of         known residues around.). Small predictions (up to 300 residues with similarity of template and target sequences) you can even run       without any paralelization. If you want to predict more structures or you have very long structures the best strategy is to use         parallel execution of FARFAR tasks.

How TO RUN PREDICTION
 1. Supply "pairs" file(s) (named like "60-75s30-45g_51_100l.txt", only numbers and letters can be changed) to pair folder. (Prepared       examples can be found in "RNA_Predictor/MyTools/SimilarityOfFastas/pairs/*/")
 2. Make sure that you have desired PDB and FASTA files (in case you provided your own file into pairs folder which coud contain RNA         structures which is not present in fastas and pdbs folders).
 3. Run prepare_rosetta_prediction.sh
 4. Run startFARFAR.sh
 5. Wait till FARFAR prediction is done (default 23:59:00, can be changed in config.py)
 6. Run extract_pdbs.sh
 7. Run concat_pdbs.sh
 
In case of any questions please contact Rastislav Galvanek at r.galvanek@outlook.sk
