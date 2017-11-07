## TROOPER (Template Based RNA Predictor)

RNA tertiary structure prediction approaches can be divided into two groups: de novo methods and template-based modeling. De novo are applicable only for small molecules while in case of medium and large size RNA molecules, template-based modeling needs to be employed. While this type of modeling is quite common in protein structure prediction field, there exist only very few tools for template-based RNA structure prediction. Therefore, we present a methodology for prediction of RNA three dimensional structure (target) utilizing a known secondary a tertiary structure of a related RNA molecule (template).

## HOW TO RUN LOCALY:
The software is still in development, but can be run on local pc. There are prepared data for single test prediction. Running this software is not recomended on local machine, expect for finding out how does it work or in case of simple predictions (software is designed to run many parallel predictions at one time, which makes the prediction of large structures reasonably fast). 

## Requirements
1. UNIX-based operating system
2. Rosetta (https://www.rosettacommons.org/software/license-and-download)
3. Emboss (http://emboss.sourceforge.net/what/) (must be accesible from command line)
4. VienaRNA (https://www.tbi.univie.ac.at/RNA/), only if you want to use version witch firstly predicts secondary structure.(must be accesible from command line)
5. PyMOL (must be accesible from command line)
6. Python v2.7 (must be accesible from command line)
7. BioPython

## INSTALATION
 1. Download run_on_local_pc branch.
 2. Install all required programs.
 3. Go to Predictor/scripts/config.py and modify path to configPathToRosetta and configPathToRNATools 

## How TO RUN PREDICTION
 There are prepared data for one prediction, so if you do not want predict your own structures pass to step 3. Also it is recomended to run this prediction to ensure that everything works fine.
 1. Supply "pairs" file(s) (named like "60-75s30-45g_51_100l.txt", only numbers and letters can be changed) to pair folder. (Prepared examples can be found in "RNA_Predictor/MyTools/SimilarityOfFastas/pairs/*/")
 2. Add template and target PDB and FASTA files (target PDB is used at the end to compute RMSD of the prediction) and template secondary structure. Use name convention as you can see for the test prediction.
 3. Run prepare_rosetta_prediction.sh (check if prediction folder was created)
 4. Run startFARFAR.sh
 5. Wait till FARFAR prediction is done (for test prediction about 1 hour) (in prediction folder check rosetta subfolder if it contains predicted substructures).
 6. Run extract_pdbs.sh
 7. Run concat_pdbs.sh (check results in result subfolder)
 
In case of any questions please contact Rastislav Galvanek at r.galvanek@outlook.sk
