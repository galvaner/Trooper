## TROOPER (Template Based RNA Predictor)

RNA tertiary structure prediction approaches can be divided into two groups: de novo methods and template-based modeling. De novo are applicable only for small molecules while in case of medium and large size RNA molecules, template-based modeling needs to be employed. While this type of modeling is quite common in protein structure prediction field, there exist only very few tools for template-based RNA structure prediction. Therefore, we present a methodology for prediction of RNA three dimensional structure (target) utilizing a known structure of a related RNA molecule (template). First, the target and template sequences are aligned. Next, sequentially similar regions in the alignment are identified and corresponding substructures are transferred from template to target. The remaining parts of the target structures are predicted using an external tool. This phase includes treatment of indels and valid linking of the transferred and predicted portions of the target structure. Our proposed method is able to predict even large ribosomal RNA structures when sufficiently similar template is available. The experiments have shown that the main impact on the quality of prediction has the sequence similarity of the template and target and number of indels. For structures with size of hundreds of nucleotides with sequence similarity with template over 50% and ratio of indels up to 50% the method is able to generate target structures up to ten RMSD with respect to the reference structure.

## TOOLS USED
1. UNIX-based operating system
2. Rosetta (https://www.rosettacommons.org/)
3. Emboss (http://emboss.sourceforge.net/what/)
4. VienaRNA (https://www.tbi.univie.ac.at/RNA/), only if you want to use version witch firstly predicts secondary structure.
5. PyMOL, only if you want to compare predicted and native structure at the end of prediction.
6. Python v2.7
7. BioPython

## INSTALATION
 1. Download RNA_PREDICTOR/Predictor folder
 2. Download Rosetta at (https://www.rosettacommons.org/software/license-and-download). It is free for academic users, I recommend to download alredy built version for UNIX - then it is only copy and run. Be prepared that aproximetly 18 GiB diskspace is needed. Place built Rosetta software to "RNA_PREDICTOR/rosetta rsetta_bin_linux_2015.38.58158_bundle" (important thig is that "rosetta" folder has to be in the same folder as "Predictor" folder and rosetta/rosetta_bin_linux_2015.38.58158_bundle" is the path where compiled rosetta software is located.) 
 3. Easiest way to run our tool is to use VO Metacentrum (https://metavo.metacentrum.cz/en/index.html) which is also free for academic users. In this case only Rosetta download is required - other tools are alredy available in Metacentrum.
 4. If you will use your own UNIX-based system, you will need to modify some bash scripts. First issue is to replace module import statements (imports Python, Emboss, BioPython, Pymol, VienaRNA to Metacentrum). 
    Secondly "qsub" command and necessary code around (This command puts generated task for FARFAR prediction to queue. Tasks will be executed later and many of them can be executed parallel. Prediction can generate many tasks, and each task can run up to one day or 100 generated samples according default configuration. The time needed for prediction depends on the length of predicted gap and number of known residues around.). Small predictions (up to 300 residues with similarity of template and target sequences over 60% ) you can even run without any paralelization. If you want to predict more structures or you have very long structures the best strategy is to use parallel execution of FARFAR tasks.

## How TO RUN PREDICTION
 1. Supply "pairs" file(s) (named like "60-75s30-45g_51_100l.txt", only numbers and letters can be changed) to pair folder. (Prepared examples can be found in "RNA_Predictor/MyTools/SimilarityOfFastas/pairs/*/")
 2. Make sure that you have desired PDB and FASTA files (in case you provided your own file into pairs folder which coud contain RNA structures which is not present in fastas and pdbs folders).
 3. Run prepare_rosetta_prediction.sh
 4. Run startFARFAR.sh
 5. Wait till FARFAR prediction is done (default 23:59:00, can be changed in config.py)
 6. Run extract_pdbs.sh
 7. Run concat_pdbs.sh
 
In case of any questions please contact Rastislav Galvanek at r.galvanek@outlook.sk
