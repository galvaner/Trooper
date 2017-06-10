__author__ = 'Rasto'

#------------------
#Script description
#Script takes on its input two numbers from and to and pdb file name
#and prepare arguments for rna_denovo_setup.py
#eg. rna_denovo_setup.py -fasta 20.fasta -s 20gap.pdb -fixed_stems -working_res 5-6 8 12-14 from-to ...
#Note: All residues which are in given pdb HAS TO be in option working_res
#!!!!CHANGE!!!!
#script takes all files from ../files/prepare_rosetta_script/ which have to be named eg 100-105.pdb (100, 105 included)
#and prepare scripts for them
#also uses file big_gap.pdb to load names into array
#also creates script for the whole chain prediction (except those residues gaps  predicted separately) an in case of necessity divides it into 300 residues - sections
#------------------

import argparse
from Bio.PDB import *
import os
import PredictSecondaryStructure



def prepare_rosetta_script(structures_number, chainID, target_length, dir, target, template ):
    secStrObject = PredictSecondaryStructure.RunSecStrPRediction(target, template, dir) #  predict secondary structure
    # "end" line HAS TO be at the file end
    # load gaps into array and last residue, also shifts beginning and ending gaps
    file = open(dir+"cut_spheres/big_gaps.txt",'r')
    line = file.readline()
    excluded_gaps_array = []
    while line != "end\n":
        excluded_gaps_array.append([int(line.split()[0])-1, int(line.split()[1])+1])
        line = file.readline()
    '''last_residue_nm = int(file.readline())'''
    print excluded_gaps_array
    # divide the whole sequence into smaller
    divide_into_smaller_seq_and_prepare_statements(excluded_gaps_array,
                                                   target_length,
                                                   300,
                                                   dir+"conserved_edited_regions.pdb",
                                                   structures_number,
                                                   "../target.fasta",
                                                   dir+"prepare_rosetta_scripts/",
                                                   dir+"prepared_statements.txt",
                                                   chainID,
                                                   0,
                                                   secStrObject)
    prepare_statements_for_cut_spheres(dir+"prepared_statements.txt",
                                       "../target.fasta",
                                       structures_number,
                                       dir+"cut_spheres/",
                                       chainID,
                                       0,
                                       secStrObject)

class SelectResidues(Select):
    def __init__(self, accept_boundaries, chain_id):
        self.accept_boundaries = accept_boundaries
        self.chain_id = chain_id

    def accept_residue(self, residue):
        is_conserved = True
        if residue.id[1] >= self.accept_boundaries[0] and residue.id[1] <= self.accept_boundaries[1]:
            return 1
        else:
            return 0

def divide_into_smaller_seq_and_prepare_statements(excluded_gaps_array,
                                                   target_length,
                                                   new_seq_length=300,
                                                   pdb_in="../files/conserved_edited_regions.pdb",
                                                   structures_to_generate="50",
                                                   fasta_name="../target.fasta",
                                                   path="../files/prepare_rosetta_scripts/",
                                                   statement_file="../files/prepared_statements.txt",
                                                   chainID="A",
                                                   modelID=0,
                                                   secStrObject = None):
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('target', pdb_in)
    model = structure[modelID]
    chain = model[chainID]
    break_ = new_seq_length
    last_break = 1
    last_res_id = 0
    section_array = []
    for residue in chain:
        if residue.id[1] >= break_ and residue.id[1] == last_res_id + 1:
            section_array.append([last_break, last_res_id])
            last_break = residue.id[1]
            break_ = break_ + new_seq_length   #decides subsection length
        last_res_id = residue.id[1]
    if last_res_id < target_length:
        last_res_id = int(target_length)
    section_array.append([last_break,last_res_id])
    print "Section array: " + str(section_array)

    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('template', pdb_in)
    io = PDBIO()
    io.set_structure(structure)
    for section in section_array:
        io.save(path + str(section[0]) + "_" + str(section[1]) + ".pdb", select=SelectResidues(section, chainID))

    file_for_statements = open(statement_file,'w')

    for section in section_array:
        first_gap = True
        script_params = "rna_denovo_setup.py -fasta " + fasta_name + " -s " + str(section[0]) + "_" + str(section[1]) + ".pdb -fixed_stems -nstruct " + str(structures_to_generate) + " -secstruct_file ../secstruForRosetta.secstr" + " -working_res " # be aware that one shell scripts extract path from this file...
        workingResPairs = []
        last_gap = section
        for gap in excluded_gaps_array:
            if gap[0] >= section[0] and gap[1] <= section[1]:
                if first_gap:
                    first_gap = False
                    if section[0] == gap[0]: # if's are here because of error later if in working residues is this pattern: 111-111
                        script_params = script_params + str(section[0]) + " "
                        workingResPairs.append([section[0], section[0]])
                    else:
                        script_params = script_params + str(section[0]) + "-" + str(gap[0]) + " "
                        workingResPairs.append([section[0], gap[0]])
                else:
                    if last_gap[1] == gap[0]:
                        script_params = script_params + str(last_gap[1]) + " "
                        workingResPairs.append([last_gap[1], last_gap[1]])
                    else:
                        script_params = script_params + str(last_gap[1]) + "-" + str(gap[0]) + " "
                        workingResPairs.append([last_gap[1], gap[0]])
                last_gap = gap
        if first_gap: #case that there are no gaps in this section (so all residues are working)
            script_params = script_params + str(section[0]) + "-" + str(section[1]) + " "
            workingResPairs.append([section[0], section[1]])
        else:         #last part of working residues
            if last_gap[1] == section[1]:
                script_params = script_params + str(last_gap[1]) + " "
                workingResPairs.append([last_gap[1], last_gap[1]])
            else:
                script_params = script_params + str(last_gap[1]) + "-" + str(section[1]) + " "
                workingResPairs.append([last_gap[1], section[1]])
        print "WORKING_RES_PAIR" + str(workingResPairs)
        addToMatchPairing = secStrObject.FixUnpairedChosenWorkingResidues(workingResPairs)  # secondary structure purpouse
        print "RESULT: " + str(addToMatchPairing)
        script_params += addToMatchPairing
        file_for_statements.write(script_params + "\n")
    file_for_statements.close()

def prepare_statements_for_cut_spheres(statement_file="../files/prepared_statements.txt",
                                       fasta="../target.fasta",
                                       structures_to_generate="100",
                                       path_to_output_cut_spheres='../files/cut_spheres/',
                                       chainID='A',
                                       modelID=0,
                                       secStrObject = None): #prepare script calls for longer gaps created in script cut spheres
    array_index = 0
    file_for_statements = open(statement_file,'a')
    for input_pdb in os.listdir(path_to_output_cut_spheres):
        if input_pdb == "big_gaps.txt" or input_pdb == "statistic.txt":
            continue
        from_ = input_pdb.split("_")[0]
        to_ = input_pdb.split("_")[1].split(".")[0]
        array_index = array_index + 1
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure('self', path_to_output_cut_spheres + input_pdb)
        model = structure[modelID]
        chain = model[chainID]
        last = -1
        current = -1
        seq_start = -1
        seq_length = 0
        first_cycle = True
        workingResPairs = []
        script_params = "rna_denovo_setup.py -fasta " + fasta + " -s " + input_pdb + " -fixed_stems -nstruct " + str(structures_to_generate) + " -secstruct_file ../secstruForRosetta.secstr" + " -working_res "
        for residue in chain:
            current = residue.id[1]
            if first_cycle:
                last = current
                seq_start = current
                first_cycle = False
                continue
            if current == last + 1:
                last = last + 1 #=current
                seq_length = seq_length + 1
                continue
            else:
                if seq_length == 0:
                    script_params = script_params + str(seq_start) + " "
                    workingResPairs.append([seq_start, seq_start])
                else:
                    script_params = script_params + str(seq_start) + "-" + str(seq_start + seq_length) + " "
                    workingResPairs.append([seq_start, seq_start + seq_length])
                seq_length = 0
                last = current
                seq_start = current
        #residues "stuck" in for cycle
        if seq_length == 0:
            script_params = script_params + str(seq_start) + " "
            workingResPairs.append([seq_start, seq_start])
        else:
            script_params = script_params + str(seq_start) + "-" + str(seq_start + seq_length) + " "
            workingResPairs.append([seq_start, seq_start + seq_length])
        script_params = script_params + from_ + "-" + to_
        workingResPairs.append([int(from_), int(to_)])
        print "WORKING_RES_PAIR" + str(workingResPairs)
        addToMatchPairing = secStrObject.FixUnpairedChosenWorkingResidues(workingResPairs) #  secondary structure purpouse
        print "RESULT: " + str(addToMatchPairing)
        script_params += addToMatchPairing
        file_for_statements.write(script_params + "\n")
    file_for_statements.close()

