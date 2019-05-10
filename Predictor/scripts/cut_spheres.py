__author__ = 'Rasto'
#------------------
#Script description:
#takes pdb from edit_conserved_regions and looks for gaps longer than 8 residues
#then for each that gap computes all residues which have an atom in sphere with radius r
#radius differs on the size of gap NOW set by BULGARIAN CONSTANT
#   <8  residues = 18
#   <15 residues = 28
#   <22 residues = 45
#   <30 residues = 60
#Changed ->gap length is multiplied by magic constant 2.3
# based on measuring of molecule, size cca 200x120x90
#TODO: play with the radii and lower gap-size bound
#output is set of pdb files named eg. 20_31.pdb(number of residues for modeling)
#------------------

import math
from Bio.PDB import *
import numpy as np

# misterious constant is constant which says how each residue is aprox long (take the middle between two residues nad to the sphere take every residues which lies in #missingRes*mysteriousConstant)
def cut_spheres(input_pdb, target_length,  path_to_output, chainID, mysterious_constant=2.5, minimum_gap_length = 5, modelID = 0):
    print "CUT_SPHERES: start"
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('self', input_pdb)
    model = structure[modelID]
    chain = model[chainID]
    gaps = [] #(from(including), to(including), gap_length, last_res_before_gap, first_res_after_gap)
    last = 0
    current = 0
    gap = 0
    first = True
    last_residue_index = 0
    gaps_for_statistics = []
    #prepare array with gaps
    for residue in chain:
        last_residue_index = residue.id[1]
        if first:
            first = False
            last = residue.id[1]
            current = residue.id[1]
            if first != 1:
                gaps_for_statistics.append(current-1)
            continue
        current = residue.id[1]
        if current != last + 1:
            gap = current - last + 1 - 2
            gaps_for_statistics.append(gap)
            if gap > minimum_gap_length:                     #set the min limit for gap size
                gaps.append((last + 1, current - 1, gap, chain[(' ', last, ' ')], chain[(' ', current, ' ')]))
        last = current

    if current < target_length:
        gaps_for_statistics.append(target_length - current)

    file = open(path_to_output+"big_gaps.txt",'w')
    for gap in gaps:
        file.write(str(gap[0]) +  " " + str(gap[1]) + "\n")
    file.write("end\n" + str(last_residue_index))
    file.close()

    middles = []
    #get the middles of sphere
    for gap_record in gaps:
        try:
                atom1_coord = gap_record[3]['C4'].get_vector()
                atom2_coord = gap_record[4]['C4'].get_vector()
        except:
                atom1_coord = gap_record[3]['P'].get_vector()
                atom2_coord = gap_record[4]['P'].get_vector()
        atom_coord_middle = ((atom1_coord[0]+atom2_coord[0])/2, (atom1_coord[1]+atom2_coord[1])/2, (atom1_coord[2]+atom2_coord[2])/2)
        middles.append(atom_coord_middle)

    #array where the residues for delete are marked
    choose_residue = []
    for i  in range(0,last_residue_index+100):
        choose_residue.append(False)
    print "Last residue number from conserved_edited_template.pdb: " + str(last_residue_index)
    print "CUT_SPHERES: end"

    for i in range(0,middles.__len__()):
        #accepted = 0 #debug
        total = 0 #debug
        #in_if = 0 #debug
        #atoms = 0 #debug
        choose_residue = []
        for j  in range(0,last_residue_index+1):
            choose_residue.append(False)
        for residue in chain:
            total = total + 1
            for atom in residue:
                coords = atom.get_vector()
                if(math.sqrt(pow((middles[i][0]-coords[0]),2) +  pow((middles[i][1]-coords[1]),2) + pow((middles[i][2]-coords[2]),2))) <= gaps[i][2]*mysterious_constant: #multiply gap length by TODO: bulgarian constant :D
                    #atoms = atoms + 1
                    #in_if = in_if + 1
                    parent_res = atom.get_parent()
                    choose_residue[parent_res.id[1]] = True
                    #accepted = accepted + 1
                    break
        #print "accepted = " + str(accepted)
        #print "total = " + str(total)
        #print "in if = " + str(in_if)
        #print "atoms = " + str(atoms)
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure('template', input_pdb)
        io = PDBIO()
        io.set_structure(structure)
        io.save(path_to_output+str(gaps[i][0])+"_"+ str(gaps[i][1]) +".pdb", select=SelectResidues(choose_residue, chainID))
    file = open(path_to_output+"statistic.txt",'w')
    file.write("GFS: " + str(gaps_for_statistics)+ "\n")
    file.write("Mean: " + str(np.mean(gaps_for_statistics))+ "\n")
    file.write("Median: " + str(np.median(gaps_for_statistics))+ "\n")
    file.write("Std:" + str(np.std(gaps_for_statistics))+ "\n")
    file.write("MAX: " + str(np.amax(gaps_for_statistics))+ "\n")
    file.write("Gaps count: " + str(len(gaps_for_statistics))+ "\n")
    file.close()


class SelectResidues(Select):
    def __init__(self, residues_bool, chain_id):
        self.residues_bool = residues_bool
        self.chain_id = chain_id

    def accept_residue(self, residue):
        is_conserved = True
        if self.residues_bool[residue.id[1]]:
            return 1
        else:
            return 0


