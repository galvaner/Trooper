__author__ = 'Rasto'
import os

def get_pdb_files(directory_path):
    file_list = []
    for file in os.listdir(directory_path):
        if file.endswith(".pdb"):
            file_list.append(file)
    return file_list

def get_range_to_copy(pdb_file_name): #strings like 1_10.out.pdb
    str = pdb_file_name.split('.')[0]
    numbers = str.split('_')
    return {'from': int(numbers[0]),'till': int(numbers[1])}

from Bio.PDB import StructureBuilder
from Bio.PDB import *

def copy_residues(pdbs_directory, chainID, modelID = 0):
    #debug check_array = [False for i in range(0,3000)]
    file_list = get_pdb_files(pdbs_directory)
    parser_pdb = PDBParser()
    structure_builder = StructureBuilder.StructureBuilder()
    structure_builder.init_structure("struct")
    structure_builder.init_model(0)
    structure_builder.init_chain("A")
    connected_structure = structure_builder.get_structure()
    markingList = __initializeMarkingList__()
    for file in file_list:
        range_ = get_range_to_copy(file)
        predicted_part_structure = parser_pdb.get_structure("predicted_part", pdbs_directory + "/" + file)
        for res in predicted_part_structure[modelID][chainID]:
            if res.id[1] >= range_['from'] and res.id[1] <= range_['till']:
                if not __isAlredyResidueInFinalStructure__(res.id[1], markingList, (range_['till'] - range_['from'] + 1)):
                    connected_structure[modelID][chainID].insert(res.id[1], res)
                    #print markingList
    #debug for res in connected_structure[modelID][chainID]:
        #debug check_array[res.id[1]] = True
    #debug for i in range(0,2999):
        #debug if not check_array[i]:
            #debug print i
    io = PDBIO()
    io.set_structure(connected_structure)
    io.save("predicted_structure.pdb")

def __initializeMarkingList__():
    list = [None] * 500  # Todo: should not be constant
    return list

def __isAlredyResidueInFinalStructure__(residueId, markExisting, predictedSubseqLength):
    if markExisting[residueId] is None:
        markExisting[residueId] = predictedSubseqLength
        return False
    elif markExisting[residueId] > predictedSubseqLength:
        markExisting[residueId] = predictedSubseqLength
        return True  # treba prepracovat
    else:
        return True
    

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--chainID')
parser.add_argument('-d', '--dir')
args = vars(parser.parse_args())
chainID = args["chainID"]
dir = args["dir"]
print "*******************************************"

import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



with cd(dir):
    print os.getcwd()
    copy_residues("predicted_parts", chainID)