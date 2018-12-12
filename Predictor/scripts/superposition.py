# The MIT License
#
# Copyright (c) 2010-2016 Anders S. Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Merge method in this script is inspired by code from https://gist.github.com/andersx/6354971


import Bio.PDB
from Bio.PDB import *


# Start the parser
def Merge(toStruct, fromStruct, inGapRes, connectorsRes, chainIdFrom, originalTemplateChainId):
  temp = "temp.pdb"
  atoms_to_be_aligned = connectorsRes
  createTempSampleStruct(temp, fromStruct, chainIdFrom, inGapRes+connectorsRes, originalTemplateChainId)
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)

  # Get the structures
  ref_structure = pdb_parser.get_structure("reference", toStruct)
  sample_structure = pdb_parser.get_structure("sample", temp)

  # Use the first model in the pdb-files for alignment
  # Change the number 0 if you want to align to another structure
  ref_model = ref_structure[0]
  sample_model = sample_structure[0]

  # Make a list of the atoms (in the structures) you wish to align.
  # In this case we use CA atoms whose index is in the specified range
  ref_atoms = []
  sample_atoms = []

  # Iterate of all chains in the model in order to find all residues
  for ref_chain in ref_model:
    # Iterate of all residues in each model in order to find proper atoms
    for ref_res in ref_chain:
      # Check if residue number ( .get_id() ) is in the list
      if ref_res.get_id()[1] in atoms_to_be_aligned:
        # Append CA atom to list
        ref_atoms.append(ref_res['P'])

  # Do the same for the sample structure
  for sample_chain in sample_model:
    for sample_res in sample_chain:
      if sample_res.get_id()[1] in atoms_to_be_aligned:
        sample_atoms.append(sample_res['P'])

  CheckIfAtomsArePaired(ref_atoms, sample_atoms)

  # Now we initiate the superimposer:
  super_imposer = Bio.PDB.Superimposer()
  super_imposer.set_atoms(ref_atoms, sample_atoms)
  super_imposer.apply(sample_model.get_atoms())

  # Print RMSD:
  print super_imposer.rms

  # Save the aligned version of 1UBQ.pdb
  io = Bio.PDB.PDBIO()
  io.set_structure(sample_structure)
  #io.set_structure(ref_structure)
  #io.save(toStruct + ".rotated.out")
  ConnectAndSaveAlignedTemplateToOriginalTemplate(toStruct, sample_structure, ref_structure, inGapRes, originalTemplateChainId)

def ConnectAndSaveAlignedTemplateToOriginalTemplate(toStruct, sample_structure, ref_structure, atomsFromSampleStrToSave, commonChainID):
    for sampleRes in sample_structure[0][commonChainID]:
        if IsNumberInList(sampleRes.id[1], atomsFromSampleStrToSave):
            ref_structure[0][commonChainID].add(sampleRes)
    ref_structure[0][commonChainID].child_list.sort(key=lambda x: x.id[1])
    io = Bio.PDB.PDBIO()
    io.set_structure(ref_structure)
    io.save(toStruct)

def IsNumberInList(number, list):
    for listElement in list:
        if  listElement == number:
            return True
    return False

def CheckIfAtomsArePaired(ref_atoms, sample_atoms): #, firstGapRes, lastGapRes):
    ref_atoms.sort(key=lambda x: x.parent.id[1], reverse=True)
    sample_atoms.sort(key=lambda x: x.parent.id[1], reverse=True)

def createTempSampleStruct(newStructureName, pdbFile, chainId, resSeq, newChainID):
  pdb_parser = PDBParser()
  structure = pdb_parser.get_structure("fullStructure", pdbFile)
  if structure[0][chainId].id != newChainID:
    structure[0][chainId].id = newChainID
  io_1 = PDBIO()
  io_1.set_structure(structure)
  io_1.save(newStructureName, select=SelectResidues(resSeq, newChainID))


class SelectResidues(Select):
    def __init__(self, residues_id, chain_id):
        self.residues_id = residues_id
        self.chain_id = chain_id

    def accept_residue(self, residue):
        if ((residue.id[1] in self.residues_id) and (residue.parent.id == self.chain_id)):
            return 1
        else:
            return 0
