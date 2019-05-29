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
def Merge(target_name, template_name, fill_gap_template_ids, template_connector_residue_ids, target_connector_residue_ids, template_chain_id, target_chain_id, create_separate_files_for_gaps, gap_name):
  '''Initialy made to superposition another templates to gaps in primary template, but Rosetta can not cope with this. Currently creating new files fit another templates'''
  temp = "temp.pdb"
  # get only the first in from tuple (its residue id in template structure)
  fill_gap_template_ids_original_ids = [i[0] for i in fill_gap_template_ids]
  if create_separate_files_for_gaps:
      createPdbFileForGap(gap_name, template_name, template_chain_id, fill_gap_template_ids_original_ids, target_chain_id, fill_gap_template_ids)
  else:
      createTempSampleStruct(temp, template_name, template_chain_id, fill_gap_template_ids_original_ids + template_connector_residue_ids, target_chain_id)

      pdb_parser = Bio.PDB.PDBParser(QUIET = True)

      # Get the structures
      ref_structure = pdb_parser.get_structure("reference", target_name)
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
          if ref_res.get_id()[1] in target_connector_residue_ids:
            # Append CA atom to list
            ref_atoms.append(ref_res['P'])

      # Do the same for the sample structure
      for sample_chain in sample_model:
        for sample_res in sample_chain:
          if sample_res.get_id()[1] in template_connector_residue_ids:
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
      ConnectAndSaveAlignedTemplateToOriginalTemplate(target_name, sample_structure, ref_structure, fill_gap_template_ids, target_chain_id)


def ConnectAndSaveAlignedTemplateToOriginalTemplate(toStruct, sample_structure, ref_structure, atomsFromSampleStrToSave, commonChainID):
    for sampleRes in sample_structure[0][commonChainID]:
        new_residue_id = CheckIsNumberInListAndReturnNewId(sampleRes.id[1], atomsFromSampleStrToSave)
        if new_residue_id != -1:
            # renumber added ids
            residue_to_copy = sampleRes
            if new_residue_id != residue_to_copy.id[1]:
                residue_to_copy.detach_parent()
                residue_to_copy.id = (' ', new_residue_id, ' ')
            ref_structure[0][commonChainID].add(residue_to_copy)

    #remove residues like "RIA", because rosetta can not work with them
    for res in ref_structure[0][commonChainID]:
        if (res.id[0] != '' and res.id[0] != ' ') or (res.id[2] != '' and res.id[0] != ' '):
            ref_structure[0][commonChainID].detach_child(res.id)

    # sort structure
    ref_structure[0][commonChainID].child_list.sort(key=lambda x: x.id[1])

    ## gaps in structure has to be longer than 1 residue (else rosetta wont find any fragment in its library)
    #last_res_id = 1
    #for res in ref_structure[0][commonChainID]:
    #    if res.id[1] == last_res_id + 2:
    #        ref_structure[0][commonChainID].detach_child(res.id)
    #    else:
    #        last_res_id = res.id[1]

    io = Bio.PDB.PDBIO()
    io.set_structure(ref_structure)
    io.save(toStruct)


def CheckIsNumberInListAndReturnNewId(number, list):
    for listElement in list:
        if listElement[0] == number:
            return listElement[1]
    return -1

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

def createPdbFileForGap(new_structure, pdbFile, chainId, resSeq, newChainID, atomsFromSampleStrToSave):
  pdb_parser = PDBParser()
  temp = 'new_str_for_gap.pdb'
  structure = pdb_parser.get_structure("fullStructure", pdbFile)
  if structure[0][chainId].id != newChainID:
    structure[0][chainId].id = newChainID
  io_1 = PDBIO()
  io_1.set_structure(structure)
  io_1.save(temp, select=SelectResidues(resSeq, newChainID))

  sample_structure = pdb_parser.get_structure("sample", temp)
  ref_structure = pdb_parser.get_structure("sample", temp)
  # detach all children, so I can build the structure from scratch
  all_ids_to_detach = []
  for res in ref_structure[0][newChainID]:
      all_ids_to_detach.append(res.id)
  for id in all_ids_to_detach:
      ref_structure[0][newChainID].detach_child(id)

  for sampleRes in sample_structure[0][newChainID]:
      new_residue_id = CheckIsNumberInListAndReturnNewId(sampleRes.id[1], atomsFromSampleStrToSave)
      if new_residue_id != -1:
          # renumber added ids
          residue_to_copy = sampleRes
          if new_residue_id != residue_to_copy.id[1]:
              residue_to_copy.detach_parent()
              residue_to_copy.id = (' ', new_residue_id, ' ')
          ref_structure[0][newChainID].add(residue_to_copy)

  # remove residues like "RIA", because rosetta can not work with them
  for res in ref_structure[0][newChainID]:
      if (res.id[0] != '' and res.id[0] != ' ') or (res.id[2] != '' and res.id[0] != ' '):
          ref_structure[0][newChainID].detach_child(res.id)

  # sort structure
  ref_structure[0][newChainID].child_list.sort(key=lambda x: x.id[1])

  io = Bio.PDB.PDBIO()

  # remove all redundant models
  model_ids = []
  zero_id = False
  for model in ref_structure:
    if model.id != 0:
        model_ids.append(model.id)
    else:
        zero_id = True
  for id in model_ids:
    if zero_id != False:
        ref_structure.detach_child(id)
    else:
        zero_id = True

  io.set_structure(ref_structure)
  io.save(new_structure)


class SelectResidues(Select):
    def __init__(self, residues_id, chain_id):
        self.residues_id = residues_id
        self.chain_id = chain_id

    def accept_residue(self, residue):
        if ((residue.id[1] in self.residues_id) and (residue.parent.id == self.chain_id)):
            return 1
        else:
            return 0
