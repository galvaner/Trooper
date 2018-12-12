import Bio.PDB
from Bio.PDB import *

chainId = 'A'
newChainID = 'B'
pdbFile = "./input.pdb"
newStructureName = "renamedChain.pdb"


pdb_parser = PDBParser()
structure = pdb_parser.get_structure("fullStructure", pdbFile)
if structure[0][chainId].id != newChainID:
  structure[0][chainId].id = newChainID
io_1 = PDBIO()
io_1.set_structure(structure)
io_1.save(newStructureName)#, select=SelectResidues(resSeq, newChainID))

