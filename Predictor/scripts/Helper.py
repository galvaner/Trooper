# XXXX_anything, I want firts 4 characters in lowercase
def TrimPDBName(fastaName):
    return fastaName[:4].lower()


# XXXX_A_anything, I want A
def GetChainID(fastaName):
    return fastaName[5:6]


def GetFastaNAmeFromFileName(fileName):
    return fileName[:6]

CONST_SEC_STR_FOLDER = '../secondary_structures/'

def __makeSecStrName__(rawName):
    return CONST_SEC_STR_FOLDER + rawName + ".secstr"

#expects fastaName like XXXX_A.anything
def LoadSecondaryStructureToString(fastaName):
    rawName = GetFastaNAmeFromFileName(fastaName)
    with open(__makeSecStrName__(rawName), 'r') as file:
        firstLine = file.next()
    return firstLine


def ListOfPairsToFiles(zeroIndexFileName, listOfPairs):
    fileZero = open(zeroIndexFileName, 'w')
    fileZero.write(">SEQ:\n")
    for pair in listOfPairs:
        fileZero.write(pair[0])
    fileZero.write('\n')
    for pair in listOfPairs:
        fileZero.write(pair[1])
    fileZero.close()

# compare order of fasta file and corresponding pdb file
def check_order_of_fasta_and_pdb(input_pdb, chainID, fasta, fasta_shift = 0):
    #global cashe
    #if fasta in cashe:
    #    return True
    from Bio import SeqIO
    from Bio.PDB import *
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('self', input_pdb)
    model = structure[0]
    chain = model[chainID]
    records = list(SeqIO.parse(fasta, "fasta"))
    fasta_seq = ""
    for r in records[0]:
        fasta_seq = fasta_seq + r
    fasta_seq = "x"+fasta_seq+"x"
    for i in range(0,fasta_shift):
        fasta_seq = "y" + fasta_seq
    b = True
    for res in chain:
        if res.id[1] >= len(fasta_seq):
            b = False
            break
        if res.resname[2] != fasta_seq[res.id[1]]:
            b = False
    return b

def modify_pdb(pdb, chainID, dir_ = './', run_on_winows = False):
    file = open("script_modify_pdb.sh", 'wb')
    file.write('#!/bin/bash')
    file.write('\ngrep "^ATOM.................'+chainID + '" ' + pdb + ' > ' + dir_ + 'template.pdb' )
    file.close()
    import subprocess
    if run_on_winows:
        subprocess.call(['C:\\Cygwin\\bin\\sh', './script_modify_pdb.sh'], shell=True)
    else:
        subprocess.call(['chmod', '+x' , 'script_modify_pdb.sh'])
        subprocess.call(['./script_modify_pdb.sh'])

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'