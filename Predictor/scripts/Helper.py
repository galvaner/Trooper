from Bio import SeqIO
from Bio.PDB import *
import os


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


# expects fastaName like XXXX_A.anything
# seondary structure is on second line
def LoadSecondaryStructure(fastaName):
    rawName = GetFastaNAmeFromFileName(fastaName)
    with open(__makeSecStrName__(rawName), 'r') as file:
        sequence_raw = file.readline()
        sequence = ""
        for res in sequence_raw:
            if res.upper() in ('A', 'G', 'C', 'U'):
                sequence += res.upper()
        secondaryStructure_raw = file.readline()
        secondaryStructure = ""
        for res in secondaryStructure_raw:
            if res.upper() in ('.', '(', ')', '[', ']'):
                secondaryStructure += res.upper()
    return {'sequence': sequence, 'sec_str': secondaryStructure}


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
def is_fast_and_pdb_ordered_correctly(input_pdb, chain_id, fasta, fasta_shift=0):
    if not os.path.isfile(input_pdb):
        print 'ERROR: file' + input_pdb + " does not exists."
        return False
    if not os.path.isfile(fasta):
        print 'ERROR: file' + fasta + " does not exists."
        return False
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('self', input_pdb)
    model = structure[0]
    chain = model[chain_id]
    records = list(SeqIO.parse(fasta, "fasta"))
    fasta_seq = ""
    for r in records[0]:
        fasta_seq = fasta_seq + r
    fasta_seq = "x"+fasta_seq+"x"
    for i in range(0, fasta_shift):
        fasta_seq = "y" + fasta_seq
    b = True
    for res in chain:
        if res.id[1] >= len(fasta_seq):
            b = False
            break
        if res.resname[2] != fasta_seq[res.id[1]]:
            b = False
    return b


# select only ATOM lines with correct chain_id to template
def select_relevant_chain_from_template_pdb(pdb, chain_id):
    import re
    with open("template.pdb", "wb") as modified_pdb:
        with open(pdb, 'r') as original_pdb:
            for line in original_pdb.readlines():
                if re.match('ATOM.................'+chain_id, line):
                    modified_pdb.write(line)



def load_fasta(fasta_file):
    from Bio import SeqIO
    FastaFile = open("../fastas/" + fasta_file.upper() + ".fasta", 'rU')
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
    FastaFile.close()
    return seq



# Inserts new inside original at pos.
def insert_char_to_str(original, new, pos):
    return original[:pos] + new + original[pos:]


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
