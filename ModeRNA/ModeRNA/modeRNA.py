__author__ = 'Rasto'

# script used by runModernaPredictionBasedOnPairs.py 
# ModeRNA has to be installed in the system - I used it only with Windows operating system

from Bio.PDB import *
from Bio.Emboss.Applications import *
from moderna import *

class molecule_structure:
    def __init__(self, pdb, chain):
        self.pdb_file=pdb
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure('self', pdb)
        self.pdb_model = structure[0]
        self.pdb_chain = self.pdb_model[chain]
        self.chainID = chain

    def create_fasta_from_pdb(self, output_file, fasta_name):
        file = open(output_file, 'w')
        file.write(">"+ fasta_name + "\n")
        self.fasta_from_pdb = ""
        for residue in self.pdb_chain:
            file.write(residue.get_resname())
            self.fasta_from_pdb+=str(residue.get_resname()[1])
        file.close()

def create_alignment(input_file_target, input_file_template, output_file):
    print "Creating alignment..."
    needle_cli = NeedleCommandline(asequence=input_file_target, \
                           bsequence=input_file_template, \
                           gapopen=10, \
                           gapextend=0.5, \
                           outfile=output_file, \
                           )
    needle_cli()

#expects that target.fasta starts with >SEQUENCE and template with >TEMPLATE
def modify_alignment(alignment_file):
    print "Modifying alignment..."
    original_alignment = open(alignment_file)
    template = ">TEMPLATE\n"
    target = ">TARGET\n"
    for line in iter(original_alignment):
        if line.startswith("TEMPLATE"):
            template += line[21:-8]
        if line.startswith("SEQUENCE"):
            target += line[21:-8]
    output = open('tta.aln', 'w')
    output.write(str.upper(target))
    output.write('\n')
    output.write(str.upper(template))
    output.close()
    original_alignment.close()
def build_model(chainID):
    print "Building model..."
    t = load_template('template.pdb', chainID)
    a = load_alignment('tta.aln')
    m = create_model(t,a, chainID)
    m.write_pdb_file('model.pdb')


#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('-chainID', '--chainID')
#args = vars(parser.parse_args())
#
#
#chainID = args["chainID"]
#struct = molecule_structure("template.pdb", chainID)
#struct.create_fasta_from_pdb("template.fasta", "TEMPLATE")
#create_alignment("target.fasta", "template.fasta", "tta_.aln")
#modify_alignment("tta_.aln")
#build_model(chainID)