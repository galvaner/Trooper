__author__ = 'Rasto'

# script used to automatically predict RNA structures by ModeRNA.
# ModeRNA has to be installed in the system - I used it only with Windows operating system
# expect my classic placement of fastas and pdbs (but in the same directory as this script)
# using modeRNA.py script

import os
from shutil import copyfile
import modeRNA
import re
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def modify_pdb(in_pdb_file, out_pdb_file, chainID):
    with open(in_pdb_file,'r') as origin_file:
            updated_pdb = open(out_pdb_file, 'w')
            for line in origin_file:
                pattern = r'^ATOM...............[AGCU].'+chainID+'.*'
                line_ = re.findall(pattern, line)
                if len(line_) != 0:
                    updated_pdb.write(line_[0]+'\n')
            updated_pdb.close()



newpath = './prediction'
if not os.path.exists(newpath):
    os.makedirs(newpath)

for file in os.listdir("./pairs"):
    if not os.path.exists("./prediction/"+file):
        os.makedirs("./prediction/"+file)
    else:
        print "Directory " + file + " alredy exists."
        continue
    pairs = open("./pairs/"+file, 'r')
    os.chdir("./prediction/"+file)
    for line in pairs:
        target = line[:4]
        target_chain = line[5:6]
        template = line[13:17]
        template_chain = line[18:19]
        target_fasta = target+'_'+target_chain
        prediction_dir = target+'_'+target_chain+'_'+template+'_'+template_chain
        print "Working on " + prediction_dir + " model"
        os.makedirs(prediction_dir)
        os.chdir(prediction_dir)
        copyfile('../../../pdbs/' + template.lower() + '.pdb', './template_.pdb')
        modify_pdb("./template_.pdb", "./template.pdb", template_chain)
        copyfile('../../../pdbs/' + target.lower() + '.pdb', './target_.pdb')
        modify_pdb("./target_.pdb", "./target.pdb", target_chain)
        copyfile('../../../fastas/' + target_fasta + '.fasta', './target.fasta')
        struct = modeRNA.molecule_structure("template.pdb", template_chain)
        struct.create_fasta_from_pdb("template.fasta", "TEMPLATE")
        modeRNA.create_alignment("target.fasta", "template.fasta", "tta_.aln")
        modeRNA.modify_alignment("tta_.aln")
        try:
            modeRNA.build_model(template_chain)
            print bcolors.OKGREEN + "Model "+ prediction_dir  +" was successfully created." + bcolors.ENDC
        except:
            e = sys.exc_info()[0]
            print(str(e))
        os.chdir('..')
    os.chdir("../..")
    pairs.close()

