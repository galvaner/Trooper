__author__ = 'Rasto'

import os
from subprocess import call
from Bio.PDB import *
import align
import sliding_window
import edit_conserved_regions
import cut_spheres
import prepare_rosetta_script
import PredictSecondaryStructure
import config

cashe = {}

def modify_pdb(pdb, chainID, dir_ = './'):
    file = open("script_modify_pdb.sh", 'wb')
    file.write('#!/bin/bash')
    file.write('\ngrep "^ATOM.................'+chainID + '" ' + pdb + ' > ' + dir_ + 'template.pdb' )
    file.close()
    import subprocess
    #####subprocess.call(['C:\\Cygwin\\bin\\sh', './script_modify_pdb.sh'], shell=True)
    subprocess.call(['chmod', '+x' , 'script_modify_pdb.sh'])
    subprocess.call(['./script_modify_pdb.sh'])


def connect_scripts(max_count_to_predict, file_name):
    ttpairs = open(file_name)
    # parse file "pairs.txt"
    parsed_tt = []
    used_target = {}
    used_template = {}
    for line in ttpairs:
        parsed_tt.append(process_pair(line))
    for line in parsed_tt:
        modify_pdb(line['template_pdb'], line['template_chainID'])
        global bad_cache
        if line['template_fasta'] in bad_cache:
            continue
        if(check_if_renumber_is_needed('./template.pdb', line['template_chainID'], line['template_fasta'])):
            global cashe
            cashe[line['template_fasta']] = True
            if max_count_to_predict <= 0:
                print "connect_scripts.py: Max prediction limit was reached!"
                break
            if not line['target'] in used_target and not line['template'] in used_target and not line['template'] in used_template and not line['target'] in used_template:
                max_count_to_predict = max_count_to_predict - 1
                print("**********template***********" + line['template'] + "***********target**********" + line['target'])
                create_directories(line['target'], line['template'])
                prepare_prediction(line)
                ######Vypnute filtrovanie rovnakych used_target[line['target']] = True
                ######Vypnute filtrovanie rovnakych used_template[line['template']] = True
                #return
            #else:
                #print line['target'] + " alredy used as target!"
        else:
            if False and try_adding_ens(line['template_pdb'], line['template_chainID'], line['template_fasta']) :
                #create_directories(line['target'], line['template'])
                #prepare_prediction(line)
                print line['template_fasta'] + " saved by n adding"
            else:
                bad_cache[line["template_fasta"]] = False
            print(line['template'] + " is not suitable as template for prediction without renumbering!")
    print 10000-max_count_to_predict

def process_pair(pair):
    target = pair.split(' ')[0].split('.')[0]
    template = pair.split(' ')[1].split('.')[0]
    template_chainID = template.split('_')[1]
    template_pdb = "../pdbs/" + template.split('_')[0].lower() + ".pdb"
    template_fasta = "../fastas/" + template + ".fasta"
    target_fasta = "../fastas/" + target + ".fasta"
    dir = "../prediction/" + target + "_" + template + "/files/"
    return {"target": target, "template":template, "template_chainID":template_chainID, "template_pdb":template_pdb, "target_fasta":target_fasta, "template_fasta":template_fasta, "dir":dir }

def create_directories(target,  template):
        directories = []
        directories.append("../prediction/" + target + "_" + template + "/files/cut_spheres")
        directories.append("../prediction/" + target + "_" + template + "/files/prepare_rosetta_scripts")
        directories.append("../prediction/" + target + "_" + template + "/results/predicted_parts")
        for dir in directories:
            if not os.path.exists(dir):
                os.makedirs(dir)
            else:
                print "connect_scripts.py: WARNING: directory " + dir + " alredy exists"


def get_target_length(target_fasta):
    from Bio import SeqIO
    records = list(SeqIO.parse(target_fasta, "fasta"))
    fasta_seq = ""
    for r in records[0]:
        fasta_seq = fasta_seq + r
    return len(fasta_seq)

bad_cache = {}
def try_adding_ens(input_pdb, chainID, fasta):
    from Bio import SeqIO
    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('self', input_pdb)
    model = structure[0]
    chain = model[chainID]
    records = list(SeqIO.parse(fasta, "fasta"))
    fasta_seq = ""
    for r in records[0]:
        fasta_seq = fasta_seq + r
    fasta_seq = "x"+fasta_seq+"x"
    fasta_shift = 0
    while fasta_shift <= 50:
        fasta_shift = fasta_shift + 1
        fasta_seq = "y" + fasta_seq
        b = True
        for res in chain:
            if res.id[1] >= len(fasta_seq):
                b = False
                break
            if res.resname[2] != fasta_seq[res.id[1]]:
                b = False
                break
        if b:
            #print "+++++good++++++++++++++++" + fasta + "+++++++++++++++++++++++++"
            return True
    return False
    #print "--------wrong-----------" + fasta + "------------------------"


def check_if_renumber_is_needed(input_pdb, chainID, fasta, fasta_shift = 0):
    global cashe
    if fasta in cashe:
        return True
    from Bio import SeqIO
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

def shell_script_call():
    import subprocess
    ###subprocess.call(['C:\\Cygwin\\bin\\sh', './script.sh'], shell=True)
    subprocess.call(['chmod', '+x' , 'script.sh'])
    subprocess.call(['./script.sh'])
def create_shell_script(target, template, pdb, chainID, dir_):
    import io
    dir = target + "_" + template
    file = open("script.sh", 'wb')
    file.write('#!/bin/bash \n grep "^SEQUENCE" "'+ dir_ +'emboss.aln" | sed -e "s/^SEQUENCE/'+target+'/;N" | sed -e "s/^SEQUENCE/'+template+'/" | sed -e "s/ [[:digit:]]*/ /g" > '+ dir_ + 'tta.aln')
    file.write('\ngrep "^ATOM.................'+chainID + '" ' + pdb + ' > ' + dir_ + 'template.pdb' )
    file.close()

def prepare_prediction(parsed_tt):
    align.compute_emboss_alignment(parsed_tt['target_fasta'], parsed_tt['template_fasta'], parsed_tt['dir']+"emboss.aln")
    create_shell_script(parsed_tt["target"], parsed_tt["template"], parsed_tt["template_pdb"], parsed_tt["template_chainID"], parsed_tt['dir'])
    shell_script_call()
    sliding_window.run_sliding_window(parsed_tt['dir']+"tta.aln", parsed_tt["target"], parsed_tt["template"], config.configSlidingWindowSize, config.configSlidingWindowMinCoverage,
                                      parsed_tt['dir']+"template.pdb", parsed_tt['template_chainID'], parsed_tt['dir'] )
    edit_conserved_regions.edit_conserved_regions(parsed_tt['template'], parsed_tt['target'], parsed_tt['dir']+"tta.aln", parsed_tt['dir']+"conserved_regions.pdb"
                           ,parsed_tt['dir']+"conserved_edited_regions.pdb", parsed_tt['dir']+"target.fasta", parsed_tt["template_chainID"])
    length = get_target_length(parsed_tt['target_fasta'])
    print("Target length = " + str(length))
    cut_spheres.cut_spheres(parsed_tt['dir']+"conserved_edited_regions.pdb", length, parsed_tt['dir']+'cut_spheres/', parsed_tt["template_chainID"], 2.7, 5)
    prepare_rosetta_script.prepare_rosetta_script(config.configNumberOFStructuresGeneratedByFARFAR, parsed_tt["template_chainID"], length, parsed_tt['dir'], parsed_tt["target"], parsed_tt["template"])


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f','--file')
args = vars(parser.parse_args())
file = args["file"]



connect_scripts(900000, file)
#create_directories("target", "template")
#create_shell_script('x','y')
#shell_script_call()