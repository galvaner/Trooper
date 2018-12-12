__author__ = 'Rasto'

import os
from subprocess import call
from Bio.PDB import *
import align
import sliding_window
import edit_conserved_regions
import cut_spheres
import prepare_rosetta_script
import config
import TemplateSelector
import argparse
import Helper
import PredictSecondaryStructure
import UseSecondTemplate

#cashe = {}
run_on_windows = False



def connect_scripts(file_name):
    ttpairs = open(file_name)
    # parse file "pairs.txt"
    parsed_tt = []
    used_target = {}
    used_template = {}
    for line in ttpairs:
        parsed_tt.append(process_pair(line))
    for line in parsed_tt:
        Helper.modify_pdb(line['template_pdb'], line['template_chainID'])
        global bad_cache
        if line['template_fasta'] in bad_cache:
            continue
        if(Helper.check_order_of_fasta_and_pdb('./template.pdb', line['template_chainID'], line['template_fasta'])):
            #global cashe
            #cashe[line['template_fasta']] = True
            if not line['target'] in used_target and not line['template'] in used_target and not line['template'] in used_template and not line['target'] in used_template:
                print("**********template***********" + line['template'] + "***********target**********" + line['target'])
                create_directories(line['target'], line['template'])
                prepare_prediction(line)
                ######Vypnute filtrovanie rovnakych used_target[line['target']] = True
                ######Vypnute filtrovanie rovnakych used_template[line['template']] = True
                #return
            #else:
                #print line['target'] + " alredy used as target!"
        else:
            #if False and try_adding_ens(line['template_pdb'], line['template_chainID'], line['template_fasta']) :
            #    #create_directories(line['target'], line['template'])
            #    #prepare_prediction(line)
            #    print line['template_fasta'] + " saved by n adding"
            #else:
            #    bad_cache[line["template_fasta"]] = False
            bad_cache[line["template_fasta"]] = False
            print(line['template'] + " is not suitable as template for prediction without renumbering!")

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


def shell_script_call():
    import subprocess
    if run_on_windows:
        subprocess.call(['C:\\Cygwin\\bin\\sh', './script.sh'], shell=True)
    else:
        subprocess.call(['chmod', '+x' , 'script.sh'])
        subprocess.call(['./script.sh'])
def create_shell_script(target, template, pdb, chainID, dir_):
    file = open("script.sh", 'wb')
    file.write('#!/bin/bash \n grep "^SEQUENCE" "'+ dir_ +'emboss.aln" | sed -e "s/^SEQUENCE/'+target+'/;N" | sed -e "s/^SEQUENCE/'+template+'/" | sed -e "s/ [[:digit:]]*/ /g" > '+ dir_ + 'tta.aln')
    file.write('\ngrep "^ATOM.................'+chainID + '" ' + pdb + ' > ' + dir_ + 'template.pdb' )
    file.close()

def prepare_prediction(parsed_tt):
    try:
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
        if config.predictSecondaryStructure:
            secStrObject = PredictSecondaryStructure.RunSecStrPRediction(parsed_tt["target"], parsed_tt["template"], dir, run_on_windows) #  predict secondary structure
        prepare_rosetta_script.prepare_rosetta_script(config.configNumberOFStructuresGeneratedByFARFAR, parsed_tt["template_chainID"], length, parsed_tt['dir'], secStrObject)
    except:
        print "prepare_prediction: ERROR catched!"


def __prepare_prediction_for_single_target__(target):
    suitable_templates = TemplateSelector.SelectTemplate(target)
    if len(suitable_templates) == 0:
        raise "No template found for " + target + " sequence"
    target_fasta = target + ".fasta"
    target_template_pairs = [target_fasta + " " + template + ".fasta" for template in suitable_templates]
    prepared_target_template_pairs = [process_pair(pair) for pair in target_template_pairs]
    for i in xrange(len(prepared_target_template_pairs)):
        parsed_tt = prepared_target_template_pairs[i]
        print("**********template***********" + line['template'] + "***********target**********" + line['target'])
        create_directories(line['target'], line['template'])
        # like prepare prediction but with possibility of multiple templates
        if config.useMultipleTemplates:
            try:
                align.compute_emboss_alignment(parsed_tt['target_fasta'], parsed_tt['template_fasta'],
                                               parsed_tt['dir'] + "emboss.aln")
                create_shell_script(parsed_tt["target"], parsed_tt["template"], parsed_tt["template_pdb"],
                                    parsed_tt["template_chainID"], parsed_tt['dir'])
                shell_script_call()
                sliding_window.run_sliding_window(parsed_tt['dir'] + "tta.aln", parsed_tt["target"],
                                                  parsed_tt["template"], config.configSlidingWindowSize,
                                                  config.configSlidingWindowMinCoverage,
                                                  parsed_tt['dir'] + "template.pdb", parsed_tt['template_chainID'],
                                                  parsed_tt['dir'])
                edit_conserved_regions.edit_conserved_regions(parsed_tt['template'], parsed_tt['target'],
                                                              parsed_tt['dir'] + "tta.aln",
                                                              parsed_tt['dir'] + "conserved_regions.pdb",
                                                              parsed_tt['dir'] + "conserved_edited_regions.pdb",
                                                              parsed_tt['dir'] + "target.fasta",
                                                              parsed_tt["template_chainID"])
                length = get_target_length(parsed_tt['target_fasta'])
                print("Target length = " + str(length))
                if config.useMultipleTemplates:
                    for j in xrange(len(suitable_templates)):
                        if j == i:
                            continue
                        UseSecondTemplate.SecondTemplateUsage(parsed_tt['dir'] + "conserved_edited_regions.pdb"
                                                              , parsed_tt['target']
                                                              , suitable_templates[j]
                                                              , parsed_tt["template_chainID"])
                cut_spheres.cut_spheres(parsed_tt['dir'] + "conserved_edited_regions.pdb", length,
                                        parsed_tt['dir'] + 'cut_spheres/', parsed_tt["template_chainID"], 2.7, 5)

                if config.predictSecondaryStructure:
                    secStrObject = PredictSecondaryStructure.RunSecStrPRediction(parsed_tt["target"],
                                                                                 parsed_tt["template"], dir,
                                                                                 run_on_windows)  # predict secondary structure

                prepare_rosetta_script.prepare_rosetta_script(config.configNumberOFStructuresGeneratedByFARFAR,
                                                              parsed_tt["template_chainID"], length, parsed_tt['dir'],
                                                              secStrObject)
            except:
                print "prepare_prediction: ERROR catched!"


if config.predictListOfTargets:
    if config.singleTargetPredictionFromCmd:
        parser = argparse.ArgumentParser()
        parser.add_argument('-t', '--target')
        parser.add_argument('-w', '--windows')
        args = vars(parser.parse_args())
        target = args["target"]
        run_on_windows = bool(args["windows"])
        __prepare_prediction_for_single_target__(target)
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--file')
        parser.add_argument('-w', '--windows')
        args = vars(parser.parse_args())
        file = args["file"]
        run_on_windows = bool(args["windows"])
        #todo: multi target prediction from file
        targets = open(file)
        for line in targets:
            __prepare_prediction_for_single_target__(line)
        file.close()
else:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file')
    parser.add_argument('-w', '--windows')
    args = vars(parser.parse_args())
    file = args["file"]
    run_on_windows = bool(args["windows"])
    connect_scripts(file)
#create_directories("target", "template")
#create_shell_script('x','y')
#shell_script_call()