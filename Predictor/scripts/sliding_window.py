__author__ = 'Rasto'

import re
from Bio.PDB import *

#this method expect data in alignment like I grep/sed-ed them or like it was given by cp_predict
def prepare_alignment_data(alignment_file_name, target_name, template_name):
    alignment_file = open(name = alignment_file_name, mode = 'r')
    template_seq = 'x'
    target_seq = 'y'
    for line in alignment_file:
        matched_temp = re.match('^'+template_name+' *([A-Za-z-]*) *', line)
        matched_trgt = re.match('^'+target_name+' *([A-Za-z-]*) *', line)
        if matched_temp:
            template_seq += matched_temp.group(1)
        if matched_trgt:
            target_seq += matched_trgt.group(1)
    alignment_file.close()
    return data(template_seq, target_seq, None, None)

class data(object):
    def __init__(self, template_fasta_alignment, target_fasta_alignment, window_size, match_lower_bound ):
        self.template_fasta_alignment = template_fasta_alignment
        self.target_fasta_alignment = target_fasta_alignment
        self.window_size = window_size
        self.match_lower_bound = match_lower_bound

def add_params_to_half_populated_data_object(data, window_size, match_lower_bound):
    data.window_size = window_size
    data.match_lower_bound = match_lower_bound
    return data

def compute_conserved_regions(data): #data has to be type class data
    length = len(data.target_fasta_alignment)
    length_template = len(data.template_fasta_alignment)
    if length != length_template:
        print "sliding_window.py: WARNING Difference in alignment lengths found!"
    conserved_parts =  [False for k in range(0,length)]
    if data.window_size % 2 == 1: #need odd number
        data.window_size = data.window_size + 1
    window_half_length = data.window_size // 2
    for i in range(window_half_length, length-window_half_length):#sliding window
        match = 0
        for j in range(i-window_half_length, i + window_half_length):
            if data.target_fasta_alignment[j] == data.template_fasta_alignment[j]:
                match = match + 1
        if match >= data.match_lower_bound:
            conserved_parts[i] = True
    template_conserved_indicator = []
    for i in range(0, length): # create mapping of bool object(indication of conserved residues) and template fasta (sth like delete "-" characters)
        if data.template_fasta_alignment[i] != "-":
            template_conserved_indicator.append(conserved_parts[i])
    #debug print "Length without \"-\":" + str(len(template_conserved_indicator)-1)
    #debug print len(template_conserved_indicator)
    return template_conserved_indicator

class SelectResidues(Select):
    def __init__(self, indicator_array, chain_id): #in indicator array bool means that residue should be preserved
        self.indicator_array = indicator_array
        self.chain_id = chain_id

    def accept_residue(self, residue):
        if residue.get_parent().get_id() != self.chain_id:
            return 0
        if self.indicator_array[residue.id[1]]:
            return 1
        else:
            return 0

def pick_conserved_res_from_template_pdb(template_pdb, template_conserved_indicator, chainID, path ): # can be added before / after edit_conserved_regions
    start = -1
    for i in range(1,len(template_conserved_indicator)):
        if start != -1 and not template_conserved_indicator[i]:
            continue
        if not template_conserved_indicator[i]:
            start = i
            continue
        if start != -1:
            print "Unconserved:" + str(start) + "-" + str(i-1)
            start = -1
    print "Unconserved:" + str(start) + "-" + str(i)

    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure('template_pdb', template_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save( path + "conserved_regions.pdb", select=SelectResidues(template_conserved_indicator, chainID))

def run_sliding_window(alignment_file_name, target_name, template_name, window_size, match_lower_bound, template_pdb_name, chainID, path):
    print "SLIDING_WINDOW: running sliding window algorithm..."
    data = prepare_alignment_data(alignment_file_name, target_name, template_name)
    data = add_params_to_half_populated_data_object(data, window_size, match_lower_bound)
    indicator_array =  compute_conserved_regions(data)
    pick_conserved_res_from_template_pdb(template_pdb_name, indicator_array, chainID, path)
    print "SLIDING_WINDOW: running sliding window algorithm finished."


