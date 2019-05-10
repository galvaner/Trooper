from align import parse_alignment

__author__ = 'Rasto'

#------------------
# Script description:
#
# Script takes pdb, which contains only conserved parts of sequence
# and targert-template-allignment(tta) and then script deletes residues
# which don't correspond in the tta, in case that in tta stands gap against residue script
# deletes also the first closest neighbor from each side of gap.
# also creates fasta file suitable for rna_denovo_setup.py (lowercase)
# Possible futere TODO:
# Current constant in gap vs. residue delete is 1.
# Try add a dependency on the gap length.
#------------------

from Bio.PDB import *


def handle_gap(residues_to_delete_list, gap_length, from_res): #deletes half of gap length from both sides of gap so if gap "contains" 5 residues 3+5+3 will be deleted
    extra_delete = (gap_length + 1) // 2
    #debugprint "RealDelete: " + str(extra_delete) + "  Gap_length: " + str(gap_length)
    #debugprint str(from_res - extra_delete) + "-" + str(from_res + gap_length + extra_delete)
    for i in range(from_res - extra_delete, from_res + gap_length + extra_delete):
        #debugprint i
        residues_to_delete_list.append(i)


def edit_conserved_regions(template_name, target_name, alignment_file,input_pdb, output_pdb, target_fasta_output,chainID="A", modelID=0 ):
    print "EDIT_CONSERVED_REGIONS: start algorithm"
    alignment_data = prepare_alignment_data(alignment_file, target_name, template_name)
    template_fasta = alignment_data.template_fasta_alignment
    target_fasta = alignment_data.target_fasta_alignment
    renumber = []
    renumber.append([0,0])
    #"pointry" na pohyb vo fastach
    target_pointer = 0
    template_pointer = 0

    #gap length dependancy
    gap_length = 0
    start_res_gap = 0
    end_res_gap = 0
    #--------------------
    residues_to_delete = []
    for i in range(1,len(target_fasta)): #index od nuly!!!!!!
        if template_fasta[i] != '-':
            template_pointer += 1
        if target_fasta[i] != '-':
            target_pointer += 1
        renumber.append([template_pointer,target_pointer])  #map template residues with target residues
        if template_fasta[i] == target_fasta[i]:
            if gap_length > 0:
                handle_gap(residues_to_delete, gap_length, start_res_gap)
                gap_length = 0
            continue
        if template_fasta[i] != target_fasta[i] and template_fasta[i] != '-' and target_fasta[i] != '-':
            if gap_length > 0:
                handle_gap(residues_to_delete, gap_length, start_res_gap)
                gap_length = 0
            residues_to_delete.append(template_pointer)

        if template_fasta[i] == '-' and target_fasta[i] != '-':
            if gap_length == 0:
                start_res_gap = template_pointer
            gap_length = gap_length + 1
        if template_fasta[i] != '-' and target_fasta[i] == '-':
            if gap_length == 0:
                start_res_gap = template_pointer
            gap_length = gap_length + 1

    residues_to_delete = list(set(residues_to_delete)) # I want to get unique values
    for i in residues_to_delete: # vymazanie pripadnych zapornych alebo moc velkych hodnot
        if i < 1 or i > template_pointer:
            residues_to_delete.remove(i)
    create_edited_pdb(residues_to_delete,renumber, input_pdb,output_pdb, chainID=chainID,model=modelID)
    create_lowercase_target_fasta_from_tta(target_fasta, target_fasta_output)
    print "EDIT_CONSERVED_REGIONS: algorithm done"


class SelectResidues(Select):
    def __init__(self, residues_id, chain_id):
        self.residues_id = residues_id
        self.chain_id = chain_id

    def accept_residue(self, residue):
        is_conserved = True
        for i in self.residues_id:
                nm = int(residue.id[1])
                if (i == nm ) | (residue.parent.id != self.chain_id):
                    is_conserved=False
                    break
        if is_conserved:
            return 1
        else:
            return 0

def create_edited_pdb(residues_to_delete, template_target_mapping,input_pdb, output_pdb,chainID, model=0):
    parser_pdb = PDBParser()
    structure_1 = parser_pdb.get_structure('template_1', input_pdb)
    io_1 = PDBIO()
    io_1.set_structure(structure_1)
    #delete marked residues
    io_1.save(output_pdb + ".temp" ,select=SelectResidues(residues_to_delete, chainID))
    parser_pdb = PDBParser()
    structure_2 = parser_pdb.get_structure('template_2', output_pdb + ".temp")
    model = structure_2[model]
    chain = model[chainID]
    # renumber residues from template  to target (according to tt_alignment) has to be done after marked residues are deleted
    # add and then subtract highest id to prevent conflicts in residues ids
    max_residue_id = __get_max_residue_id__(chain)
    for residue in chain:
        for i in range(1, len(template_target_mapping)):
            if residue.id[1] ==  template_target_mapping[i][0]:
                residue.id = (' ', template_target_mapping[i][1] + max_residue_id, ' ')
                break
    for residue in chain:
        residue.id = (' ', residue.id[1] - max_residue_id, ' ')
    io_2 = PDBIO()
    io_2.set_structure(structure_2)
    io_2.save(output_pdb)


def __get_max_residue_id__(chain):
    return max(residue.id[1] for residue in chain)


def create_lowercase_target_fasta_from_tta(target_fasta, name):
    #make target FASTA
    fasta_out = ''
    for res in target_fasta:
        if res != '-':
            fasta_out += res
    ff = open(name, 'w')
    ff.write('>SEQUENCE\n')
    fasta_out = fasta_out.lower()
    print "Target fasta from tta pdb length: " + str(len(fasta_out)-1)
    ff.write(fasta_out[1:])
    ff.close()


class data(object):
    def __init__(self, template_fasta_alignment, target_fasta_alignment, window_size, match_lower_bound ):
        self.template_fasta_alignment = template_fasta_alignment
        self.target_fasta_alignment = target_fasta_alignment
        self.window_size = window_size
        self.match_lower_bound = match_lower_bound


def prepare_alignment_data(alignment_file_name, target_name, template_name):
    import re
    alignment_file = open(name = alignment_file_name, mode = 'r')
    template_seq = 'x'
    target_seq = 'x'
    for line in alignment_file:
        matched_temp = re.match('^'+template_name+' *([A-Za-z-]*) *', line)
        matched_trgt = re.match('^'+target_name+' *([A-Za-z-]*) *', line)
        if matched_temp:
            template_seq += matched_temp.group(1)
        if matched_trgt:
            target_seq += matched_trgt.group(1)
    alignment_file.close()
    print "Template alignment length: " + str(len(template_seq)-1)
    print "Target alignment length: " + str(len(target_seq)-1)
    return data(template_seq, target_seq, None, None)

