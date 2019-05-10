# Script with emboss alignment functionality
import re


def compute_global_alignment(input_file_target, input_file_template, output_file):
    from Bio.Emboss.Applications import NeedleCommandline
    needle_cli = NeedleCommandline(asequence=input_file_target, \
                               bsequence=input_file_template, \
                               gapopen=10, \
                               gapextend=0.5, \
                               outfile=output_file, \
                               )
    needle_cli()


def compute_local_alignment(input_file_target, input_file_template, output_file):
    from Bio.Emboss.Applications import WaterCommandline
    water_cli = WaterCommandline(  cmd='water',\
                                   asequence=input_file_target, \
                                   bsequence=input_file_template, \
                                   gapopen=10, \
                                   gapextend=0.5, \
                                   outfile=output_file, \
                                   )
    water_cli()


# TODO hlavicka vstupneho FASTA suboru musi mat spravny format, aby parsovanie fungovalo - prerobit, aby to bolo nezavisle???
def parse_alignment(alignment_file_name, target_name, template_name, alnWithLeadingX=True):
    """Get alignment from emboss.water or emboss.needle aln file and its parameters. Return filled Alignment object. Returned alignments starts with x as dummy character, so the indexing matches"""
    import re
    alignment_file = open(name=alignment_file_name, mode='r')
    template_seq = ''
    target_seq = ''
    if alnWithLeadingX:
        template_seq = 'x'
        target_seq = 'x'
    target_start = -1
    target_end = -1
    template_start = -1
    template_end = -1
    for line in alignment_file:
        matched_temp = re.match('^'+template_name+' *([0-9]+) *([A-Za-z-]*) *([0-9]+)$', line)
        matched_trgt = re.match('^'+target_name+' *([0-9]+) *([A-Za-z-]*) *([0-9]+)$', line)
        if matched_temp:
            template_seq += matched_temp.group(2)
            if template_start == -1:
                template_start = int(matched_temp.group(1))
            template_end = int(matched_temp.group(3))
        if matched_trgt:
            target_seq += matched_trgt.group(2)
            if target_start == -1:
                target_start = int(matched_trgt.group(1))
            target_end = int(matched_trgt.group(3))
    alignment_file.close()
    template_fasta_alignment = template_seq
    target_fasta_alignment = target_seq
    return Alignment(template_fasta_alignment, template_start, template_end, target_fasta_alignment, target_start, target_end)


class Alignment(object):
    def __init__(self, template_fasta_alignment, template_start_index, template_end_index, target_fasta_alignment, target_start_index, target_end_index):
        self.template_fasta_alignment = template_fasta_alignment
        self.target_fasta_alignment = target_fasta_alignment
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        self.template_start_index = template_start_index
        self.template_end_index = template_end_index


def get_similarity_from_alignment(file_name):
    """Get similarity and gap from emboss.water or emboss.needle aln file. Returns list of floats [similarity, gap]."""
    file = open(name=file_name, mode='r')
    similarity = 'no_match'
    gap = 'no_match'
    for line in file:
        matched = re.match('^# Similarity:[^\(]*\((.*)%\)', line)
        matched_gap = re.match('^# Gaps:[^\(]*\((.*)%\)', line)
        if matched:
            similarity = matched.group(1)
        if matched_gap:
            gap = matched_gap.group(1)
    file.close()
    return [float(similarity), float(gap)]