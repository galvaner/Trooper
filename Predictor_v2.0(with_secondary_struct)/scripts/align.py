__author__ = 'Rasto'

from Bio.Emboss.Applications import *

def compute_emboss_alignment(input_file_target, input_file_template, output_file):
    #debug print input_file_1 + "       " +  input_file_2
    needle_cli = NeedleCommandline(asequence=input_file_target, \
                               bsequence=input_file_template, \
                               gapopen=10, \
                               gapextend=0.5, \
                               outfile=output_file, \
                               )
    needle_cli()
    #debug print "align.py: FILE ../files/emboss.aln CREATED"
    #debug aln = AlignIO.read("needle_fname", "emboss")
    #debug print aln

