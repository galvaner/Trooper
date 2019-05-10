from Bio.Emboss.Applications import *
import os
from Bio import AlignIO
import re


class MyEmboss:
    '''Compute global or local fasta files alignment and provide functions to get results (aligned sequences, aln similarity percentage, aln gap percentage)'''
    CONST_FASTAS_DIR = "./fastas/"  # "./fasta_secstr/"
    CONST_FASTA_FILE_TYPE = ".fasta"
    CONST_TEMP_ALN_FILE_NAME = "aln.emboss"

    def __init__(self, targetFastaName, templateFastaName, typeOfAlignemnt = 'global', target_fasta_file="", template_fasta_file="",gapOpen=10, gapExtend=0.5):
        self.targetFastaName = targetFastaName
        self.templateFastaName = templateFastaName
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
        if target_fasta_file == "":
            self.inputTargetFile = self.CONST_FASTAS_DIR + targetFastaName + self.CONST_FASTA_FILE_TYPE
        else:
            self.inputTargetFile = target_fasta_file
        if template_fasta_file == "":
            self.inputTemplateFile = self.CONST_FASTAS_DIR + templateFastaName + self.CONST_FASTA_FILE_TYPE
        else:
            self.inputTemplateFile = template_fasta_file
        if typeOfAlignemnt == 'local':
            self.__computeLocalAlignment__()
        else:
            self.__computeGlobalAlignment__()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.CONST_TEMP_ALN_FILE_NAME)

    def __computeGlobalAlignment__(self):
        try:
            needle_cli = NeedleCommandline(asequence=self.inputTargetFile, \
                                   bsequence=self.inputTemplateFile, \
                                   gapopen=self.gapOpen, \
                                   gapextend=self.gapExtend, \
                                   outfile=self.CONST_TEMP_ALN_FILE_NAME \
                                   )
            needle_cli()
        except:
            #print "Needle returned exception - dummy file created (target: " + self.inputTargetFile + "), (template: " + self.inputTemplateFile + ") "
            self.__createDummyNeedleResult__()


    def __computeLocalAlignment__(self):
        try:
            from Bio.Emboss.Applications import WaterCommandline
            water_cli = WaterCommandline(cmd='water', \
                                         asequence=self.inputTargetFile, \
                                         bsequence=self.inputTemplateFile, \
                                         gapopen=10, \
                                         gapextend=0.5, \
                                         outfile=self.CONST_TEMP_ALN_FILE_NAME, \
                                         )
            water_cli()
        except:
            print "Needle returned exception - dummy file created (target: " + self.inputTargetFile + "), (template: " + self.inputTemplateFile + ") "
            self.__createDummyNeedleResult__()


    def __createDummyNeedleResult__(self):
        with open(self.CONST_TEMP_ALN_FILE_NAME, 'w') as aln:
            aln.write('Length: 12\
                        # Identity:       4/12 (00.0%)\
                        # Similarity:     4/12 (00.0%)\
                        # Gaps:           4/12 (00.0%)\
                        # Score: 4.0')


    def __regexpForSimilarityOrGap__(self, startWord):
        with open(self.CONST_TEMP_ALN_FILE_NAME) as aln:
            for line in aln:
                m = re.search(r'# ' + startWord + r':.*\(([0-9.]*)%\)', line)
                if m is not None:
                    return m.group(1)


    def GetSimilarity(self):
        similarity = self.__regexpForSimilarityOrGap__('Similarity')
        return similarity


    def GetGaps(self):
        gaps = self.__regexpForSimilarityOrGap__('Gaps')
        return gaps


    def GetAlignment(self):
        align = AlignIO.read(self.CONST_TEMP_ALN_FILE_NAME, "emboss")
        target = list(align[0].seq)
        template = list(align[1].seq)
        return {'target': target, 'template': template}

    def ParseSemiGlobalAlignment(self):
        """Get similarity and indexes to template in alignment a-sequence must be target, b-sequence must be template"""
        import re
        alignment_file = open(name=self.CONST_TEMP_ALN_FILE_NAME, mode='r')
        a_seq = ''
        b_seq = ''
        b_start, b_end, a_start, a_end = 0,0,0,0
        matched_a_seq = None
        matched_b_seq = None
        for line in alignment_file:
            if matched_a_seq is None:
                matched_a_seq = re.match('^[A-Za-z]* *([0-9]+) *([A-Za-z-]*) *([0-9]+)$', line)
                continue
            if matched_b_seq is None:
                matched_b_seq = re.match('^[A-Za-z]* *([0-9]+) *([A-Za-z-]*) *([0-9]+)$', line)
            if matched_a_seq and matched_b_seq:
                if int(matched_a_seq.group(1)) > 0 and int(matched_a_seq.group(1)) != int(matched_a_seq.group(3)):
                    a_seq += matched_a_seq.group(2)
                    b_seq += matched_b_seq.group(2)
                    if int(matched_a_seq.group(1)) == 1:
                        b_start = int(matched_b_seq.group(1))
                    if int(matched_a_seq.group(1)) != int(matched_a_seq.group(3)):
                        b_end = int(matched_b_seq.group(3))
                matched_b_seq = None
                matched_a_seq = None
        alignment_file.close()

        if a_seq == '' or b_seq == '':
            return a_seq, b_seq, 0, 0, 0

        leading_dash_count = 0
        while a_seq[leading_dash_count] == '-':
            leading_dash_count += 1
        a_seq = a_seq[leading_dash_count:]
        b_seq = b_seq[leading_dash_count:]
        b_start = b_start + leading_dash_count

        trailing_dash_count = 0
        while a_seq[-trailing_dash_count-1] == '-':
            trailing_dash_count += 1
        if trailing_dash_count > 0:
            a_seq = a_seq[:-trailing_dash_count]
            b_seq = b_seq[:-trailing_dash_count]
            b_end = b_end - trailing_dash_count

        # get similarity
        count_equals = 0
        for i in xrange(0, len(a_seq)):
            if a_seq[i] == b_seq[i]:
                count_equals += 1
        similarity = float(count_equals)/len(a_seq)

        return a_seq, b_seq, b_start, b_end, similarity





