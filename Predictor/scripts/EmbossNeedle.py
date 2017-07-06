from Bio.Emboss.Applications import *
import os
from Bio import AlignIO
import re


class MyEmboss:

    CONST_FASTAS_DIR = "../fastas/"  # "./fasta_secstr/"
    CONST_FASTA_FILE_TYPE = ".fasta"
    CONST_TEMP_ALN_FILE_NAME = "aln.emboss"

    def __init__(self, targetFastaName, templateFastaName):
        self.inputTargetFile = self.CONST_FASTAS_DIR + targetFastaName + self.CONST_FASTA_FILE_TYPE
        self.inputTemplateFile = self.CONST_FASTAS_DIR + templateFastaName + self.CONST_FASTA_FILE_TYPE
        self.__computeEmbossAlignment__()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.CONST_TEMP_ALN_FILE_NAME)

    def __computeEmbossAlignment__(self):
        try:
            needle_cli = NeedleCommandline(asequence=self.inputTargetFile, \
                                   bsequence=self.inputTemplateFile, \
                                   gapopen=10, \
                                   gapextend=0.5, \
                                   outfile=self.CONST_TEMP_ALN_FILE_NAME, \
                                   )
            needle_cli()
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