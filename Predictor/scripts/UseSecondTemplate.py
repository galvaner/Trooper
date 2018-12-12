# try to fill longer gap by using another template
# vystup je do input.pdb
# ci sa da template pouzit zalezi len na FASTe - v pripade, ze nesedi s PDB-ckom, sa mozu vyzskytnut problemy
# ToDo: ako s koncami? Vyextrahujem z FASTY -> nebezpecenstvo chyby ale inak sa to neda urobit
# ToDo: fosforom sa indexuje v superposition.py
import EmbossNeedle
import Helper
import superposition
from Bio.PDB import *


class SecondTemplateUsage:
    longGaps = []
    enrichedGaps = []
    pdb_directory = "pdbs/"
    modelID = 0
    CONST_FASTAS_DIR = "../fastas/"
    CONST_FASTA_FILE_TYPE = ".fasta"

    def __init__(self, originalReadyForRosettaStruct, target, secondaryTemplate, primaryTemplateChainId):
        self.target = target
        self.originalTemplateChainId = primaryTemplateChainId
        self.template = secondaryTemplate
        self.chainID = Helper.GetChainID(secondaryTemplate)
        self.originalReadyForRosettaStruct = originalReadyForRosettaStruct
        with EmbossNeedle.MyEmboss(self.target, self.template) as myEmboss:
            self.aln = myEmboss.GetAlignment()
        self.longGaps = self.__findLongGaps__(20)  # TODO: set the minimal "big gap" size
        self.enrichedGaps = self.__tryToFillLongGapsFromAnotherTemplate__()
        print self.enrichedGaps

    def __findLongGaps__(self,  minimalGapLength):
        # "identify and return long gaps"
        parser_pdb = PDBParser()
        predicted_part_structure = parser_pdb.get_structure("origReadyForFARFARStr", self.originalReadyForRosettaStruct)
        gaps = []
        lastResId = 0
        predicted_part_structure[self.modelID][self.originalTemplateChainId].child_list.sort(key=lambda x: x.id[1])
        for res in predicted_part_structure[self.modelID][self.originalTemplateChainId]:
            if res.id[1] > lastResId + 1:
                gap = [lastResId + 1, res.id[1], res.id[1] - (lastResId + 1)]  # from(including), till(including + 1), gapLength
                if gap[2] >= minimalGapLength:
                    gaps.append(gap)
            lastResId = res.id[1]
        fastaLength = self.__getTargetFastaLength__()
        if lastResId < fastaLength:
            gap = [lastResId + 1, fastaLength, fastaLength - (lastResId + 1)]  # from(including), till(including + 1), gapLength
            if gap[2] >= minimalGapLength:
                gaps.append(gap)
        return gaps

    def __getResiduesWhichCanBeUSedFromSecondTemplate__(self):
        self.numberedAnotatedTemplate = self.__annotateSecTemplateAlignment__()
        propsToBeUsedFromTemplateSelection = []
        for gap in self.longGaps:
            gapProps = self.__prepareTemplateGapProps__(gap[0], gap[1])
            if self.__decideIfTemplateCanBeUsedForGap__(gapProps, gap[2]):
                propsToBeUsedFromTemplateSelection.append(gapProps)
        return propsToBeUsedFromTemplateSelection

    def __prepareTemplateGapProps__(self, startGap, endGap):
        # map the gap to alignment
        startGapInAln = None
        endGapInAln = None
        for i in range(0, len(self.numberedAnotatedTemplate), 1):
            if self.numberedAnotatedTemplate[i][2] == startGap:
                startGapInAln = i
            if self.numberedAnotatedTemplate[i][2] == endGap:
                endGapInAln = i
                break
        gapCoverageBySecondTempalte = 0
        # check how many residues in gap are conserved
        templateIndexesToFillTheGap = []
        for i in range(startGapInAln, endGapInAln, 1):
            if self.numberedAnotatedTemplate[i][3]:
                gapCoverageBySecondTempalte += 1
                templateIndexesToFillTheGap.append(self.numberedAnotatedTemplate[i][1])
        # check if there is a chance to connect (seq before and after gap has to be conserved)
        beforeGap = 0
        afterGap = 0
        templateIndexesToConnectTheGap = []
        if startGapInAln == 0:
            isStart = True
        else:
            isStart = False
        if endGapInAln == (len(self.numberedAnotatedTemplate)):
            isEnd = True
        else:
            isEnd = False
        for i in range(startGapInAln - 4, startGapInAln):
            if i < 0:
                continue
            if self.numberedAnotatedTemplate[i][3]:
                beforeGap += 1
                templateIndexesToConnectTheGap.append(self.numberedAnotatedTemplate[i][1])
        for i in range(endGapInAln, endGapInAln + 4):
            if i > (len(self.numberedAnotatedTemplate) - 1):
                continue
            if self.numberedAnotatedTemplate[i][3]:
                afterGap += 1
                templateIndexesToConnectTheGap.append(self.numberedAnotatedTemplate[i][1])
        return {"isStart": isStart, "isEnd": isEnd, "beforeGap": beforeGap, "afterGap": afterGap,
                "templateIndexesToConnectTheGap": templateIndexesToConnectTheGap,
                "templateIndexesToFillTheGap": templateIndexesToFillTheGap,
                "gapCoverage": gapCoverageBySecondTempalte}

    def __decideIfTemplateCanBeUsedForGap__(self, gapProperties, gapLength):
        # important "if" - decides the conditions for acceptance of secondary template
        if (((gapProperties["gapCoverage"] / gapLength) < 0.5)  # TODO: set the boundary for template acceptance
                or (gapProperties["beforeGap"] == 0 and not gapProperties["isStart"])
                or (gapProperties["afterGap"] == 0 and not gapProperties["isEnd"])
                or ((gapProperties["beforeGap"] + gapProperties["afterGap"]) <= 3)):
            return False
        else:
            return True

    def __tryToFillLongGapsFromAnotherTemplate__(self):
        # if possible fill the empty gaps
        #self.__coverLongGapsFromSecondaryTemplate__()
        residuesToAdd = self.__getResiduesWhichCanBeUSedFromSecondTemplate__()
        self.callSuperpositionMethod(residuesToAdd)

    def callSuperpositionMethod(self, residuesToAdd):
            for props in residuesToAdd:
                try:
                    superposition.Merge(toStruct=self.originalReadyForRosettaStruct,
                                        fromStruct="../pdbs/" + Helper.TrimPDBName(self.template) + ".pdb",
                                        inGapRes=props['templateIndexesToFillTheGap'],
                                        connectorsRes=props['templateIndexesToConnectTheGap'],
                                        chainIdFrom=self.chainID,
                                        originalTemplateChainId=self.originalTemplateChainId)
                except:
                    print "EXCEPTION: Template fragment could not be aded to original template structure."

    def __annotateSecTemplateAlignment__(self):
        numberedTarget = self.__numberAlignment__(self.aln['target'])
        numberedTemplate = self.__numberAlignment__(self.aln['template'])
        numberedAnotatedTemplate = [] # [template residue, template residue index, target residue index, is successfully aligned to target]
        for i in range(0, len(self.aln['template']), 1):
            if numberedTarget[i][0] != numberedTemplate[i][0]:
                numberedAnotatedTemplate.append([numberedTemplate[i][0], numberedTemplate[i][1], numberedTarget[i][1], False])
            else:
                numberedAnotatedTemplate.append([numberedTemplate[i][0], numberedTemplate[i][1], numberedTarget[i][1], True])
        return numberedAnotatedTemplate

    def __numberAlignment__(self, alnSeq):  # add the position as it was without "-"
        currentNumber = 1
        numberedAln = []
        for residue in alnSeq:
            if residue != '-':
                numberedAln.append([residue, currentNumber])
                currentNumber += 1
            else:
                numberedAln.append(['-', '-'])
        return numberedAln

    def __getTargetFastaLength__(self):
        from Bio import SeqIO
        FastaFile = open(self.CONST_FASTAS_DIR + self.target.upper() + self.CONST_FASTA_FILE_TYPE, 'rU')
        for rec in SeqIO.parse(FastaFile, 'fasta'):
            name = rec.id
            seq = rec.seq
            seqLen = len(rec)
        FastaFile.close()
        return seqLen


#scu = SecondTemplateUsage("input.pdb", "2O45_A", "2O45_A", 'B')
