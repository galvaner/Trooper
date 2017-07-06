import EmbossNeedle
import Helper
import os


class SecStrPredictor:

    def __init__(self, target, template, outputDirectory):  # target, template like XXXX_A
        self.target = target
        self.template = template
        self.PREDICTED_TARGET_SECSTR = outputDirectory + self.target + '_predicted_sec_struct.secstr'
        self.TEMP_FOR_RNA_FOLD = outputDirectory + "template_secstr_fasta.temp"
        self.RNA_FOLD_SCRIPT = outputDirectory + "RNAfold_run.sh"
        self.SEC_STRU_FOR_ROSETTA = outputDirectory + "secstruForRosetta.secstr"
        self.__predictSecondaryStructure__()

    def __predictSecondaryStructure__(self):
        with EmbossNeedle.MyEmboss(self.target, self.template) as myEmboss:
            self.aln = myEmboss.GetAlignment()
        templateNumberedAln = self.__numberAlignment__(self.aln['template'])
        self.enriBothAln = self.__enrichTemplateAlnWithSecondaryStructureAndTargetAln__(templateNumberedAln)  # contains: [ templateResidue, indexOfTemplResBeforeAln, secStr, targetResidue]
        self.__checkIfBracketsArePaired__()
        self.listPreparedForPrediction = self.__prepareForSecPrediction__()  # modifies self.
        Helper.ListOfPairsToFiles(self.TEMP_FOR_RNA_FOLD, self.listPreparedForPrediction)
        self.__runRnaFold__()
        #os.remove(self.TEMP_FOR_RNA_FOLD)
        self.__modifySecstructForRosetta__()

    def __runRnaFold__(self):
        file = open(self.RNA_FOLD_SCRIPT, 'wb')
        file.write('#!/bin/bash \n')
        file.write('RNAfold -C < ' + self.TEMP_FOR_RNA_FOLD + ' > ' + self.PREDICTED_TARGET_SECSTR)
        file.close()
        #self.__windowsBashScriptCall__()
        #os.remove(self.RNA_FOLD_SCRIPT)
        self.__unixBashScriptCall__()
        #os.remove("SEQ__ss.ps") # RnaFold 'temp' file


    def __windowsBashScriptCall__(self):
        import subprocess
        subprocess.call(['C:\\cygwin64\\bin\\bash', './' + self.RNA_FOLD_SCRIPT], shell=True)

    def __unixBashScriptCall__(self):
        import subprocess
        subprocess.call(['chmod', '+x', './' + self.RNA_FOLD_SCRIPT])
        subprocess.call(['./' + self.RNA_FOLD_SCRIPT])

    def __prepareForSecPrediction__(self):
        result = []
        self.__makeSecStrRNAfoldFriendly__()
        for item in self.enriBothAln:
            if item[0] != '-' and item[3] != '-':
                result.append([item[3], item[2]])
                continue
            if item[0] == '-' and item[3] != '-':
                result.append([item[3], '.'])
                continue
            if item[0] != '-' and item[3] == '-':
                continue
            print "error in __prepareForSecPrediction__"
        return result

    def __makeSecStrRNAfoldFriendly__(self):
        for i in xrange(0, len(self.enriBothAln), 1):
            if self.enriBothAln[i][2] not in ['(', ')', '.']:
                self.enriBothAln[i][2] = '.'
        for i in xrange(0, len(self.enriBothAln), 1):
            if self.enriBothAln[i][3] == '-':  # item[0] can not be '-', if gap in target aln
                if self.enriBothAln[i][2] == '(':
                    self.enriBothAln[i][2] = '.'
                    self.__removeRightBracket__(i)
                if self.enriBothAln[i][2] == ')':
                    self.enriBothAln[i][2] = '.'
                    self.__removeLeftBracket__(i)

    def __removeRightBracket__(self, currentPosition):
        bracketCounter = 0
        for i in xrange(currentPosition, len(self.enriBothAln), 1):
            if self.enriBothAln[i][2] == '(':
                bracketCounter += 1
            if self.enriBothAln[i][2] == ')':
                bracketCounter -= 1
            if bracketCounter == -1:
                self.enriBothAln[i][2] = '.'
                break

    def __removeLeftBracket__(self, currentPosition):
        bracketCounter = 0
        for i in xrange(currentPosition, -1, -1):
            if self.enriBothAln[i][2] == '(':
                bracketCounter -= 1
            if self.enriBothAln[i][2] == ')':
                bracketCounter += 1
            if bracketCounter == -1:
                self.enriBothAln[i][2] = '.'
                break

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

    def __enrichTemplateAlnWithSecondaryStructureAndTargetAln__(self, numberedAln):  # returns: [templateResidue, indexOfTemplResBeforeAln, secStr, targetResidue]
        secStr = Helper.LoadSecondaryStructureToString(self.template)
        self.__checkIfBracketsArePairedGiveStr__(secStr)
        enrichedSeqence = []
        __bracketSequence = ""
        __notBrackedSequence = ""
        positionInAln = 0
        for pair in numberedAln:
            if pair[0] != '-':
                enrichedSeqence.append([pair[0], pair[1], secStr[pair[1] - 1], self.aln['target'][positionInAln]])  # [A, 5, (]
                __bracketSequence += secStr[pair[1] - 1]
                __notBrackedSequence = __notBrackedSequence + str(pair[1]) + " "
            else:
                enrichedSeqence.append(['-', '-', '-', self.aln['target'][positionInAln]])
            positionInAln += 1
        #print  __notBrackedSequence
        self.__checkIfBracketsArePairedGiveStr__(__bracketSequence)
        return enrichedSeqence

    def __checkIfBracketsArePaired__(self):
            bracketCounter = 0
            for i in xrange(0, len(self.enriBothAln), 1):
                if self.enriBothAln[i][2] == '(':
                    bracketCounter += 1
                if self.enriBothAln[i][2] == ')':
                    bracketCounter -= 1
                if bracketCounter < 0:
                    raise Exception("Inpaired brckets")
            if bracketCounter != 0:
                raise Exception("Inpaired brckets")

    def __checkIfBracketsArePairedGiveStr__(self, structure):
        bracketCounter = 0
        for i in xrange(0, len(structure), 1):
            if structure[i] == '(':
                bracketCounter += 1
            if structure[i] == ')':
                bracketCounter -= 1
            if bracketCounter < 0:
                raise Exception("Inpaired brckets")
        if bracketCounter != 0:
            raise Exception("Inpaired brckets")

    def __modifySecstructForRosetta__(self):
        file = open(self.PREDICTED_TARGET_SECSTR, 'r')
        fastaPart = ""
        secStruPart = ""
        for i, line in enumerate(file):
            if  i == 1:
                fastaPart = line.lower()
            if i == 2:
                secStruPart = line
        file.close()
        file = open(self.SEC_STRU_FOR_ROSETTA, 'w')
        file.write(secStruPart[:len(fastaPart)] + '\n')
        file.write(fastaPart)
        file.close()

    def FixUnpairedChosenWorkingResidues(self, rangePairList): # expects list of lists of pairs chosen residues (edge values included)
        file = open(self.SEC_STRU_FOR_ROSETTA, 'r')
        for i, line in enumerate(file):
            if i == 0:
                self.predictedSecStruOfTarget = line
        file.close()
        self.__numberMatchingParenthesis__()
        print self.numberedPredictedSecStruct
        workingResiduesToBeAdded = []
        for pair in rangePairList:
            print pair
            for i in xrange(pair[0], pair[1]):
                index = i - 1  # numbering of residues is indexed from 1
                if not self.__isPairedBrackedAlsoChosenAsWorkingRes__(self.numberedPredictedSecStruct[index], rangePairList):
                    workingResiduesToBeAdded.append(int(self.numberedPredictedSecStruct[index]) + 1)  # numbering of residues is indexed from 1
        returnvalue = self.__convertToStringOfWorkingResidues__(workingResiduesToBeAdded)
        print returnvalue
        return returnvalue

    def __isPairedBrackedAlsoChosenAsWorkingRes__(self, i, listOfPairs):
        if i == '-':
            return True
        for pair in listOfPairs:
            if ((pair[0] - 1) <= int(i)) and ((pair[1] - 1) >= int(i)):  # numbering of residues is indexed from 1 TODO odcitavam spravne jednicku???
                return True
        return False

    def __convertToStringOfWorkingResidues__(self, workingResToBeAddedSet):
        outputString = ""
        for i in workingResToBeAddedSet:
            outputString += " " + str(i)
        return outputString

    def __numberMatchingParenthesis__(self):
        self.numberedPredictedSecStruct = []
        for i in xrange(0, len(self.predictedSecStruOfTarget), 1):
            if self.predictedSecStruOfTarget[i] == '.':
                self.numberedPredictedSecStruct.append('-')
            if self.predictedSecStruOfTarget[i] == '(':
                self.numberedPredictedSecStruct.append(str(self.__matchRightBracket__(i)))
            if self.predictedSecStruOfTarget[i] == ')':
                self.numberedPredictedSecStruct.append(str(self.__matchLeftBracket__(i)))

    def __matchRightBracket__(self, currentPosition):
        bracketCounter = 0
        for i in xrange(currentPosition, len(self.predictedSecStruOfTarget), 1):
            if self.predictedSecStruOfTarget[i] == '(':
                bracketCounter += 1
            if self.predictedSecStruOfTarget[i] == ')':
                bracketCounter -= 1
            if bracketCounter == 0:
                return i
        raise Exception("Inpaired brckets")

    def __matchLeftBracket__(self, currentPosition):
        bracketCounter = 0
        for i in xrange(currentPosition, -1, -1):
            if self.predictedSecStruOfTarget[i] == '(':
                bracketCounter -= 1
            if self.predictedSecStruOfTarget[i] == ')':
                bracketCounter += 1
            if bracketCounter == 0:
                return i
        raise Exception("Inpaired brckets")




def RunSecStrPRediction(target, template, outputDirectory):  # "2GIS_A", "3V7E_D"
    myPred = SecStrPredictor(target, template, outputDirectory)
    return myPred

# object = RunSecStrPRediction("3V7E_C", "2GIS_A", "./")
# a = object.FixUnpairedChosenWorkingResidues([[1, 36] ,[116, 126]])
# b = object.FixUnpairedChosenWorkingResidues([[5, 36], [116, 121], [37, 115]])
# print a