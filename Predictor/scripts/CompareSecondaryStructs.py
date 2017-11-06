# neporovnavam [,],{,}
# script is run from "scripts" folder

import glob
import os

class ComparationStruct:
    def __init__(self):
        self.falsePositives = 0
        self.falseNegatives = 0
        self.truePositives = 0
        self.trueNegatives = 0

class SecStruComparator:
    def compare(self, original, prediction):
        originalString = self.__getBracketStringFromDefaultFormat__(original)
        predictionString = self.__getBracketStringFromFARFARFormat__(prediction)
        originalBP = self.__findBasePairs__(originalString)
        predictionBP = self.__findBasePairs__(predictionString)
        comparation = self.__compareBasePairs__(originalBP, predictionBP)
        return comparation

    # returns base pairs as dictionary keys e.g. bp[10 20] from string of brackets
    def __findBasePairs__(self, inputString):
        basePairs = {}
        for i in range(0, len(inputString), 1):
            if inputString[i] == '(':
                rightIndex = self.__getBracketRightPairIndex__(i, inputString)
                basePairs[str(i) + " " + str(rightIndex)] = "base pair"
        return basePairs

    # returns right index of bracket paired to the given left one
    def __getBracketRightPairIndex__(self, currentPosition, secStru):
        bracketCounter = 0
        for i in xrange(currentPosition, len(secStru), 1):
            if secStru[i] == '(':
                bracketCounter += 1
            if secStru[i] == ')':
                bracketCounter -= 1
            if bracketCounter == 0:
                return i
        raise Exception("Inpaired brckets")

    # ignore first line, base pairs are in second one
    def __getBracketStringFromDefaultFormat__(self, filePath):
        basePairsStr = ""
        structFile = open(filePath, 'rU')
        lineCounter = 0
        for line in structFile:
            lineCounter += 1
            if lineCounter > 2:
                raise Exception("Third row in input file " + filePath + " found.")
            if lineCounter == 1:
                continue
            basePairsStr = line
        return basePairsStr

    # ignore second line, base pairs are in  first one
    def __getBracketStringFromFARFARFormat__(self, filePath):
        basePairsStr = ""
        structFile = open(filePath, 'rU')
        lineCounter = 0
        for line in structFile:
            lineCounter += 1
            if lineCounter > 2:
                raise Exception("Third row in input file " + filePath + " found.")
            if lineCounter == 2:
                continue
            basePairsStr = line
        return basePairsStr

    # compares base pairs and return true/false positives/negatives
    def __compareBasePairs__(self, originalStr, predictedStr):
        comp = ComparationStruct()
        for predictedStructureBP in predictedStr.keys():
            if originalStr.has_key(predictedStructureBP):
                comp.truePositives += 1
            if not originalStr.has_key(predictedStructureBP):
                comp.falsePositives += 1
        for originalStructureBP in originalStr.keys():
            if not predictedStr.has_key(originalStructureBP):
                comp.falseNegatives += 1
        # ToDo: trueNegatives
        return comp

# loop through all predicted structures and compare them to experimentally created secondary structures
def RunComparassionInMetacentrum():
    cmp = SecStruComparator()
    results_list = [('structure', 'TruePos', 'FalseNeg', 'FalsePos', 'TrueNeg')]
    results_dict = {} # dictionary won't be complete (each prediction can be unique (RNA_Fold nondeterministic algorithm))
    for predictionPath in glob.glob('../prediction/*/*/rosetta/secstruForRosetta.secstr'):
        path = predictionPath.replace('\\', '/')
        splitedName = str.split(path, '/')
        structureName = splitedName[3][0:6]
        originalPath = "../secondary_structures/" + structureName + ".secstr"
        try:
            result = cmp.compare(originalPath, path)
            results_list.append((structureName, result.truePositives, result.falseNegatives, result.falsePositives, result.trueNegatives))
            results_dict[structureName] = (result.truePositives, result.falseNegatives, result.falsePositives, result.trueNegatives)
        except:
            print "Exception occured when comparing", structureName, "secondary structure."
    file = open('SecStruComparation.txt', 'w')
    file.write('structure TruePos FalseNeg FalsePos TrueNeg\n')
    # Use dictionary to make it more easy. I think diferences between the best and the worst prediction are minimal
    for rec in results_dict:
        file.write(rec + " " + str(results_dict[rec][0]) + " " + str(results_dict[rec][1]) + " " + str(results_dict[rec][2]) + "\n")
    file.close()

def TestInLocalFolder():
    cmp = SecStruComparator()
    # 1EHZ_A  tends to have bad results in tertiary structures prediction as a target
    # result = cmp.compare('./1EHZ_A.secstr', './secstruForRosetta_1EHZ_A.secstr')
    # 2GDI_X tends to have bad results in tertiary structures prediction as a target
    result = cmp.compare('./2GDI_X.secstr', './secstruForRosetta_2GDI_X.secstr')
    print 'True Positives:', result.truePositives
    print 'False Negatives', result.falseNegatives
    print 'False Positives', result.falsePositives
    print 'True Negatives: (not rated)', result.trueNegatives

#TestInLocalFolder()
RunComparassionInMetacentrum()