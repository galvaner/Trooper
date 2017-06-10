import os
import argparse
import NameModifications
import EmbossNeedle

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastaName') #like "46XY_A"
args = vars(parser.parse_args())
targetFastaName = args["fastaName"]
x = 0
for fastaFile in os.listdir('../fastas/'):
    potentialTempalteFastaFormatName = NameModifications.GetFastaNAmeFromFileName(fastaFile)
    with EmbossNeedle.MyEmboss(targetFastaName, potentialTempalteFastaFormatName) as instance:
        similarity = instance.GetSimilarity()

        try:
            if float(similarity) > 60:
                print '******************************'
            if float(similarity) > 60:
                print targetFastaName + " " + potentialTempalteFastaFormatName + ": " + similarity + '\n'
        except:
            x = x + 1
print "Number of exceptions: " + str(x)