#get secondary structure from PDBs
import sys
import json
outputPath = './secondary_structures/'
f = sys.stdin #open('json.temp', 'r')
jsonText = json.load(f)
for chain in jsonText['chains'].iteritems():
    outFileName = unicode.upper(jsonText['metadata']['str_id']) + "_" + chain[0][6] + ".secstr"
    outFile = open(outputPath + outFileName, 'w+')
    outFile.write(chain[1]['bseq'] + '\n')
    outFile.write(chain[1]['sstr'])
    outFile.close()