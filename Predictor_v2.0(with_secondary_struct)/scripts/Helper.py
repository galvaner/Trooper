# XXXX_anything, I want firts 4 characters in lowercase
def TrimPDBName(fastaName):
    return fastaName[:4].lower()


# XXXX_A_anything, I want A
def GetChainID(fastaName):
    return fastaName[5:6]


def GetFastaNAmeFromFileName(fileName):
    return fileName[:6]

CONST_SEC_STR_FOLDER = '../secondary_structures/'

def __makeSecStrName__(rawName):
    return CONST_SEC_STR_FOLDER + rawName + ".secstr"

#expects fastaName like XXXX_A.anything
def LoadSecondaryStructureToString(fastaName):
    rawName = GetFastaNAmeFromFileName(fastaName)
    with open(__makeSecStrName__(rawName), 'r') as file:
        firstLine = file.next()
    return firstLine


def ListOfPairsToFiles(zeroIndexFileName, listOfPairs):
    fileZero = open(zeroIndexFileName, 'w')
    fileZero.write(">SEQ:\n")
    for pair in listOfPairs:
        fileZero.write(pair[0])
    fileZero.write('\n')
    for pair in listOfPairs:
        fileZero.write(pair[1])
    fileZero.close()
