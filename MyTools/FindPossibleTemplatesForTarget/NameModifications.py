# XXXX_anything, I want firts 4 characters in lowercase
def TrimPDBName(fastaName):
    return fastaName[:4].lower()

# XXXX_A_anything, I want A
def GetChainID(fastaName):
    return fastaName[5:6]

def GetFastaNAmeFromFileName(fileName):
    return fileName[:6]