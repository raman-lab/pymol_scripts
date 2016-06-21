import glob
import sys
import re 

def getVars(inFile):
    # read in file and only keep important lines
    allLines = [line.rstrip('\n') for line in open(inFile)]
    inputs = []
    vPairInputs = []
    for entry in allLines:
        if re.search('[a-zA-Z]', entry) and '\t' not in entry and '#' not in entry:
            inputs.append(entry)
            inLineList = entry.split(': ')
        if '\t' in entry and '#' not in entry:
            vPairInputs.append(entry)

    # get variables from input
    pdb = inputs[0].split(': ')[1]
    lig = inputs[1].split(': ')[1]
    const = inputs[2].split(': ')[1]
    print const
    alignedArea = inputs[3].split(': ')[1]
    if alignedArea.lower() == 'no':
        alignedArea = 'Ax'
    cAtomsStr = inputs[4].split(': ')[1]
    cAtomsList = cAtomsStr.split(', ')

    vAtomsDict = {}
    for i,entry in enumerate(vPairInputs):
        print entry
        vAtomsTmpStr = entry.split(': ')[1]
        vAtomsTmpList = vAtomsTmpStr.split(', ')
        vAtomsDict[i] = vAtomsTmpList

    return pdb, lig, const, alignedArea, cAtomsList, vAtomsDict

def reformatConf(confFile):
    # reformat conformer file so that output is Rosetta compatible
    with open(confFile,'r') as f:
        lines = f.readlines()
    with open(confFile,'w') as f:
        for line in lines:
            lineStr = ''.join(line)
            if 'HETATM' in lineStr:
                lineStr = lineStr.replace('UNK   900', 'LG1 X   1')
                f.write(lineStr)

def alignLigs(pdb,lig,const,alignedArea,cAtomsList,vAtomsDict):
    # define pymol variables
    pdbM = pdb.split('.')[0]
    base = '//X/LG1`1'
    fill = '////'
    idM = base[4:base.find('`')]
    
    # loop over all conformers and all pairings
    confFileList = glob.glob(lig + '_*.pdb')
    for i,entry in enumerate(confFileList):
        reformatConf(entry)
        for j in range(0,len(vAtomsDict)):
            confName = entry.split('.')[0]
            
            # load files
            cmd.reinitialize()
            cmd.load(pdb)
            cmd.select('ogLig', 'resname ' + idM)
            cmd.load(entry)

            # align atoms
            print confName + '_P' + str(j) + '_C' + str(i)
            if const == 'pdbAtoms':
                if len(cAtomsList) == 3:
                    cmd.pair_fit(confName+fill+vAtomsDict[j][0], pdbM+fill+cAtomsList[0],
                                 confName+fill+vAtomsDict[j][1], pdbM+fill+cAtomsList[1],
                                 confName+fill+vAtomsDict[j][2], pdbM+fill+cAtomsList[2])
                elif len(cAtomsList) == 4:
                    cmd.pair_fit(confName+fill+vAtomsDict[j][0], pdbM+fill+cAtomsList[0],
                                 confName+fill+vAtomsDict[j][1], pdbM+fill+cAtomsList[1],
                                 confName+fill+vAtomsDict[j][2], pdbM+fill+cAtomsList[2],
                                 confName+fill+vAtomsDict[j][3], pdbM+fill+cAtomsList[3])
                elif len(cAtomsList) == 5:
                    cmd.pair_fit(confName+fill+vAtomsDict[j][0], pdbM+fill+cAtomsList[0],
                                 confName+fill+vAtomsDict[j][1], pdbM+fill+cAtomsList[1],
                                 confName+fill+vAtomsDict[j][2], pdbM+fill+cAtomsList[2],
                                 confName+fill+vAtomsDict[j][3], pdbM+fill+cAtomsList[3],
                                 confName+fill+vAtomsDict[j][4], pdbM+fill+cAtomsList[4])
                elif len(cAtomsList) == 6:
                    cmd.pair_fit(confName+fill+vAtomsDict[j][0], pdbM+fill+cAtomsList[0],
                                 confName+fill+vAtomsDict[j][1], pdbM+fill+cAtomsList[1],
                                 confName+fill+vAtomsDict[j][2], pdbM+fill+cAtomsList[2],
                                 confName+fill+vAtomsDict[j][3], pdbM+fill+cAtomsList[3],
                                 confName+fill+vAtomsDict[j][4], pdbM+fill+cAtomsList[4],
                                 confName+fill+vAtomsDict[j][5], pdbM+fill+cAtomsList[5])
                    
                    
            elif const == 'ligAtoms':
                if len(cAtomsList) == 3:
                    cmd.pair_fit(confName+fill+cAtomsList[0], pdbM+fill+vAtomsDict[j][0],
                                 confName+fill+cAtomsList[1], pdbM+fill+vAtomsDict[j][1],
                                 confName+fill+cAtomsList[2], pdbM+fill+vAtomsDict[j][2])
                elif len(cAtomsList) == 4:
                    cmd.pair_fit(confName+fill+cAtomsList[0], pdbM+fill+vAtomsDict[j][0],
                                 confName+fill+cAtomsList[1], pdbM+fill+vAtomsDict[j][1],
                                 confName+fill+cAtomsList[2], pdbM+fill+vAtomsDict[j][2],
                                 confName+fill+cAtomsList[3], pdbM+fill+vAtomsDict[j][3])
                elif len(cAtomsList) == 5:
                    cmd.pair_fit(confName+fill+cAtomsList[0], pdbM+fill+vAtomsDict[j][0],
                                 confName+fill+cAtomsList[1], pdbM+fill+vAtomsDict[j][1],
                                 confName+fill+cAtomsList[2], pdbM+fill+vAtomsDict[j][2],
                                 confName+fill+cAtomsList[3], pdbM+fill+vAtomsDict[j][3],
                                 confName+fill+cAtomsList[4], pdbM+fill+vAtomsDict[j][5])
                elif len(cAtomsList) == 6:
                    cmd.pair_fit(confName+fill+cAtomsList[0], pdbM+fill+vAtomsDict[j][0],
                                 confName+fill+cAtomsList[1], pdbM+fill+vAtomsDict[j][1],
                                 confName+fill+cAtomsList[2], pdbM+fill+vAtomsDict[j][2],
                                 confName+fill+cAtomsList[3], pdbM+fill+vAtomsDict[j][3],
                                 confName+fill+cAtomsList[4], pdbM+fill+vAtomsDict[j][4],
                                 confName+fill+cAtomsList[5], pdbM+fill+vAtomsDict[j][5])
            print ''

            # generate new pdb file with orginal ligand replaced
            cmd.remove('ogLig')
            cmd.set('pdb_use_ter_records', 0)
            cmd.save(pdbM + '_' + lig + '_' + alignedArea + '_P' + str(j) + '_C' + str(i) + '.pdb')


# execute
print 'Starting'
print ''
inFile = sys.argv[-1]
pdb, lig, const, alignedArea, cAtomsList, vAtomsDict = getVars(inFile)
alignLigs(pdb, lig, const, alignedArea, cAtomsList, vAtomsDict)
print 'Done'                   
            
