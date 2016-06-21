import glob
import sys
import re
import random

def getVars(inFile):
    # read in file and only keep important lines
    allLines = [line.rstrip('\n') for line in open(inFile)]
    inputs = []
    vPairInputs = []
    for entry in allLines:
        if re.search('[a-zA-Z]', entry) and '\t' not in entry and '#' not in entry:
            inputs.append(entry)
            inLineList = entry.split(': ')
        if '\t' in entry:
            vPairInputs.append(entry)

    # get variables from input
    pdb = inputs[0].split(': ')[1]
    lig = inputs[1].split(': ')[1]
    atomNum = inputs[2].split(': ')[1]
    comboNum = inputs[3].split(': ')[1]
    return pdb, lig, atomNum, comboNum

def proteinLigandAtoms(pdb):
    # get all atoms of ogLig
    ogLigAtomList = []
    base = '//X/LG1`1'
    idM = base[4:base.find('`')]
    with open(pdb, 'r') as pdbf:
        for line in pdbf:
                if line.split(None, 1)[0] == 'HETATM' and line.split(None, 4)[3] == idM:
                        platm = line.split()[2]
                        ogLigAtomList.append(platm)
    return ogLigAtomList

def testLigandAtoms(ligFile):
    # get all atoms of test ligand
    lig_list = []
    with open(ligFile, 'r') as ligf:
        for line in ligf:
                if line.split(None, 1)[0] == 'HETATM':
                        atm = line.split()[2]
                        lig_list.append(atm)
    return lig_list

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

def alignLigs(pdb, lig, atomNum, comboNum, ogLigAtomList):
    # define pymol variables
    pdbM = pdb.split('.')[0]
    base = '//X/LG1`1'
    fill = '////'
    idM = base[4:base.find('`')]
    
    # loop over all conformers
    confFileList = glob.glob(lig + '_*.pdb')
    for i,entry in enumerate(confFileList):
        reformatConf(entry)
        
        # get atomNum random conformer atoms 
        lig_list = testLigandAtoms(entry)
        confAtoms = random.sample(lig_list,int(atomNum))

        # generate comboNum random pairings per conformer
        for j in range(0,int(comboNum)):
            confName = entry.split('.')[0]
            ogAtoms = random.sample(ogLigAtomList, int(atomNum))
                                    
            # load files
            cmd.reinitialize()
            cmd.load(pdb)
            cmd.select('ogLig', 'resname ' + idM)
            cmd.load(entry)

            # align atoms
            print lig + '_P' + str(j) + '_C' + str(i)
            print 'Conformer Atoms Used: ' + ', '.join(confAtoms)
            print 'Original Ligand Atoms Used: ' +', '.join(ogAtoms)
            if int(atomNum) == 3:
                cmd.pair_fit(confName+fill+confAtoms[0], pdbM+fill+ogAtoms[0],
                             confName+fill+confAtoms[1], pdbM+fill+ogAtoms[1],
                             confName+fill+confAtoms[2], pdbM+fill+ogAtoms[2])
            elif int(atomNum) == 4:
                cmd.pair_fit(confName+fill+confAtoms[0], pdbM+fill+ogAtoms[0],
                             confName+fill+confAtoms[1], pdbM+fill+ogAtoms[1],
                             confName+fill+confAtoms[2], pdbM+fill+ogAtoms[2],
                             confName+fill+confAtoms[3], pdbM+fill+ogAtoms[3])
            elif int(atomNum) == 5:
                cmd.pair_fit(confName+fill+confAtoms[0], pdbM+fill+ogAtoms[0],
                             confName+fill+confAtoms[1], pdbM+fill+ogAtoms[1],
                             confName+fill+confAtoms[2], pdbM+fill+ogAtoms[2],
                             confName+fill+confAtoms[3], pdbM+fill+ogAtoms[3],
                             confName+fill+confAtoms[4], pdbM+fill+ogAtoms[4])
            elif int(atomNum) == 4:
                cmd.pair_fit(confName+fill+confAtoms[0], pdbM+fill+ogAtoms[0],
                             confName+fill+confAtoms[1], pdbM+fill+ogAtoms[1],
                             confName+fill+confAtoms[2], pdbM+fill+ogAtoms[2],
                             confName+fill+confAtoms[3], pdbM+fill+ogAtoms[3],
                             confName+fill+confAtoms[4], pdbM+fill+ogAtoms[4],
                             confName+fill+confAtoms[5], pdbM+fill+ogAtoms[5])

            print ''

            # generate new pdb file with orginal ligand replaced
            cmd.remove('ogLig')
            cmd.set('pdb_use_ter_records', 0)
            cmd.save(pdm + '_' + lig + 'rdm_P' + str(j) + '_C' + str(i) + '.pdb')


# execute
print 'Starting'
print ''
inFile = sys.argv[-1]
pdb, lig, atomNum, comboNum = getVars(inFile)
ogLigAtomList = proteinLigandAtoms(pdb)
alignLigs(pdb, lig, atomNum, comboNum, ogLigAtomList)
print 'Done'                   
            
