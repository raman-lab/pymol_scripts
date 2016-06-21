# maePairFit

import os
import glob

ligName = 'Piper_Conf'
originalLig = 'Piper'
base = '///UNK`900/'
atoms = ['C1', 'C3', 'C5']
cmd.reinitialize()
cmd.load(originalLig + '.pdb')

fileList = glob.glob(ligName + '*.pdb')
print fileList
for i in range(1, len(fileList)):
    namepdb = fileList[i]
    cmd.load(fileList[i])
    nameT = namepdb.split('.')[0]

    print nameT
    cmd.pair_fit(nameT+base+atoms[0], originalLig+base+atoms[0],
                 nameT+base+atoms[1], originalLig+base+atoms[1],
                 nameT+base+atoms[2], originalLig+base+atoms[2])
    print ''
