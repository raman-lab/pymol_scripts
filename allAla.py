# replace LG into all Alanine Protein
import os
import glob
import sys

cmd.reinitialize()
ALApdb = sys.argv[-1]


pdbList = glob.glob('*_P*_C*_4AC0.pdb')
pdb = ', '.join(pdbList)

for i,entry in enumerate(pdbList):
    ALApdbS = ALApdb.split('.')[0]
    cmd.load(ALApdb)
    cmd.remove('organic')

    pName = entry.split('.')[0]
    print str(i) + ' ' + pName
    cmd.load(pName + '.pdb')
    cmd.create('LG1', 'organic')
    cmd.remove(pName)
    cmd.set('pdb_use_ter_records', 0)
    cmd.save('ALA_'+pName+'_ALL.pdb')
    cmd.reinitialize()
