from ccdc import io 
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from time import process_time
start=process_time()

LIST_OF_ELEMENT = ['La',
                   'Ce',
                   'Pr',
                   'Nd',
                   'Sm',
                   'Eu',
                   'Gd',
                   'Tb',
                   'Dy',
                   'Ho',
                   'Er',
                   'Tm',
                   'Yb',
                   'Lu'
                   ]



countno3d = 0
mol_reader = MoleculeReader('csd')
for mol in mol_reader:

    mol_name=mol
    molecule_name=mol_name.identifier

    #identify whether molecule has 3D info 
    # print(mol_name.is_3d)
    if mol_name.is_3d == False: 
        countno3d += 1 
        # print(countno3d)
        with open("no3d.txt","a") as f:
                f.write(molecule_name)
                f.write("\n")
        continue
    
    print(mol_name.smiles)
    # print(mol_name.atoms)

    for ELEMENT in LIST_OF_ELEMENT:


        max_cn_for_mol=0
        for atom in mol_name.atoms:
            if atom.atomic_symbol==ELEMENT:
                if len(atom.neighbours)>max_cn_for_mol:
                    max_cn_for_mol=len(atom.neighbours)

        if max_cn_for_mol>2:
            with open('subset1_'+ELEMENT+'.txt','a') as f:
                        f.write(mol_name.identifier + '\n')

mol_reader.close()

#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))