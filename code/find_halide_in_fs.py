#find halide with at least 1 organic ligand in first shell in csd

import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from time import process_time

start = process_time()
LIST_OF_ELEMENT = ['La',
                    'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
                    ]

halide_list = ['F','Cl','Br','I']
folder_path_list = ['org_halide','error_adding_H']
for folder_path in folder_path_list:
    # 如果文件夹不存在，创建它
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# writerA = pd.ExcelWriter('org_halide_in_fs.xlsx')
# writerB = pd.ExcelWriter('analysis_org_halide_in_fs.xlsx')
# A = np.zeros(shape=(4000, 10), dtype=object)
# B = np.zeros(shape=(10, 5), dtype=object)  # Updated shape to hold average and std

for ELEMENT in LIST_OF_ELEMENT:


    file_path_list = ['org_halide/org_halide_'+ELEMENT+'.txt','error_adding_H/error_adding_H_'+ELEMENT+'.txt']
    for file_path in file_path_list:
        # 如果同名文件存在，删除它
        if os.path.exists(file_path):
            os.remove(file_path)
    count = -1
    filepath = './subset2_2023/subset2_' + ELEMENT + '.txt'
    mol_reader = MoleculeReader(filepath, format='identifiers')

    CN_values_halide = []
    count_halide=0
    for mol in mol_reader:
        count += 1
        # A[count, 0] = mol.identifier

        try:
            mol.assign_bond_types()
            mol.add_hydrogens()
        except RuntimeError as e:
            print(f"Error adding hydrogens to molecule {mol.identifier}: {e}")

            with open('error_adding_H/error_adding_H_'+ELEMENT+'.txt','a') as f:
                f.write(mol.identifier + '\n')
            continue

        find_halide = False
        find_ELEMENT = False
        mol_is_organic = False
        for com in mol.components:
            for atom in com.atoms:
                if atom.atomic_symbol == ELEMENT:
                    find_ELEMENT = True
                    counted_nitrate_groups = set()
                    CN = len(atom.neighbours)
                    for atom1 in atom.neighbours:
                        if atom1.atomic_symbol in halide_list:
                            
                            find_halide = True

                   #identify is_organic
                    for atom1 in com.atoms:
                        if atom1.atomic_symbol == 'C':
                            mol_is_organic = True                
                if find_ELEMENT: break
            if find_ELEMENT: break
    if mol_is_organic and find_halide:
        count_halide+=1
        CN_values_halide.append(CN)
        with open('org_halide/org_halide_'+ELEMENT+'.txt','a') as f:
            f.write(mol.identifier + '\n')

    mol_reader.close()

end = process_time()
print('Running time: %s Seconds' % (end - start))