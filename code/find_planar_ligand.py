#this script is used to find Ln complexes with planar ligand in csd

#we can evaluate how planar a ligand is
import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from time import process_time


def find_planar_ligand(ELEMENT,filepath,B,writerB):
    mol_reader = MoleculeReader(filepath, format='identifiers')
    for mol in mol_reader:
        try:
            mol.assign_bond_types()
            mol.add_hydrogens()
        except RuntimeError as e:
            print(f"Error adding hydrogens to molecule {mol.identifier}: {e}")
            continue

        find_ELEMENT = False
        # mol_is_organic = False
        for com in mol.components:
            for atom in com.atoms:
                if atom.atomic_symbol == ELEMENT:
                    find_ELEMENT = True
                    
                if find_ELEMENT: break
            if find_ELEMENT: break

def main():
    start = process_time()
    LIST_OF_ELEMENT = ['La',
                        'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
                        ]
    folder_path_list = ['planar_ligand']
    for folder_path in folder_path_list:
        # if filefolder does not exist, then create it
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    writerB = pd.ExcelWriter('analysis_of_planar_ligand_in_fs.xlsx')
    # A = np.zeros(shape=(4000, 10), dtype=object)
    B = np.zeros(shape=(10, 5), dtype=object)  # Updated shape to hold average and std

    
    for ELEMENT in LIST_OF_ELEMENT:
        filepath = './subset2_2023/subset2_' + ELEMENT + '.txt'
        file_path_list = ['planar_ligand/pl_'+ELEMENT+'.txt']
        for file_path in file_path_list:
            # if file has same name with other file ,then delete it
            if os.path.exists(file_path):
                os.remove(file_path)

        find_planar_ligand(ELEMENT,filepath,B,writerB)

    
    writerB.save()
    #processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))

if __name__ == "__main__":
    main()