#find ratio of nitrate/water in first coordination shell in subset2n_w

import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from time import process_time


def match_atomic_symbols(neighbours, target_list):
    neighbour_symbols = [neighbour.atomic_symbol for neighbour in neighbours]
    target_count = {symbol: target_list.count(symbol) for symbol in target_list}
    neighbour_count = {symbol: neighbour_symbols.count(symbol) for symbol in set(neighbour_symbols)}
    return target_count == neighbour_count

    
def find_nitrate_and_water(ELEMENT,filepath,A,B,writerB):

    count = -1
    list_oh = [ELEMENT, 'H']
    list_h2o = [ELEMENT, 'H', 'H']

    mol_reader = MoleculeReader(filepath, format='identifiers')

    CN_values = []

    ratio_values = []

    for mol in mol_reader:
        count += 1
        A[count, 0] = mol.identifier

        try:
            mol.assign_bond_types()
            mol.add_hydrogens()
        except RuntimeError as e:
            print(f"Error adding hydrogens to molecule {mol.identifier}: {e}")
            continue

        monodentate_nitrate = 0
        bidentate_nitrate = 0
        h2o = 0
        oh = 0
        find_ELEMENT = False
        mol_is_organic = False
        for com in mol.components:
            for atom in com.atoms:
                if atom.atomic_symbol == ELEMENT:
                    find_ELEMENT = True
                    counted_nitrate_groups = set()
                    CN = len(atom.neighbours)
                    CN_values.append(CN)
                    for atom1 in atom.neighbours:
                        if atom1.atomic_symbol == 'O':
                            if match_atomic_symbols(atom1.neighbours, list_oh):
                                oh += 1
                            elif match_atomic_symbols(atom1.neighbours, list_h2o):
                                h2o += 1
                                
                            for atom2 in atom1.neighbours:
                                if atom2.atomic_symbol == 'N':  # atom2 is N
                                    count_o = 0
                                    count_donor_o = 0
                                    oxygen_ids = []
                                    for atom3 in atom2.neighbours:
                                        if (atom3.atomic_symbol == 'O'):
                                            count_o += 1
                                            oxygen_ids.append(atom3.index)
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol == ELEMENT:
                                                    count_donor_o += 1
                                    if count_o == 3:
                                        nitrate_id = tuple(sorted(oxygen_ids + [atom2.index]))
                                        if nitrate_id not in counted_nitrate_groups:
                                            if count_donor_o == 2:
                                                bidentate_nitrate += 1
                                            elif count_donor_o == 1:
                                                monodentate_nitrate += 1
                                            counted_nitrate_groups.add(nitrate_id)
                    #identify is_organic
                    for atom1 in com.atoms:
                        if atom1.atomic_symbol == 'C':
                            mol_is_organic = True                
                if find_ELEMENT: break
            if find_ELEMENT: break


        A[count,1] = CN                                     #CN distribution
        A[count,2] = h2o                                    #Number of water distribution
        A[count,3] = bidentate_nitrate+monodentate_nitrate  #if A[count,3] == 3, first shell is neutral
        A[count,4] = bidentate_nitrate 
        A[count,5] = monodentate_nitrate   

        ratio_values.append((bidentate_nitrate+monodentate_nitrate)/h2o)
    mol_reader.close()

    B[0,0]= np.mean(ratio_values) 
    B[0,1] = np.std(ratio_values)  

    
    # pd.DataFrame(A).to_excel(writerA, sheet_name=ELEMENT, float_format='%.5f')
    pd.DataFrame(B).to_excel(writerB, sheet_name=ELEMENT, float_format='%.5f')       
    

def main():
    start = process_time()
    LIST_OF_ELEMENT = [
        'La', 'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
                        ]
    folder_path_list = ['water','nitrate','org_water','org_nitrate','org_nitrate_no_water','org_water_and_nitrate','org_water_no_nitrate']
    for folder_path in folder_path_list:
        # 如果文件夹不存在，创建它
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    # writerA = pd.ExcelWriter('org_n_and_w_in_fs.xlsx')
    writerB = pd.ExcelWriter('ratio_org_n_and_w_in_fs.xlsx')
    A = np.zeros(shape=(4000, 10), dtype=object)
    B = np.zeros(shape=(10, 5), dtype=object)  # Updated shape to hold average and std

    
    for ELEMENT in LIST_OF_ELEMENT:
        filepath = './org_water_and_nitrate/org_water_and_nitrate_'+ELEMENT+'.txt'

        find_nitrate_and_water(ELEMENT,filepath,A,B,writerB)

    # writerA.save()
    writerB.save()
    #processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))


if __name__ == "__main__":
    main()