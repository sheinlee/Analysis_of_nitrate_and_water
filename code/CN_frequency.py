import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from collections import Counter
from time import process_time


def find_cn_frequency(ELEMENT, key_word,writerA):
    file_path = f'data/{key_word}/{key_word}_{ELEMENT}.txt'
    mol_reader = MoleculeReader(file_path, format='identifiers')
    CN_values = []
    count = 0
    for mol in mol_reader:
        count += 1
        # A[count, 0] = mol.identifier

        try:
            mol.assign_bond_types()
            mol.add_hydrogens()
        except RuntimeError as e:
            print(f"Error adding hydrogens to molecule {mol.identifier}: {e}")
            continue
        find_ELEMENT = False
        for com in mol.components:
            for atom in com.atoms:
                if atom.atomic_symbol == ELEMENT:
                    find_ELEMENT = True
                    counted_nitrate_groups = set()
                    CN = len(atom.neighbours)
                    CN_values.append(CN)
                if find_ELEMENT: break
            if find_ELEMENT: break
    mol_reader.close()
    cn_counter = Counter(CN_values)

    # Prepare data for the DataFrame
    freq_data = {'Coordination Number': list(cn_counter.keys()), 'Frequency': list(cn_counter.values())}
    freq_df = pd.DataFrame(freq_data)

    # Write the frequency data to the Excel file
    freq_df.to_excel(writerA, sheet_name=ELEMENT, index=False)

    # Write the additional data A to the Excel file
    # pd.DataFrame(A).to_excel(writerA, sheet_name=f'{ELEMENT}_Data', float_format='%.5f')

    return cn_counter

            
def main():
    start = process_time()
    LIST_OF_ELEMENT = [
        'La', 
        # 'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
    ]
    # folder_path_list = ['water', 'nitrate', 'org_water', 'org_nitrate', 'org_nitrate_no_water', 'org_water_and_nitrate', 'org_water_no_nitrate']
    key_word_list = ['water', 'nitrate', 'org_water', 'org_nitrate', 'org_nitrate_no_water', 'org_water_and_nitrate', 'org_water_no_nitrate']

    # for folder_path in folder_path_list:
    #     # 如果文件夹不存在，创建它
    #     if not os.path.exists(folder_path):
    #         os.makedirs(folder_path)

    for key_word in key_word_list:
        writerA = pd.ExcelWriter(f'analysis_CN_frequency_{key_word}.xlsx')
        for ELEMENT in LIST_OF_ELEMENT:
            
            A = np.zeros(shape=(30, 4), dtype=object)
            cn_frequency = find_cn_frequency(ELEMENT, key_word,writerA)
        writerA.save()

    # processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))


if __name__ == "__main__":
    main()
