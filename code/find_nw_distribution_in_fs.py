#this script is used to find nitrate distribution in subset2n and water distribution in subset2w

import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from collections import Counter
from time import process_time

def match_atomic_symbols(neighbours, target_list):
    neighbour_symbols = [neighbour.atomic_symbol for neighbour in neighbours]
    target_count = {symbol: target_list.count(symbol) for symbol in target_list}
    neighbour_count = {symbol: neighbour_symbols.count(symbol) for symbol in set(neighbour_symbols)}
    return target_count == neighbour_count

def find_distribution_nitrate_water_infs(ELEMENT, key, A, writerA, data_dict):
    filepath = f'data/{key}/{key}_' + ELEMENT + '.txt'
    count = -1
    list_oh = [ELEMENT, 'H']
    list_h2o = [ELEMENT, 'H', 'H']

    mol_reader = MoleculeReader(filepath, format='identifiers')

    num_nitrate = []
    num_water = []

    for mol in mol_reader:
        count += 1

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
                                        if atom3.atomic_symbol == 'O':
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
                    # identify is_organic
                    for atom1 in com.atoms:
                        if atom1.atomic_symbol == 'C':
                            mol_is_organic = True
                if find_ELEMENT:
                    break
            if find_ELEMENT:
                break

        num_nitrate.append(monodentate_nitrate + bidentate_nitrate)
        num_water.append(h2o + oh)

    if key == 'nitrate':
        nn_counter = Counter(num_nitrate)
        for coord_num, freq in nn_counter.items():
            if coord_num in data_dict[key]:
                data_dict[key][coord_num][ELEMENT] = freq
            else:
                data_dict[key][coord_num] = {ELEMENT: freq}

    if key == 'water':
        nw_counter = Counter(num_water)
        for coord_num, freq in nw_counter.items():
            if coord_num in data_dict[key]:
                data_dict[key][coord_num][ELEMENT] = freq
            else:
                data_dict[key][coord_num] = {ELEMENT: freq}

def main():
    start = process_time()
    LIST_OF_ELEMENT = [
        'La', 'Ce',
          'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
    ]

    writerA = pd.ExcelWriter('distribution_nw_in_fs.xlsx')

    A = np.zeros(shape=(100, 100), dtype=object)
    key_list = ['nitrate', 'water']

    data_dict = {'nitrate': {}, 'water': {}}

    for ELEMENT in LIST_OF_ELEMENT:
        for key in key_list:
            find_distribution_nitrate_water_infs(ELEMENT, key, A, writerA, data_dict)

    for key in key_list:
        combined_data = []
        for coord_num, elements_freq in data_dict[key].items():
            row = {'Coordination Number': coord_num}
            row.update(elements_freq)
            combined_data.append(row)
        df = pd.DataFrame(combined_data)
        df.to_excel(writerA, sheet_name=key, index=False)

    writerA.save()

    # processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))

if __name__ == "__main__":
    main()

