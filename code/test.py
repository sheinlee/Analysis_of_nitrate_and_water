from ccdc import io 
import os
from ccdc.io import MoleculeReader,CrystalReader
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time
import math

mol_reader=MoleculeReader('csd')
mol_name=mol_reader.molecule('AQAKIN')
ELEMENT='La'

#calculate distance between two atoms based on their coordinates
def calculate_the_distance(atom1_coordinates,atom2_coordinates):
    distance = math.sqrt((atom1_coordinates.x-atom2_coordinates.x)**2+(atom1_coordinates.y-atom2_coordinates.y)**2+(atom1_coordinates.z-atom2_coordinates.z)**2)
    return distance

find_element = False
for com in mol_name.components:
    for atom in com.atoms:
        if atom.atomic_symbol == ELEMENT: 
            find_element = True
            counted_nitrate_groups = set()
            ln_coord = atom.coordinates
            # print(atom.atomic_symbol)
            # print(atom.coordinates)
            # for atom1 in atom.neighbours:
                # print(atom1.atomic_symbol)
                # print(atom1.coordinates)
            for bond in com.bonds:
                if atom in bond.atoms:
                    print(bond.ideal_bond_length)
                    print(bond.length)

            #identify nitrate
            for atom1 in atom.neighbours:
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
                                    # bidentate_nitrate += 1
                                    print('bidentate')
                                    for atom3 in atom2.neighbours:
                                        if (atom3.atomic_symbol == 'O'):
                                            print(calculate_the_distance(atom3.coordinates,atom.coordinates))

                                elif count_donor_o == 1:
                                    # monodentate_nitrate += 1
                                    print('monodentate')
                                    for atom3 in atom2.neighbours:
                                        if (atom3.atomic_symbol == 'O'):
                                            print(calculate_the_distance(atom3.coordinates,atom.coordinates))

                                counted_nitrate_groups.add(nitrate_id)

        if find_element:break
    if find_element: break


print(len(mol_name.atoms))
# for com in mol_name.components:
#     print(len(com.atoms))
#     print(com.atoms)
#     average_distance=0
#     for atom in com.atoms:
#         average_distance+=calculate_the_distance(atom.coordinates,ln_coord)
#     average_distance=average_distance/len(com.atoms)
#     print(average_distance)

find_element = False
for com in mol_name.components:
    for atom in com.atoms:
        if atom.atomic_symbol == ELEMENT: 
            find_element = True
            print(com.atoms)
            print('com.formal_charge:',com.formal_charge)

            # for atom1 in atom.neighbours:
            #     print(atom1.atomic_symbol)
            #     print(atom1.partial_charge)
            #     print(atom1.formal_charge)
            # fs_net_charge1 = sum(neighbour.partial_charge for neighbour in atom.neighbours)
            # fs_net_charge2 = sum(neighbour.formal_charge for neighbour in atom.neighbours)
        if find_element:break
    if find_element: break

# print(fs_net_charge1)
# print(fs_net_charge2)
                
