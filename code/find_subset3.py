# this script is used for finding subset3. 
# Subset3 : no nitrate in first shell.


from ccdc import io 
from ccdc.io import MoleculeReader, CrystalReader,EntryReader
import numpy as np
import pandas as pd
from itertools import compress
import os
from time import process_time

def Count_Ln(molecule,element):
    countln=0
    for atom in molecule.atoms:
        if atom.atomic_symbol==element:
            countln+=1
    return countln

def Count_Metal(molecule):
    countmetal=0
    for atom in molecule.atoms:
        if atom.is_metal:
            countmetal+=1
    return countmetal

def is_mono(ELEMENT,cry_name):
    # ELEMENT='Dy'
    countb_is_mono=False
    cry_reader=CrystalReader('csd')

    cry_name=cry_reader.crystal(cry_name)


    asy_cry_name=cry_name.asymmetric_unit_molecule
    for com in asy_cry_name.components:
        print(com.atoms)


    if len(asy_cry_name.components)==1:
        if Count_Ln(asy_cry_name,ELEMENT)==1 and Count_Metal(asy_cry_name) == 1:
            countb_is_mono=True
    else:
        
        for com1 in asy_cry_name.components:
            if Count_Metal(com1) == 1 and Count_Ln(com1,ELEMENT) == 1:
                countb_is_mono=True
                



    if countb_is_mono: 
        print('this cry is mono!')
    else:
        print('this cry is not mono!') 

    return countb_is_mono

def read_file(filepath):
    with open(filepath, 'r') as file:
        return set(file.readlines())

def write_file(filepath, lines):
    with open(filepath, 'w') as file:
        file.writelines(lines)




def find_subset3(ELEMENT):
    filepath1 = 'data/subset2_2023/subset2_' + ELEMENT + '.txt'    
    filepath2 = 'data/nitrate/nitrate_' + ELEMENT + '.txt'
    filepath3 = f'data/subset3/subset3_{ELEMENT}.txt'
    filepath4 = f'data/subset3/subset3_mono_{ELEMENT}.txt'

    #if filepath3/4 exists, delete it.
    if os.path.exists(filepath3):
        os.remove(filepath3)
    if os.path.exists(filepath4):
        os.remove(filepath4)

    def extract_unique_elements():
        
        # 读取文件内容
        set1 = read_file(filepath1)
        set2 = read_file(filepath2)

        # 找出文件1中存在但文件2中不存在的部分
        unique_lines = set1 - set2

        # 将这些部分写入文件3
        write_file(filepath3, unique_lines)

    extract_unique_elements()

    mol_reader = MoleculeReader(filepath3, format='identifiers')
    for mol in mol_reader:

        #verify is mono or not
        if is_mono(ELEMENT,mol.identifier):
            with open(filepath4, 'a') as file:
                file.write(mol.identifier+ '\n')

    mol_reader.close()
    return 0

def main():
    start = process_time()
    LIST_OF_ELEMENT = [
        'La', 'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
                    ]       
    
    folder_path = 'data/subset3'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path) 
        
    for ELEMENT in LIST_OF_ELEMENT:

        find_subset3(ELEMENT)    

    #processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))

if __name__ == "__main__":
    main()