from ccdc import io 
from ccdc.io import MoleculeReader, CrystalReader, EntryReader
import numpy as np
import pandas as pd
from itertools import compress
import re
from time import process_time

start = process_time()

def Count_Ln(molecule, element):
    countln = 0
    for atom in molecule.atoms:
        if atom.atomic_symbol == element:
            countln += 1
    return countln

def Count_Metal(molecule):
    countmetal = 0
    for atom in molecule.atoms:
        if atom.is_metal:
            countmetal += 1
    return countmetal

LIST_OF_ELEMENT = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
cry_reader = CrystalReader('csd')

# Read refcodes from file
with open('data/search3.txt', 'r') as file:
    refcodes = file.read().split()

results = []

for refcode in refcodes:
    countb_is_mono = False
    try:
        cry_name = cry_reader.crystal(refcode)
        asy_cry_name = cry_name.asymmetric_unit_molecule
        
        if len(asy_cry_name.components) == 1:
            if any(Count_Ln(asy_cry_name, element) == 1 for element in LIST_OF_ELEMENT) and Count_Metal(asy_cry_name) == 1:
                countb_is_mono = True
        else:
            for com1 in asy_cry_name.components:
                if Count_Metal(com1) == 1 and any(Count_Ln(com1, element) == 1 for element in LIST_OF_ELEMENT):
                    countb_is_mono = True
        
        results.append((refcode, countb_is_mono, cry_name.has_disorder))
        
        if countb_is_mono: 
            print(f'{refcode}: this cry is mono!')
        else:
            print(f'{refcode}: this cry is not mono!')
        
    except Exception as e:
        print(f'{refcode}: Error processing crystal - {e}')
        results.append((refcode, 'Error', str(e)))

# Print or save the results as needed
results_df = pd.DataFrame(results, columns=['Refcode', 'IsMono', 'HasDisorder'])
print(results_df)

# Optionally, save the results to a CSV file
results_df.to_csv('mononuclear_results.csv', index=False)

print(f'Process completed in {process_time() - start:.2f} seconds.')
