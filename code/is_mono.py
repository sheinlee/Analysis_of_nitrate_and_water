#test crystal is mono or not
from ccdc import io 
from ccdc.io import MoleculeReader, CrystalReader,EntryReader
import numpy as np
import pandas as pd
from itertools import compress
import re
from time import process_time


start=process_time()

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

# mol_reader=MoleculeReader('csd')
# mol_name=mol_reader.molecule('AQAKIN')


ELEMENT='Dy'
countb_is_mono=False
cry_reader=CrystalReader('csd')
csd_reader=io.EntryReader('csd')
cry_name=cry_reader.crystal('UPEMOU')
# entry_name=csd_reader.entry(cry_name)
# if entry_name.has_3d_structure : print('no 3D structure')

asy_cry_name=cry_name.asymmetric_unit_molecule
for com in asy_cry_name.components:
    print(com.atoms)


if len(asy_cry_name.components)==1:
    if Count_Ln(asy_cry_name,ELEMENT)==1 and Count_Metal(asy_cry_name) == 1:
        countb_is_mono=True
else:
    for com1 in asy_cry_name.components:
        if Count_Metal(com1)==1 and Count_Ln(com1,ELEMENT)==1:
            countb_is_mono=True



if countb_is_mono: 
    print('this cry is mono!')
else:
    print('this cry is not mono!') 

print(cry_name.disorder)
print(cry_name.disordered_molecule)   
print(cry_name.has_disorder)