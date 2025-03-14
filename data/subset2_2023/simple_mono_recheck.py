#this script is the final version to get subset2 (mononuclear subset)

from ccdc.io import CrystalReader, EntryReader
import re
from time import process_time
start=process_time()
def count_center_metal(molecule, element):
    return sum(1 for atom in molecule.atoms if atom.atomic_symbol == element)

def count_metal(molecule):
    return sum(1 for atom in molecule.atoms if atom.is_metal)



LIST_OF_ELEMENTS = [
    'La', 
    'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
]

Ln_ELEMENTS = [
    'La', 
    'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
]

def count_ln(molecule):
    return sum(1 for atom in molecule.atoms if atom.atomic_symbol in Ln_ELEMENTS)

def process_elements(elements):
    for element in elements:
        element_pattern = re.compile(rf'{element}(\d+)')
        filepath = f'../subset1/subset1_{element}.txt'

        with CrystalReader(filepath, format='identifiers') as cry_reader:
            for cry in cry_reader:
                cry_name = cry_reader.crystal(cry.identifier)
                entry_name = EntryReader('CSD').entry(cry.identifier)

                if not entry_name.has_3d_structure:
                    continue

                asy_cry_name = cry_name.asymmetric_unit_molecule
                if any(count_center_metal(com, element)+count_ln(com) == 2 for com in asy_cry_name.components):
                    formula = entry_name.formula
                    match = element_pattern.search(formula)
                    if match and int(match.group(1)) == 1:
                        with open(f'subset2_{element}.txt', 'a') as f:
                            f.write(cry_name.identifier + '\n')

process_elements(LIST_OF_ELEMENTS)

end=process_time()
print('Running time: %s Seconds'%(end-start))