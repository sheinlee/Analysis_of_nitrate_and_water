from ccdc.io import MoleculeReader
import pandas as pd
from time import process_time

start = process_time()

LIST_OF_ELEMENTS = [
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
]

counts = []

for element in LIST_OF_ELEMENTS:
    filepath = f'subset2_{element}.txt'
    with MoleculeReader(filepath, format='identifiers') as mol_reader:
        counts.append(len(mol_reader))

df = pd.DataFrame([counts], columns=LIST_OF_ELEMENTS)
df.to_excel('num_in_subset2.xlsx', sheet_name='subset2', float_format='%.5f', index=False)

# Record processing time
end = process_time()
print(f'Running time: {end - start} Seconds')
