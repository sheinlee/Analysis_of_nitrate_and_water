#find crystal with nitrate in subset3. in most cases, the nitrate will be at second shell
 
import ccdc.search
import ccdc.io
from ccdc.io import MoleculeReader, CrystalReader,EntryReader
import numpy as np
import pandas as pd
from itertools import compress
import os
from time import process_time

start = process_time()
LIST_OF_ELEMENT = [
'La', 
'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
            ]       

entry_reader = ccdc.io.EntryReader('CSD')
ADIDUQ01 = entry_reader.molecule('ADIDUQ01')
# for com in ADIDUQ01.components:
#     print(com.atoms)
nitrate = ADIDUQ01.components[7]
print(nitrate.atoms)
for ELEMENT in LIST_OF_ELEMENT:
    filepath4 = f'data/subset3/subset3_mono_{ELEMENT}.txt'
    # filepath4 = f'data/nitrate/nitrate_{ELEMENT}.txt'
    nitrate_substructure = ccdc.search.MoleculeSubstructure(nitrate)
    nitrate_search = ccdc.search.SubstructureSearch()
    sub_id = nitrate_search.add_substructure(nitrate_substructure)
    nitrate_hits = nitrate_search.search(ccdc.io.MoleculeReader(filepath4, format='identifiers'))
    print(len(nitrate_hits))
    print(nitrate_hits)
    for hit in nitrate_hits:

        with open(f'data/subset3/nitrate_hits_{ELEMENT}.txt','a') as f:
            f.write(hit.identifier+'\n')
        

#processing time
end = process_time()
print('Running time: %s Seconds' % (end - start))