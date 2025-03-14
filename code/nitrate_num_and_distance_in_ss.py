#this script is used to find how many and how far those nitrates in second shell
#file path =  data/subset3/merged_nitrate_hits.txt


from ccdc.io import MoleculeReader
import math
import pandas as pd
import os
from time import process_time

def is_nitrate(atomN):
    flag = False
    if all(atom.atomic_symbol == 'O' for atom in atomN.neighbours) and len(atomN.neighbours) == 3:
        flag = True
    return flag

def distance_calculation(ln_coordinates_list, N_coordinates_list):
    distance_list = []
    for lncrd in ln_coordinates_list:
        ln_x = lncrd.x
        ln_y = lncrd.y
        ln_z = lncrd.z
        for Ncrd in N_coordinates_list:
            N_x = Ncrd.x
            N_y = Ncrd.y
            N_z = Ncrd.z
            distance = math.sqrt((ln_x - N_x)**2 + (ln_y - N_y)**2 + (ln_z - N_z)**2)

            #condition distance between 4-15
            if distance<=15:
                distance_list.append(distance)
    return distance_list

def nitrate_num_and_distance_in_ss(ELEMENT):
    filepath = 'data/subset3/merged_nitrate_hits.txt'
    mol_reader = MoleculeReader(filepath, format='identifiers')
    identifiers = []
    lengths = [] 
    all_dis_lists = []

    #record dist_list
    list_dis_lists = []
    
    for mol in mol_reader:
        # try:
        #     mol.assign_bond_types()
        #     mol.add_hydrogens()
        # except RuntimeError as e:
        #     print(f"Error adding hydrogens to molecule {mol.identifier}: {e}")
        #     continue
        
        lncrds = []
        Ncrds = []
        has_element = False
        for atom in mol.atoms:
            if atom.atomic_symbol == ELEMENT:
                has_element = True
                lncrd = atom.coordinates
                lncrds.append(lncrd)
        for com in mol.components:
            if any(atom.atomic_symbol == ELEMENT for atom in com.atoms): continue 
            if has_element:
                for atom in com.atoms:
                    if atom.atomic_symbol == 'N':
                        if is_nitrate(atom):
                            Ncrd = atom.coordinates
                            Ncrds.append(Ncrd)
        dis_list = distance_calculation(lncrds, Ncrds)
        if len(dis_list) > 0:
            # identifiers.extend([mol.identifier] * len(dis_list))
            identifiers.append(mol.identifier)
            lengths.append(len(dis_list))
            all_dis_lists.extend(dis_list)
            list_dis_lists.append(dis_list)
    
    return identifiers, lengths, all_dis_lists,list_dis_lists

def main():
    start = process_time()
    LIST_OF_ELEMENT = [
        'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
    ]
    
    folder_path = 'data/subset3'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    all_identifiers = []
    all_lengths = []
    all_dis_lists = []
    all_list_dis_lists = []
    
    for ELEMENT in LIST_OF_ELEMENT:
        identifiers, lengths, dis_list ,list_dis_list= nitrate_num_and_distance_in_ss(ELEMENT)
        all_identifiers.extend(identifiers)
        all_lengths.extend(lengths)
        all_dis_lists.extend(dis_list)
        all_list_dis_lists.extend(list_dis_list)
    
    #test
    print(len(all_identifiers))
    print(len(all_dis_lists))
    print(len(all_lengths))
    print(len(all_list_dis_lists))
    # Create DataFrames for identifiers, lengths, and distances
    lengths_df = pd.DataFrame({
        'Identifier': all_identifiers,
        'Length': all_lengths,
        'dis_list': all_list_dis_lists
    })
    
    dis_list_df = pd.DataFrame({
        # 'Identifier': all_identifiers,
        'Distance': all_dis_lists
    })
    
    # Save to Excel
    lengths_path = os.path.join(folder_path, 'lengths15.xlsx')
    dis_list_path = os.path.join(folder_path, 'distances15.xlsx')
    
    lengths_df.to_excel(lengths_path, index=False)
    dis_list_df.to_excel(dis_list_path, index=False)
    
    # Processing time
    end = process_time()
    print('Running time: %s Seconds' % (end - start))

if __name__ == "__main__":
    main()




    # countc = 0
    # cry_reader = CrystalReader(filepath, format='identifiers')
    # for cry in cry_reader:
    #     asy_cry_name = cry.asymmetric_unit_molecule
    #     lncrds = []
    #     Ncrds = []
    #     has_element = False
    #     for atom in asy_cry_name.atoms:
    #             if atom.atomic_symbol == ELEMENT:
    #                 has_element = True
    #                 lncrd = atom.coordinates
    #                 lncrds.append(lncrd)
    #     for com in asy_cry_name.components:
    #         if any(atom.atomic_symbol == ELEMENT for atom in com.atoms): continue 
    #         if has_element:
    #             for atom in com.atoms:
    #                 if atom.atomic_symbol == 'N':
                        
    #                     if is_nitrate(atom):
    #                         Ncrd = atom.coordinates #Ncrd.x, Ncrd.y, Ncrd.z
    #                         Ncrds.append(Ncrd)
    #     dis_list = distance_calculation(lncrds,Ncrds)
    #     if len(dis_list)>0:
    #         print(cry.identifier) 
    #         print(dis_list)
    #         for com in asy_cry_name.components:
    #             print(com.atoms)

    #     if countc ==3: break
    #     countc += 1   