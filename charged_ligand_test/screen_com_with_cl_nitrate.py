#this script is used to screen out complex with at least one charged ligand in first coordination shell.
#target: subset 2; nitrate; water; subset3(either nitrate or water)
#step 1: remove ln from ln mono complex
#step 2: create new_mol to store the ligand
#step 3: obtain the total charge for each ligand
#step 4: identify the organic ligand with at least one C_C bond
#step 5: distribution analysis for the ligand


#display average charge for organic ligand for each complex
#with/without organic ligand in complex

from ccdc.io import MoleculeReader
import re
from time import process_time
import os
from ccdc.molecule import Molecule
from ccdc import io
import pandas as pd
from collections import defaultdict
import math


def contains_carbon_carbon_bond(molecule):
    """
    判断分子中是否包含碳-碳键。
    """
    for bond in molecule.bonds:
        atoms = bond.atoms
        if atoms[0].atomic_symbol == 'C' and atoms[1].atomic_symbol == 'C':
            return True
    return False

def is_nitrate(atomN):
    flag = False
    if all(atom.atomic_symbol == 'O' for atom in atomN.neighbours) and len(atomN.neighbours) == 3:
        flag = True
    return flag

# def charge_of_ligand_calculation(mol, ELEMENT):
#     # 打开CSD数据库
#     mol_reader = io.MoleculeReader('csd')

#     # # 读取分子
#     # mol = mol_reader.molecule('ACUHAI')
#     # ELEMENT = 'La'  # 目标金属元素
        
#     # 获取最重的组分，通常是主要的配体
#     main_component = mol.heaviest_component

#     # 创建一个新的分子对象
#     new_mol = Molecule(identifier=mol.identifier)

#     # 记录需要移除的原子和键
#     atoms_to_remove = set()
#     bonds_to_remove = set()

#     # 遍历原始分子的所有原子
#     for atom in main_component.atoms:
#         if atom.atomic_symbol == ELEMENT:
#             # 记录与金属原子相连的键和原子
#             for bond in atom.bonds:
#                 bonds_to_remove.add(bond)
#                 # connected_atom = bond.atoms(atom)
#                 # atoms_to_remove.add(connected_atom)
#             # 记录金属原子本身
#             atoms_to_remove.add(atom)

#     # 从原始分子中移除记录的键和原子
#     for bond in bonds_to_remove:
#         main_component.remove_bond(bond)
#     for atom in atoms_to_remove:
#         if atom in main_component.atoms:  # 确保原子仍在分子中
#             main_component.remove_atom(atom)

#     # 将剩余的原子和键添加到新分子中
#     for atom in main_component.atoms:
#         new_mol.add_atom(atom)
#     for bond in main_component.bonds:
#         new_mol.add_bond(bond.bond_type, bond.atoms[0], bond.atoms[1])

#     contains_organic_ligand = False
#     sum_formal_charge_for_organic_ligand = 0
#     num_organic_ligand = 0
#     for com in new_mol.components:
#         # print(len(com.atoms))
#         # print(com.atoms)
#         if not com.atoms:
#             continue  # 跳过空的组分
        
#         try:
#             com.assign_bond_types()
#             com.set_formal_charges()
#         except RuntimeError:
#             print(f"⚠️ 跳过无法分配键类型的分子 {mol.identifier}")
#             continue  # 如果 assign_bond_types() 失败，则跳过该组分

#         if contains_carbon_carbon_bond(com):
#             contains_organic_ligand = True
#             sum_formal_charge_for_organic_ligand += com.formal_charge
#             num_organic_ligand += 1
#     if num_organic_ligand == 0:
#         average_formal_charge_for_organic_ligand = 0
#     else:
#         average_formal_charge_for_organic_ligand = sum_formal_charge_for_organic_ligand / num_organic_ligand

#     return average_formal_charge_for_organic_ligand, contains_organic_ligand



# def screen_com_with_cl(ELEMENT):
#     """ 处理单个金属元素的配体统计 """
#     filepath = f'../data/subset2_2023/subset2_{ELEMENT}.txt'
#     mol_reader = MoleculeReader(filepath, format='identifiers')

#     count_with_organic = 0
#     count_without_organic = 0
#     organic_ligand_charges = []

#     for mol in mol_reader:
#         try:
            
#             charge, is_org = charge_of_ligand_calculation(mol, ELEMENT)

#             if is_org:
#                 count_with_organic += 1
#                 organic_ligand_charges.append(charge)
#             else:
#                 count_without_organic += 1

#         except RuntimeError:
#             print(f"⚠️ 跳过无效的 CSD 代码: {mol.identifier}")
#             continue  # 跳过该分子，继续下一个

#     return count_with_organic, count_without_organic, organic_ligand_charges

def charge_of_ligand_calculation(mol, ELEMENT):
    # 打开CSD数据库
    mol_reader = io.MoleculeReader('csd')

    # 获取最重的组分，通常是主要的配体
    main_component = mol.heaviest_component

    # 创建一个新的分子对象
    new_mol = Molecule(identifier=mol.identifier)

    # 记录需要移除的原子和键
    atoms_to_remove = set()
    bonds_to_remove = set()

    # 遍历原始分子的所有原子
    for atom in main_component.atoms:
        if atom.atomic_symbol == ELEMENT:
            # 记录与金属原子相连的键和原子
            for bond in atom.bonds:
                bonds_to_remove.add(bond)
            # 记录金属原子本身
            atoms_to_remove.add(atom)

    # 从原始分子中移除记录的键和原子
    for bond in bonds_to_remove:
        main_component.remove_bond(bond)
    for atom in atoms_to_remove:
        if atom in main_component.atoms:  # 确保原子仍在分子中
            main_component.remove_atom(atom)

    # 将剩余的原子和键添加到新分子中
    for atom in main_component.atoms:
        new_mol.add_atom(atom)
    for bond in main_component.bonds:
        new_mol.add_bond(bond.bond_type, bond.atoms[0], bond.atoms[1])

    contains_organic_ligand = False
    organic_ligand_charges = []  # 存储所有有机配体的电荷

    for com in new_mol.components:
        if not com.atoms:
            continue  # 跳过空的组分
        
        try:
            com.assign_bond_types()
            com.set_formal_charges()
        except RuntimeError:
            print(f"⚠️ 跳过无法分配键类型的分子 {mol.identifier}")
            continue  # 如果 assign_bond_types() 失败，则跳过该组分

        if contains_carbon_carbon_bond(com):
            contains_organic_ligand = True
            organic_ligand_charges.append(com.formal_charge)  # 记录每个有机配体的电荷

    return organic_ligand_charges, contains_organic_ligand


def screen_com_with_cl(ELEMENT):
    """ 处理单个金属元素的配体统计 """
    filepath = f'../data/nitrate/nitrate_{ELEMENT}.txt'
    mol_reader = MoleculeReader(filepath, format='identifiers')

    count_with_organic = 0
    count_without_organic = 0
    organic_ligand_charges = []

    for mol in mol_reader:
        try:
            ligand_charges, is_org = charge_of_ligand_calculation(mol, ELEMENT)

            if is_org:
                count_with_organic += 1
                organic_ligand_charges.extend(ligand_charges)  # 收集所有有机配体的电荷
            else:
                count_without_organic += 1

        except RuntimeError:
            print(f"⚠️ 跳过无效的 CSD 代码: {mol.identifier}")
            continue  # 跳过该分子，继续下一个

    return count_with_organic, count_without_organic, organic_ligand_charges


from tqdm import tqdm  # 引入 tqdm 进度条库

def main():
    LIST_OF_ELEMENT = [
        'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
    ]

    results = defaultdict(dict)
    charge_distribution = []

    print("Processing nitrate...")
    for ELEMENT in tqdm(LIST_OF_ELEMENT, desc="Processing Elements", unit="element"):
        count_with_organic, count_without_organic, organic_ligand_charges = screen_com_with_cl(ELEMENT)

        results[ELEMENT]['With Organic Ligand'] = count_with_organic
        results[ELEMENT]['Without Organic Ligand'] = count_without_organic

        if organic_ligand_charges:
            results[ELEMENT]['Average Charge of Organic Ligands'] = sum(organic_ligand_charges) / len(organic_ligand_charges)

            # 统计电荷分布并存储到 charge_distribution 列表
            for charge in organic_ligand_charges:
                rounded_charge = round(charge)  # 取最接近的整数
                floor_charge = math.floor(charge)

                if -5 <= floor_charge <= 5:
                    charge_distribution.append([ELEMENT, floor_charge, 1])  # 每行为 [元素, 电荷, 计数]
        else:
            results[ELEMENT]['Average Charge of Organic Ligands'] = None

    # 转换为 Pandas DataFrame 并存储到 Excel
    df_results = pd.DataFrame.from_dict(results, orient='index')

    # 统计所有元素的电荷分布
    df_distribution = pd.DataFrame(charge_distribution, columns=['Element', 'Charge', 'Count'])
    df_distribution = df_distribution.groupby(['Element', 'Charge']).sum().reset_index()  # 按元素和电荷统计

    output_path = 'nitrate_analysis.xlsx'
    with pd.ExcelWriter(output_path) as writer:
        df_results.to_excel(writer, sheet_name='Ligand Analysis')
        df_distribution.to_excel(writer, sheet_name='Charge Distribution', index=False)

    print(f'✅ 数据已保存至 {output_path}')




if __name__ == "__main__":
    main()