#This script is to find the inorganic anions such as nitrate, sulfate, perchlorate, chloride in the Ln neutral ligand complexes.
#The script will read the CSD database and find the inorganic anions in the ligand complexes.

#this script will read subset2.
#file path = ../data/subset2_2023/subset2_Ln.txt

from ccdc.io import MoleculeReader
import re
from time import process_time
import os
from ccdc.molecule import Molecule
from ccdc import io
import pandas as pd
from collections import defaultdict
import math
import numpy as np

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

def is_perchlorate(atomCl):
    flag = False
    if all(atom.atomic_symbol == 'O' for atom in atomCl.neighbours) and len(atomCl.neighbours) == 4:
        flag = True
    return flag

def is_sulfate(atomS):
    flag = False
    if all(atom.atomic_symbol == 'O' for atom in atomS.neighbours) and len(atomS.neighbours) == 4:
        flag = True
    return flag

def is_chloride(atomCl):
    flag = False
    if len(atomCl.neighbours) == 0:
        flag = True
    return flag

def is_bromide(atomBr):
    flag = False
    if len(atomBr.neighbours) == 0:
        flag = True
    return flag

def is_iodide(atomI):
    flag = False
    if len(atomI.neighbours) == 0:
        flag = True
    return flag

def identify_num_of_ligand_and_anions_in_fs(mol, ELEMENT):
    # 打开CSD数据库
    mol_reader = io.MoleculeReader('csd')

    contains_organic_ligand = False
    organic_ligand_charges = []  # 存储所有有机配体的电荷
    num_organic_ligand = 0  # 记录有机配体的数量
    num_nitrate_in_fs = 0  # 记录 fs 中的硝酸根离子数量
    num_nitrate_in_ss = 0  # 记录 ss 中的硝酸根离子数量
    num_sulfate_in_fs = 0
    num_chloride_in_fs = 0
    num_perchlorate_in_fs = 0
    num_bromide_in_fs = 0
    num_iodide_in_fs = 0
    coordination_number = 0

    for com in mol.components:
        for atom in com.atoms:
            if atom.atomic_symbol == 'N' and is_nitrate(atom):
                num_nitrate_in_ss += 1
            if atom.atomic_symbol == ELEMENT:
                coordination_number = len(atom.neighbours)

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
            num_organic_ligand += 1

        for atom in com.atoms:
            if atom.atomic_symbol == 'N' and is_nitrate(atom):
                num_nitrate_in_fs += 1
            elif atom.atomic_symbol == 'S' and is_sulfate(atom):
                num_sulfate_in_fs += 1
            elif atom.atomic_symbol == 'Cl':
                if is_chloride(atom):
                    num_chloride_in_fs += 1 
                if is_perchlorate(atom):
                    num_perchlorate_in_fs += 1
            elif atom.atomic_symbol == 'Br' and is_bromide(atom):
                num_bromide_in_fs += 1
            elif atom.atomic_symbol == 'I' and is_iodide(atom):
                num_iodide_in_fs += 1

    return contains_organic_ligand, num_organic_ligand, num_nitrate_in_fs, num_nitrate_in_ss, num_sulfate_in_fs, num_perchlorate_in_fs, num_chloride_in_fs, num_bromide_in_fs, num_iodide_in_fs, coordination_number

def classify_complex(num_nitrate_fs, num_nitrate_ss, coordination_number, contains_organic_ligand):
    """
    依据硝酸根的分布情况和配位数，分类复合物：
    - 3:1 -> 仅在 SS 存在硝酸根
    - 2:1 -> FS: 1-2 个硝酸根, SS: 1-2 个硝酸根
    - 1:1 -> FS: 2-3 个硝酸根, SS: 0-1 个
    - other_inorganic -> 不含有机配体
    - other_organic -> 含有有机配体但不符合上述规则
    """
    if not contains_organic_ligand:
        return "other_inorganic", coordination_number  # 归入无机复合物

    if num_nitrate_fs == 0 and num_nitrate_ss > 0:
        return "3:1", coordination_number
    elif 1 <= num_nitrate_fs <= 2 and 1 <= num_nitrate_ss <= 2:
        return "2:1", coordination_number
    elif 2 <= num_nitrate_fs <= 3 and num_nitrate_ss <= 1:
        return "1:1", coordination_number
    else:
        return "other_organic", coordination_number  # 归入其他有机复合物

def screen_com(ELEMENT, cn_summary):
    """ 处理单个金属元素的配体统计，并分类 """
    filepath = f'../data/subset2_2023/nitrate_{ELEMENT}.txt'
    mol_reader = MoleculeReader(filepath, format='identifiers')

    results = []
    cn_values = {
        "3:1": [], "2:1": [], "1:1": [],
        "other_organic": [], "other_inorganic": []
    }  # 存储各类别 Coordination Number

    for mol in mol_reader:
        try:
            is_org, num_organic_ligand, num_nitrate_fs, num_nitrate_ss, num_sulfate, num_perchloriate, num_chloride, num_bromide, num_iondide,cn = identify_num_of_ligand_and_anions_in_fs(mol, ELEMENT)
            complex_type, coordination_number = classify_complex(num_nitrate_fs, num_nitrate_ss, cn, is_org)
            
            results.append({
                "CSD_Code": mol.identifier,
                "Contains_Organic_Ligand": is_org,
                "Num_Organic_Ligands": num_organic_ligand,
                "Num_Nitrate_FS": num_nitrate_fs,
                "Num_Nitrate_SS": num_nitrate_ss,
                "Num_Sulfate_FS": num_sulfate,
                "Num_Perchlorate_FS": num_perchloriate,
                "Num_Chloride_FS": num_chloride,
                "Num_Bromide_FS": num_bromide,
                "Num_Iodide_FS": num_iondide,
                "Coordination_Number": coordination_number,
                "Complex_Type": complex_type
            })

            # 记录 Coordination Number
            cn_values[complex_type].append(coordination_number)

        except RuntimeError:
            print(f"⚠️ 跳过无效的 CSD 代码: {mol.identifier}")
            continue

    # 计算每个类别的平均 Coordination Number
    avg_cn = {key: np.mean(values) if values else None for key, values in cn_values.items()}

    # 打印 Coordination Number 统计信息
    print(f"\n📊 {ELEMENT} Coordination Number 平均值:")
    for key, value in avg_cn.items():
        print(f"  - {key}: {value:.2f}" if value is not None else f"  - {key}: 无数据")

    # 存储 CN 平均值到全局字典
    cn_summary[ELEMENT] = avg_cn

    # 转换为 DataFrame 并保存 CSV 详细数据
    df = pd.DataFrame(results)
    output_path = f'../data/anions_analysis/{ELEMENT}.csv'
    df.to_csv(output_path, index=False)
    print(f"✅ {ELEMENT} 结果已保存至 {output_path}")

def main():
    start = process_time()
    
    LIST_OF_ELEMENT = [
        'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu'
    ]
    
    folder_path = '../data/nitrate_org_ligand_nitrate/'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    cn_summary = {}  # 存储所有元素的 CN 统计结果

    for element in LIST_OF_ELEMENT:
        print(f"🔍 处理 {element} ...")
        screen_com(element, cn_summary)

    # 转换 CN 统计结果为 DataFrame
    cn_df = pd.DataFrame.from_dict(cn_summary, orient='index')

    # 保存到 Excel
    excel_path = f'../data/nitrate_org_ligand_nitrate/coordination_numbers.xlsx'
    cn_df.to_excel(excel_path)
    print(f"📊 所有 Coordination Number 统计结果已保存至 {excel_path}")

    end = process_time()
    print('⌛ 运行时间: %s 秒' % (end - start))


if __name__ == "__main__":
    main()