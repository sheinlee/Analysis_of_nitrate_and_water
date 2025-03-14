from ccdc import io
from ccdc.molecule import Molecule

def contains_carbon_carbon_bond(molecule):
    """
    判断分子中是否包含碳-碳键。
    """
    for bond in molecule.bonds:
        atoms = bond.atoms
        if atoms[0].atomic_symbol == 'C' and atoms[1].atomic_symbol == 'C':
            return True
    return False

# 打开CSD数据库
mol_reader = io.MoleculeReader('csd')

# 读取分子
mol = mol_reader.molecule('ACUHAI')
ELEMENT = 'La'  # 目标金属元素
    
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
            # connected_atom = bond.atoms(atom)
            # atoms_to_remove.add(connected_atom)
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
    print(len(com.atoms))
    print(com.atoms)
    
    com.assign_bond_types()
    with io.MoleculeWriter('com.mol2') as writer:
        writer.write(com)
    com.set_formal_charges()
    print(com.formal_charge)
    print(contains_carbon_carbon_bond(com))
    print('==============')
    # break

# # 将新分子保存为 .mol2 文件
# output_filename = 'ligand.mol2'
# with io.MoleculeWriter(output_filename) as writer:
#     writer.write(new_mol)

# print(f'配体已保存为 {output_filename}')


from ccdc.io import EntryReader
csd_reader = EntryReader('CSD')
mol = csd_reader.entry('ADOROE')
print(mol.chemical_name)