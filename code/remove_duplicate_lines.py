def remove_duplicate_lines(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    unique_lines = set(lines)

    with open(filepath, 'w') as file:
        file.writelines(unique_lines)

# 示例使用：
# filepath = 'data/subset2_2023/subset2_example_element.txt'  # 替换为实际的文件路径
# remove_duplicate_lines(filepath)
LIST_OF_ELEMENT = [
    'La', 'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
            ]  
for ELEMENT in LIST_OF_ELEMENT:
    filepath = f'data/subset3/nitrate_hits_{ELEMENT}.txt'
    remove_duplicate_lines(filepath)

