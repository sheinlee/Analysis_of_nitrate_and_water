LIST_OF_ELEMENT = [
    'La', 'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
            ]  

def read_file(filepath):
    with open(filepath, 'r') as file:
        return file.readlines()

def write_file(filepath, lines):
    with open(filepath, 'w') as file:
        file.writelines(lines)

def merge_files(output_filepath):
    all_lines = set()
    for ELEMENT in LIST_OF_ELEMENT:
        filepath = f'data/subset3/nitrate_hits_{ELEMENT}.txt'
        lines = read_file(filepath)
        all_lines.update(lines)
    
    write_file(output_filepath, all_lines)

output_filepath = 'data/subset3/merged_nitrate_hits.txt'  # 合并后的输出文件路径

merge_files(output_filepath)
