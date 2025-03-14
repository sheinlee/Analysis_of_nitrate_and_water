import difflib

def read_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.readlines()

def compare_files(file1_path, file2_path, add_output, red_output):
    file1_lines = read_file(file1_path)
    file2_lines = read_file(file2_path)
    
    d = difflib.Differ()
    diff = list(d.compare(file1_lines, file2_lines))
    
    with open(add_output, 'w', encoding='utf-8') as add_file, open(red_output, 'w', encoding='utf-8') as red_file:
        for line in diff:
            if line.startswith('+ '):
                add_file.write(line[2:])  # 去掉前面的'+ '
            elif line.startswith('- '):
                red_file.write(line[2:])  # 去掉前面的'- '

# 使用方法
file1_path = 'subset2_cn2_La.txt'
file2_path = 'subset2_La.txt'
add_output = 'add_diff.txt'
red_output = 'red_diff.txt'

compare_files(file1_path, file2_path, add_output, red_output)
