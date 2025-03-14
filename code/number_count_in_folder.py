#count number of lines in file
#you need to specify the folder name

import os
import pandas as pd

def analyze_folder(folder_path):
    stats = []
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        
        # 仅处理文本文件
        if os.path.isfile(file_path) and filename.endswith('.txt'):
            with open(file_path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
                non_empty_lines = [line for line in lines if line.strip()]
                line_count = len(non_empty_lines)
                word_count = sum(len(line.split()) for line in non_empty_lines)
                stats.append({
                    'Filename': filename,
                    'Line Count': line_count,
                    'Word Count': word_count
                })
    return pd.DataFrame(stats)

# 要分析的两个文件夹路径
folder_path1 = './org_water_no_nitrate'
folder_path2 = './org_nitrate_no_water'
folder_path3 = './org_water_and_nitrate'
folder_path4 = './org_nitrate'
folder_path5 = './org_water'
folder_path6 = './water'
folder_path7 = './nitrate'

# 分析每个文件夹
df1 = analyze_folder(folder_path1)
df2 = analyze_folder(folder_path2)
df3 = analyze_folder(folder_path3)
df4 = analyze_folder(folder_path4)
df5 = analyze_folder(folder_path5)
df6 = analyze_folder(folder_path6)
df7 = analyze_folder(folder_path7)

# 保存到xlsx文件
output_file = 'analysis_results.xlsx'
with pd.ExcelWriter(output_file) as writer:
    df1.to_excel(writer, sheet_name='org_water_no_nitrate', index=False)
    df2.to_excel(writer, sheet_name='org_nitrate_no_water', index=False)
    df3.to_excel(writer, sheet_name='org_water_and_nitrate', index=False)
    df4.to_excel(writer, sheet_name='org_nitrate', index=False)
    df5.to_excel(writer, sheet_name='org_water', index=False)
    df6.to_excel(writer, sheet_name='water', index=False)
    df7.to_excel(writer, sheet_name='nitrate', index=False)

print(f'分析结果已保存到 {output_file}')
