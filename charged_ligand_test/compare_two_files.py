def find_common_lines(file_a, file_b, output_file):
    # 读取 a.txt 和 b.txt
    with open(file_a, 'r', encoding='utf-8') as fa:
        set_a = set(fa.read().splitlines())  # 读取所有行并转换为集合

    with open(file_b, 'r', encoding='utf-8') as fb:
        set_b = set(fb.read().splitlines())  # 读取所有行并转换为集合

    # 找出重叠部分
    common_lines = set_a.intersection(set_b)

    # 写入 c.txt
    with open(output_file, 'w', encoding='utf-8') as fc:
        for line in sorted(common_lines):  # 可选：排序后写入
            fc.write(line + '\n')

    print(f"✅ 共有 {len(common_lines)} 行，已保存至 {output_file}")

# 使用函数
find_common_lines('../data/subset2_2023/subset2_Ce.txt', 'ce4.txt', 'ce_combine.txt')
find_common_lines('../data/subset2_2023/subset2_Eu.txt', 'eu2.txt', 'eu_combine.txt')