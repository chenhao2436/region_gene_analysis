#脚本用来提取snpindex的区间的，输入snpindex文件和染色体号及起止位置就提取出来了
def extract_variants(file_path, chromosome, start_pos, end_pos):
    variants = []
    with open(file_path, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            file_chromosome = data[0]
            position = int(data[1])
            if file_chromosome == chromosome and start_pos <= position <= end_pos:
                variants.append(line)
    return variants

def save_variants_to_file(variants, output_file):
    with open(output_file, 'w') as file:
        for variant in variants:
            file.write(variant)

file_path = '/data/wangxx/wangxx-data/2023qF9/2023qiuF9/00.mergeRawFq/F9LL/total.vcf.format'
chromosome = 'Chr01'
start_pos = 78500000
end_pos = 82000000
output_file = 'output.txt'

variants = extract_variants(file_path, chromosome, start_pos, end_pos)
save_variants_to_file(variants, output_file)
