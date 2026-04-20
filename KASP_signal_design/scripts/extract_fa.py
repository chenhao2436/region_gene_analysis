import pyfaidx
import argparse

def extract_sequences(input_file, ref_genome_path, output_file):
    # 加载参考基因组
    genome = pyfaidx.Fasta(ref_genome_path)
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # 写入表头
        f_out.write("ID\tRef\tAlt\tSequence\n")
        
        for line in f_in:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
                
            chrom, pos, ref, alt = parts
            pos = int(pos)
            
            # 尝试不同的染色体命名格式
            chrom_formats = [
                chrom.lower().replace("chr", ""),  # "10"
                chrom.replace("Chr", "chr"),       # "chr10"
                chrom                             # 原始格式
            ]
            
            seq_region = None
            for chrom_fmt in chrom_formats:
                try:
                    # 获取前后100bp范围
                    start = max(1, pos - 100)
                    end = pos + 100
                    
                    # 提取序列
                    seq_region = genome[chrom_fmt][start-1:end].seq
                    break
                except KeyError:
                    continue
            
            if seq_region is None:
                print(f"错误: 无法找到染色体 {chrom} 的序列")
                continue
                
            # 验证参考碱基
            if len(seq_region) > 100:
                ref_in_genome = seq_region[100]
                if ref_in_genome != ref:
                    print(f"警告: {chrom}:{pos} 参考基因组碱基为 {ref_in_genome} (非 {ref})")
            
            # 构建变异标记序列
            left_seq = seq_region[:100]
            right_seq = seq_region[101:201] if len(seq_region) > 101 else seq_region[101:]
            marked_seq = f"{left_seq}[{ref}/{alt}]{right_seq}"
            
            # 生成ID
            variant_id = f"{chrom}:{pos}"
            
            # 写入结果
            f_out.write(f"{variant_id}\t{ref}\t{alt}\t{marked_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='提取SNP位点前后100bp序列')
    parser.add_argument('--input', required=True, help='输入SNP文件路径')
    parser.add_argument('--output', required=True, help='输出文件路径')
    parser.add_argument('--ref', required=True, help='参考基因组FASTA文件路径')
    
    args = parser.parse_args()
    
    extract_sequences(args.input, args.ref, args.output)
    print(f"处理完成！结果已保存至 {args.output}")
