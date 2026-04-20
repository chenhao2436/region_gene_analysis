# 步骤4: 扩展区域（前后各100bp）
bedtools slop -i Chr07.bed -g 62.genome.Chr.fa.fai -b 100 > isolated_snps_200bp.bed

# 步骤5: 提取序列
bedtools getfasta -fi 62.genome.Chr.fa -bed isolated_snps_200bp.bed -fo isolated_snps_flanking.fasta

# 步骤6: 在序列中标记SNP位置
awk '
BEGIN {
  # 读取SNP信息
  while (getline < "isolated_snps.txt") {
    chrom = $1
    pos = $2
    ref = $3
    alt = $4
    snp_info[chrom":"pos] = ref ">" alt
  }
}
/^>/ {
  # 提取染色体和位置信息
  match($0, /:([0-9]+)-([0-9]+)/, coords)
  start = coords[1]
  end = coords[2]
  snp_pos = start + 100  # SNP在序列中的位置
  
  # 提取染色体名称
  match($0, /^>(.*):/, chr_match)
  chrom_name = chr_match[1]
  
  # 获取SNP信息
  snp_key = chrom_name ":" snp_pos
  snp_change = snp_info[snp_key]
  
  # 打印新的头信息
  print ">" chrom_name ":" snp_pos " (SNP@101) " snp_change
  next
}
{
  # 处理序列，在SNP位置添加方括号
  seq = $0
  # 确保序列长度足够
  if (length(seq) >= 101) {
    left = substr(seq, 1, 100)
    snp_base = substr(seq, 101, 1)
    right = substr(seq, 102)
    print left "[" snp_base "]" right
  } else {
    # 处理边界情况
    print seq
  }
}
' isolated_snps_flanking.fasta > final_kasp_sequences.fasta
