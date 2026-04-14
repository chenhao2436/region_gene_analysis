#我的目标：给出染色体，找到信号最强的位点，对位点进行1MB的扩宽，然后找到这1MB内的基因，
以及vcf文件和注释文件，然后还有功能注释文件
这需要我提前把vcf文件注释上，还需要和拟南芥的功能进行比对
#!/bin/bash
# 定义输入文件和参数
Chromosome="chr1"  # 替换为你感兴趣的染色体
vcf_file="input.vcf"  # 替换为你的VCF文件路径
annotation_file="annotation.txt"  # 替换为你的注释文件路径
output_file="output.txt"  # 输出文件路径
REGION_SIZE=1000000  # 1MB
#1. 找到信号最强的位点
strongest_site=$(awk -v chrom="$Chromosome" '$1 == chrom {print $0}' "$vcf_file" | sort -k5,5nr | head -n 1)
if [ -z "$strongest_site" ]; then
    echo "未找到染色体 $Chromosome 上的位点。"
    exit 1
fi
#2. 扩宽位点的1MB区间，并找到该区间内的基因
position=$(echo "$strongest_site" | awk '{print $2}')
start=$((position - REGION_SIZE))
end=$((position + REGION_SIZE))
if [ $start -lt 0 ]; then
    start=0
fi
#3. 提取区间内的gff文件
awk -v chrom="$Chromosome" -v start="$start" -v end="$end" '$1 == chrom && $4 >= start && $5 <= end {print $0}' "$annotation_file" > genes_in_region.txt
#4. 提取区间内的vcf文件
awk -v chrom="$Chromosome" -v start="$start" -v end="$end" '$1 == chrom && $2 >= start && $2 <= end {print $0}' "$vcf_file" > variants_in_region.vcf
