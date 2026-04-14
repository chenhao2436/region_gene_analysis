# gene_region_vcf 说明

生成时间: 2026-03-27

## 作用

`gene_region_vcf` 负责把交集基因 BED 对应的区域从 joint VCF 中提取出来，并用外部 GTF 构建 Annovar 注释数据库，输出变异注释结果。

## 输入

- 交集基因 BED
- 交集基因 TSV
- joint VCF / VCF.GZ
- 参考基因组 FASTA
- 外部 GTF

## 输出

- `gene_region_vcf.<label>.vcf.gz`
- `gene_region_vcf.<label>.vcf.gz.tbi`
- `gene_region_vcf.<label>.vcf`
- `gene_region_vcf.<label>.avinput`
- `gene_region_vcf.<label>.variant_function`
- `gene_region_vcf.<label>.exonic_variant_function`
- `gene_region_vcf.<label>.genes.tsv`
- `gene_region_vcf.<label>.command.log`

## 默认 GTF

- `/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract/output/gene_region_vcf/at_least_2/Canz.genome.gtf`

## 说明

- 脚本不再内部做 GFF3->GTF 转换。
- Annovar 仍然按 `gtfToGenePred -> retrieve_seq_from_fasta.pl -> convert2annovar.pl -format vcf4 -> annotate_variation.pl` 的顺序执行。