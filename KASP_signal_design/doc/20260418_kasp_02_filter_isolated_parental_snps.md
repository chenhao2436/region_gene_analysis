# 20260418_kasp_02_filter_isolated_parental_snps.py

## 运行场景
用于对候选位点做 KASP 落地前的硬过滤, 保留双亲纯合差异、双等位 SNP、且前后 50 bp 没有其他变异的孤立位点。

## 脚本目标
- 逐侧按排序依次检查候选位点
- 先检查前 20 个, 不足时自动补位
- 过滤 INDEL 和多等位位点
- 保留双亲纯合且不同的位点
- 排除前后 50 bp 存在其他变异的位点

## 输入文件
- `--candidates`: `01_ranked_candidates.tsv`
- `--vcf`: bgzip + tabix 索引的 VCF 文件

## 输出文件
- `--output-checked`: `02_checked_sites.tsv`
- `--output-pass`: `02_pass_sites.tsv`

## 参数解释
- `--parent1`: 双亲 1 样本名
- `--parent2`: 双亲 2 样本名
- `--high-pool`: 高池样本名
- `--low-pool`: 低池样本名
- `--flank-bp`: 邻域不能有其他变异的范围, 默认 50
- `--target-markers-per-side`: 每侧目标位点数, 默认 2
- `--initial-rank-limit`: 优先尝试的前筛数量, 默认 20

## 使用示例
```bash
python 20260418_kasp_02_filter_isolated_parental_snps.py \
  --candidates 01_ranked_candidates.tsv \
  --vcf sample.vcf.gz \
  --parent1 C3_Parent \
  --parent2 C6_Parent \
  --high-pool Hunchi_H \
  --low-pool Hunchi_L \
  --output-checked 02_checked_sites.tsv \
  --output-pass 02_pass_sites.tsv
```

## 结果表字段说明
`02_checked_sites.tsv` 和 `02_pass_sites.tsv` 共用字段:
- `side`: `left` 或 `right`
- `chrom`: 染色体
- `pos`: 位点坐标
- `ref/alt`: 目标位点等位基因
- `delta_snp_index`: 候选排序依据
- `parent1_gt/parent2_gt`: 两个亲本基因型
- `high_pool_gt/low_pool_gt`: 高低池基因型
- `filter_status`: `PASS` 或 `FAIL`
- `filter_reason`: 过滤原因
- `checked_rank`: 当前侧排序名次

## 常见报错与排查
- `bcftools query 失败`: 检查 VCF、索引和样本名是否正确
- `target_site_not_found`: 候选位点在 VCF 中查不到
- `multiallelic_site`: 该位点 ALT 含多个等位基因
- `indel_or_non_snp_site`: 该位点不是单碱基 SNP
- `parent1_not_homozygous` / `parent2_not_homozygous`: 亲本不纯合
- `parents_share_same_genotype`: 两个亲本基因型相同
- `nearby_variant_present`: 目标位点前后 50 bp 内存在其他变异
