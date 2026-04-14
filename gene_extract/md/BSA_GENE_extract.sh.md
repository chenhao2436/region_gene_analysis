# BSA_GENE_extract.sh 说明

生成时间: 2026-03-27

这是总控 shell。你只需要维护脚本顶部的变量，然后运行一次，就会依次调用：

1. `gene_extract`
2. `gene_common`
3. `gene_region_vcf`
4. `gene_out`

## 当前关键输入

- 5 个群体的 `plot.tsv`
- `Canz.genome.gtf` 作为外部 GTF
- `Canz.genome.fa`
- `hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz`
- `AT_Canz_reciprocal_best_hits.txt`
- `TAIR_gene.annotate.unique_v1.2`

## 关键参数

- `DELTA_THRESHOLD=0.4`
- `FLANK_BP=1000000`
- `MIN_SUPPORT=2`
- `ALL_SUPPORT=5`