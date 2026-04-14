# gene_extract 说明

生成时间: 2026-03-26

`gene_extract` 用来读取多个群体的 `plot.tsv`，找出满足 `|DELTA_SNP_INDEX| > 0.4` 的显著窗口，再把每个窗口中心上下游各延伸 `1 Mb`，最后从 `GFF3` 中提取候选基因。

输出:

- `gene_extract.significant_windows.tsv`
- `gene_extract.sample_regions.tsv`
- `gene_extract.window_genes.tsv`
- `gene_extract.sample_candidate_genes.tsv`
- `gene_extract.genes.bed`
- `gene_extract.genes.tsv`

使用:

```bash
python gene_extract \
  --plot-tsv a.plot.tsv b.plot.tsv c.plot.tsv d.plot.tsv e.plot.tsv \
  --gff3 /data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.gff3 \
  --output-dir ./output/gene_extract \
  --delta-threshold 0.4 \
  --flank 1000000
```
