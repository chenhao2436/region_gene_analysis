# BSA Gene Extract Workflow

生成时间: 2026-03-26

这套脚本用于把 5 个群体的 BSA `plot.tsv` 结果继续往下游推进，目标是：

1. 按 `|DELTA_SNP_INDEX| > 0.4` 提取显著信号窗口。
2. 以每个显著窗口中心上下游各延伸 `1 Mb`。
3. 从 `Canz.genome.gff3` 中提取这些区间内的基因。
4. 对多个群体的候选基因求交集。
5. 提取交集基因本身范围内的所有变异。
6. 用 Annovar 做变异注释。
7. 再把辣椒基因和拟南芥功能注释整合成便于阅读的表格。

建议执行:

```bash
bash /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract/scripts/BSA_GENE_extract.sh
```
