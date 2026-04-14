# gene_common 说明

生成时间: 2026-03-26

`gene_common` 用来把 `gene_extract.sample_candidate_genes.tsv` 中的基因按群体聚合，输出两套交集结果：

- 至少 2 个群体都出现的基因
- 5 个群体都出现的基因

输出:

- `gene_common.support_matrix.tsv`
- `gene_common.at_least_2.tsv`
- `gene_common.at_least_2.bed`
- `gene_common.all_5.tsv`
- `gene_common.all_5.bed`
