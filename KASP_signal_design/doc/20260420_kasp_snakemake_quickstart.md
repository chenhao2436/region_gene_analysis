# 20260420 KASP Snakemake Quickstart

## 1. 目标
这套流程把 `scripts/` 里的 4 个脚本交给 Snakemake 调度，统一管理配置、输出和日志。

## 2. 为什么这样配环境
- `base` 保持干净
- 新建 `kasp-snakemake` 只装 `snakemake`
- 真正跑 `bcftools/tabix/pyfaidx` 的是 Snakemake 自动创建的 rule 环境

## 3. 创建驱动环境
```bash
conda create -n kasp-snakemake -c conda-forge -c bioconda snakemake -y
```

检查版本：
```bash
conda run -n kasp-snakemake snakemake --version
```

## 4. 配置文件位置
- 主配置：`data/config.yaml`
- 区间表：`data/regions.tsv`
- rule 环境：`data/kasp_rule_env.yaml`

## 5. dry-run
先只展开流程，不实际运行：
```bash
conda run -n kasp-snakemake snakemake -n
```

## 6. 正式运行
```bash
conda run -n kasp-snakemake snakemake --use-conda -j 1
```

只跑最终目标也可以：
```bash
conda run -n kasp-snakemake snakemake --use-conda -j 1 results/chr05_peak/04_kasp_sequences.tsv
```

## 7. 输出在哪里
- 正式结果：`results/chr05_peak/`
- 运行日志：`temp/chr05_peak/`

## 8. 主要输出
- `01_ranked_candidates.tsv`
- `01_rank_summary.tsv`
- `02_checked_sites.tsv`
- `02_pass_sites.tsv`
- `03_extract_fa_input.tsv`
- `03_recommended_markers.tsv`
- `04_kasp_sequences.tsv`

## 9. 常见报错
- `MissingInputException`: 检查 `data/config.yaml` 路径是否真实存在
- `bcftools query 失败`: 检查 VCF 索引和样本名
- `03_extract_fa_input.tsv` 为空: 说明没有位点通过过滤
- `CondaEnvException`: 检查服务器 conda 是否可用、镜像源是否通

## 10. GitHub + 服务器更新
本机推送：
```bash
git add .
git commit -m "update kasp snakemake workflow"
git push
```

服务器更新：
```bash
git pull
```
