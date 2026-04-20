# KASP Signal Design Workflow

本仓库用于从 BSA 的 `snpindex + VCF + reference fasta` 中筛选 KASP 标记位点。

## 目录结构
- `data/`: 配置、区间表、Snakemake rule 环境文件
- `doc/`: 使用说明与脚本文档
- `results/`: 正式输出
- `scripts/`: 核心执行脚本
- `temp/`: 运行日志与临时排查文件

## 主入口
```bash
conda run -n kasp-snakemake snakemake --use-conda -j 1
```

## 第一次使用
先看：
- `doc/20260420_kasp_snakemake_quickstart.md`

## GitHub 更新流程
本机修改后：
```bash
git add .
git commit -m "update kasp workflow"
git push
```

服务器更新：
```bash
git pull
```
