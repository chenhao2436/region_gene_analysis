# Region Gene Analysis

[English](README_EN.md) | **简体中文**

一个面向植物 BSA / 群体遗传定位研究的分析脚本集合，覆盖 SNP-index / ED 计算与绘图、BSA 信号区域基因提取、跨群体共有峰筛选、区域 VCF 提取与变异注释、候选基因功能整合，以及若干 GVCF 联合分型和数据审计工具。

> **项目状态：研究工具箱 / Alpha。** 该仓库由多个研究阶段逐步积累而来，并不是一个统一配置的一键流水线。多个历史 Bash 调度脚本仍包含原服务器的绝对路径；推荐优先直接调用提供通用命令行参数的 Python 子模块。Annovar、参考基因组、GFF3/GTF、VCF 等外部资源需要用户自行准备。

## 主要模块

| 模块 | 主要用途 | 主要输入 | 主要输出 |
|---|---|---|---|
| BSA 计算与绘图 | F2 SNP-index / ED 计算及静态图 | VCF、亲本和混池样本列 | SNP-index / ED 表、PDF |
| `gene_extract` | 根据显著 BSA 窗口扩展区间并提取基因 | 多群体 plot.tsv、GFF3 | 窗口、区间、候选基因表 |
| `gene_common` | 按群体支持数汇总交集基因 | 候选基因表 | 至少 N 个群体支持、全部群体支持的 BED/TSV |
| `gene_region_vcf` | 提取候选基因区域变异并运行 Annovar 注释 | BED、joint VCF、FASTA、GTF | 区域 VCF、avinput、功能注释结果 |
| `gene_out` | 整合辣椒基因与拟南芥 RBH/TAIR 注释和变异注释 | 基因交集、RBH、TAIR、Annovar 结果 | 候选基因综合表 |
| `gene_extract_SharePeak` | 平滑 BSA 信号、检测峰并跨样本聚类 | 多群体 plot.tsv、GFF3 | Shared Peak 簇、支持区间、候选基因 |
| `interact_html` | 生成交互式 BSA HTML 图 | 多群体 plot.tsv | 交互式 HTML |
| GVCF/VCF 工具 | GVCF 列表整理、联合分型和文件检查 | GVCF、参考基因组 | joint VCF、审计报告 |

各模块有独立输入和参数，不要求必须全部运行。

## 仓库结构概览

```text
region_gene_analysis/
├── F2_BSA.py
├── bsa.fig.py
├── auto_bsa_*.sh
├── audit_bsa_sample_chain.sh
├── gene_extract/
│   ├── scripts/              # gene_extract / gene_common / gene_region_vcf / gene_out
│   └── md/                   # 模块说明
├── gene_extract_SharePeak/
│   ├── scripts/              # Shared Peak 检测和聚类
│   └── md/
├── interact_html/
│   ├── scripts/              # 交互式 BSA 图
│   └── md/
├── ploidy30_bsa_bundle_20260318_upload_v2/
│   └── gvcf_plidy_30_BSA.py
└── gvcf2vcf_genome_gpt_3m6d.sh
```

仓库还保留了部分实验记录、备份、历史脚本和压缩包。请以本文列出的当前模块入口为准，不要将历史副本当作独立产品入口。

## 软件环境

Python 环境文件覆盖当前主要 Python 分析模块的常用依赖：

```bash
conda env create -f environment.yml
conda activate region-gene-analysis
```

核心 Python 依赖：

- Python 3.10+
- NumPy
- SciPy（Shared Peak 峰值检测）
- Matplotlib（静态绘图）
- Plotly（交互式绘图）

部分模块还需要外部命令行软件或数据库：

| 工具/资源 | 用途 |
|---|---|
| bcftools | VCF 区间提取、样本和基因型查询 |
| tabix / bgzip | VCF 压缩索引与区域访问 |
| GATK 4 | GVCF 联合分型辅助脚本 |
| Annovar | 变异功能注释；需用户自行安装并遵守其许可 |
| UCSC `gtfToGenePred` | GTF 转换 |
| 参考基因组 FASTA | 序列、坐标和注释一致性基础 |
| GFF3 / GTF | 基因坐标与注释 |
| TAIR 注释表及 RBH 文件 | 辣椒候选基因与拟南芥功能映射 |

上述外部工具未全部由 Conda 环境文件安装。运行对应模块前，请分别确认工具在 PATH 中可用，并核对实际版本。

## 推荐使用流程：从 BSA 信号到候选基因

### 1. 准备 BSA plot.tsv

由上游 BSA 流程生成各群体的 `plot.tsv`。同一次比较应尽量使用相同参考基因组、染色体命名和一致的窗口参数。

### 2. 根据阈值提取候选区域基因

`gene_extract/scripts/gene_extract` 的主要参数包括：

- `--plot-tsv`：一个或多个群体的 plot.tsv。
- `--gff3`：参考基因组 GFF3。
- `--output-dir`：输出目录。
- `--delta-threshold`：DELTA_SNP_INDEX 阈值（当前历史说明示例为 0.4）。
- `--flank`：信号窗口两侧扩展长度（示例为 1,000,000 bp）。

示例：

```bash
python gene_extract/scripts/gene_extract \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv sample_168.plot.tsv \
  --gff3 reference.gff3 \
  --output-dir ./results/gene_extract \
  --delta-threshold 0.4 \
  --flank 1000000
```

主要输出包括显著窗口、样本区域、窗口基因、样本候选基因和 BED/TSV 文件。具体字段见 `gene_extract/README.md` 和 `gene_extract/md/`。

### 3. 汇总跨群体支持基因

```bash
python gene_extract/scripts/gene_common \
  --sample-candidate-tsv ./results/gene_extract/gene_extract.sample_candidate_genes.tsv \
  --output-dir ./results/gene_common \
  --min-support 2 \
  --all-support 5
```

参数含义：`--min-support` 定义进入共有候选集的最少群体数；`--all-support` 定义全部群体支持的数量。请根据实际输入群体数量设置，不要固定假设每次都是 5 个群体。

### 4. Shared Peak 备选策略

`gene_extract_SharePeak/scripts/gene_share_peak` 使用平滑和峰值检测识别正、负信号峰，再依据染色体、方向和距离进行聚类。

```bash
python gene_extract_SharePeak/scripts/gene_share_peak \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv sample_168.plot.tsv \
  --gff3 reference.gff3 \
  --output-dir ./results/shared_peak \
  --smooth-window 10 \
  --min-width 5 \
  --prominence 0.05 \
  --distance-window 3000000 \
  --flank 1000000 \
  --min-support 2
```

该模块与固定 DELTA 阈值的 `gene_extract` 是两种不同的信号提取策略，结果不应在没有方法说明的情况下混为一谈。

### 5. 交互式 BSA 图

```bash
python interact_html/scripts/interact_html \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv \
  --output-dir ./results/interactive
```

可选参数 `--sample-labels` 用于指定显示标签；`--cluster-tsv` 和 `--min-cluster-support` 用于在支持该功能的当前脚本版本中叠加 Shared Peak 簇信息。输出通常包括 combined 和 subplots HTML。请使用 `--help` 核对具体脚本版本的参数。

## 区域 VCF 提取与 Annovar 注释

`gene_extract/scripts/gene_region_vcf` 使用 bcftools 按候选基因 BED 提取 joint VCF 区域，并调用外部 Annovar 工具链进行注释。所需输入：

- 候选基因 BED 和对应 TSV。
- bgzip 压缩并索引的 joint VCF。
- 与 VCF 匹配的参考 FASTA。
- 与基因组版本匹配的 GTF。
- Annovar 安装目录及可执行脚本。
- `bcftools`、`tabix`、`gtfToGenePred`。

示例：

```bash
python gene_extract/scripts/gene_region_vcf \
  --gene-bed ./results/gene_common/gene_common.at_least_2.bed \
  --gene-table ./results/gene_common/gene_common.at_least_2.tsv \
  --vcf joint.vcf.gz \
  --reference-fasta reference.fa \
  --gtf reference.gtf \
  --annovar-dir /path/to/annovar \
  --bcftools-bin bcftools \
  --tabix-bin tabix \
  --gtfToGenePred-bin gtfToGenePred \
  --label at_least_2 \
  --output-dir ./results/variant_annotation
```

脚本当前源码保留了原服务器的 GTF、Annovar 和 bcftools 默认路径。建议始终显式提供 `--gtf`、`--annovar-dir` 和 `--bcftools-bin`，不要依赖这些默认值。

## 重要限制

- 本仓库是多模块研究工具箱，不是单一一键工作流；每个模块应独立运行并检查输出。
- 根目录及子目录中的部分 Bash 总控脚本仍写有 `/data2/chenh/` 等绝对路径。它们在其他机器上可能无法直接运行。
- `gene_region_vcf` 依赖 Annovar 等外部工具；Annovar 不随仓库分发。
- TAIR 注释、辣椒-拟南芥 RBH 对照、参考基因组和注释文件需用户准备，且要确认版本匹配。
- BSA 阈值、窗口大小、峰宽、prominence、距离聚类和支持群体数均会影响结果。它们是分析参数，不是通用生物学结论。
- 仓库中存在历史脚本、备份文件、实验记录和旧版配置，可能与当前主入口不同。
- 当前交付没有对每个子模块逐一进行真实数据端到端测试；README 描述的是静态源码审查确认的功能和已记录用法。

## 使用前检查清单

1. 参考基因组 FASTA、VCF、GFF3/GTF 必须来自同一组装版本。
2. 样本名必须与 VCF header 一致。
3. VCF 应具备正确索引；区域访问前检查 `.tbi` / `.csi`。
4. 检查 plot.tsv 的列顺序，尤其是 DELTA_SNP_INDEX 和 SNP 数量列。
5. 在小区域或少量样本上先跑通，再进行全基因组或大样本分析。
6. 保留完整命令、参数、工具版本和日志，以便复现。

## 路线图

- [ ] 将硬编码 Bash 包装脚本改为通用参数或配置文件。
- [ ] 为各模块提供统一的示例配置和输入文件模板。
- [ ] 为核心 Python 模块增加合成测试数据和自动化测试。
- [ ] 清理历史缓存与 Python 字节码文件的版本管理。
- [ ] 为各模块记录软件版本、默认参数和方法学边界。

## 问题反馈

提交问题时，请标明具体模块、运行命令、软件版本、输入文件列定义和日志。请勿将大型测序数据、未授权数据或私密服务器配置提交到公开仓库。
