# gene_share_peak 说明

生成时间: 2026-04-10

## 原理与作用
`gene_share_peak` 是本流程负责执行“信号提纯”和“共线性聚类”的主干 Python 脚本。
1. **信号平滑预处理：** 利用移动平均算法（Moving Average）在检测之前预先过滤掉数据的尖峰杂音，确保只关注具备平滑曲线宏观偏移的大范围隆起段。
2. **多层峰值提取：** 对输入的每个群体绘图表（`plot.tsv`），分别提取出在全基因组上的“信号局部极大值点”（Positive Peaks）和“极小值倒装点”（Negative Peaks）。抛弃绝对值一刀切，转而使用 `prominence` 显著突出度配合 `width` 底层全宽作为判定依据，彻底隔绝噪点。
3. **空间聚类（Cluster 合并）：** 将这 5 张表里剩余的精选大峰汇总在一条时间线上，针对同一染色体和方向，利用距离凝聚算法（Single Linkage Clustering），将相距近在咫尺（如 3MB 以内）的各个散点融合成一个组块（共有峰簇 / Shared Peak Cluster）。
4. **坐标外放与收割：** 统计出这块“共有区域”里涵盖了哪几个核心群体的样本，并通过求算该跨度极值边界向左右侧翼延展特定长度（如 1MB），最终在这个总宽域内从 GFF3 里面网罗出最终对应的基因库。

## 使用方法

通常此脚本不需要单独运行，它已配置在调度脚本当中。如果您希望独立调试，格式如下：

```bash
python gene_share_peak \
  --plot-tsv a.plot.tsv b.plot.tsv ... \
  --gff3 Canz.genome.gff3 \
  --output-dir ./output/gene_extract_SharePeak \
  --smooth-window 10 \
  --min-width 5 \
  --prominence 0.05 \
  --distance-window 3000000 \
  --flank 1000000 \
  --min-support 2
```

## 输入文件
- **plot.tsv**：由上游跑完的各 F2 群体关联指数文件，必需包含 `#CHROM` `POS` 及最重要的是包含算好的差列 `DELTA_SNP_INDEX`（第6列）。
- **GFF3**：标准的物种全基因组注释文件，此脚本仅对设定为 `gene` 的 feature 类型做相交检索。

## 输出文件与结果解读

> [!IMPORTANT]
> 重点注意：共有峰（Shared Peak）是指各样本原本的峰点（最高点）因为相邻满足距离组合在一起形成的新单位，也就是本输出中的 Cluster_XXXX 概念。

### 1. share_peak.all_clusters.tsv (全体峰簇统计全览表)
该文件记录了通过距离聚类组装出的**所有大大小小的簇**（包括单个样本产生的独立簇）。
- `cluster_id`：标识该共有簇组的唯一代号（例如：Cluster_0012_Chr01_positive）
- `sign`：指明该峰簇的大势方向（`positive`意味着上升、野生/突变型偏好等）
- `support_count`：该峰簇**涵盖的唯一样本群体数量**。这是您看**“有多少共有的峰？”**的依据！比如数值为 `2`，说明该山包趋势是 2 个样本共建的。
- `min_peak_pos`与`max_peak_pos`：该共有簇体构成点中最边缘顶点的位置收口。
- `peak_details`：展示簇中每个点来自于什么样本、在哪个准确位点以及它所携带的具体纵坐标振幅值。

### 2. share_peak.support_ge_X.bed (具有跨样本支撑的区域)
此为 BED 三列/五列常规结构区间集，它剔除了那种“孤零零只出现在某一个群体里面的特异独峰”，仅仅截取出您关注的 `min-support >= 2` 的群体交集区域段，您可用于喂给后端的 `gene_region_vcf` 和 Annovar 继续提取多态性。
- 第四列为所属源的簇名名称。
- 第五列为共有支持数量（您可用于做更硬核把关，提取 =5 的极端共同区间）。

### 3. share_peak.candidate_genes.support_ge_X.tsv (跨群体核心基因宝库)
对上述的高置信、有群体共撑片段进行了扫描，所捕获捞出的所有基因都会列在这张表格中。因为基因有时体型巨大或多个区间重叠，该表以**基因为单位**统计信息：
- `total_support_samples_count`：表示一共有几个多样本簇击中了并背书了这条基因。
- `total_support_samples`：将这些样本的具体标识直接罗列展示。
- `overall_direction_label`：用于标识这条基因所在的区域峰走势，如果这基因恰好在各样本的“只升”、“只降”中，就会被标示为 `positive_only` 或 `negative_only`，假如跨群样本出现反拉或倒错则为 `mixed`。
