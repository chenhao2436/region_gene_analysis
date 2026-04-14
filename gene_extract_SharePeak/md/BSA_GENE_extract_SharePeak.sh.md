# BSA_GENE_extract_SharePeak.sh 说明

生成时间: 2026-04-10

## 原理与作用

这是一个用于基于**趋势（Trend）**特征提取跨样本共有峰（Shared Peaks）的总控 Shell 脚本。
它负责统一调度后端的 Python 核心处理逻辑 `gene_share_peak`，并将关键配置参数（聚类距离、隆起显著度、延展长度）集中在文件顶部以便于整体控制。本流程摒弃了原有单纯依赖 `|DELTA_SNP_INDEX| > 0.4` 的一刀切硬阈值方案，转而使用局部信号凸起提取加空间聚类的算法，使得只要呈现共有起伏趋势（即便绝对值波动小）也能被鉴定和保留。

## 使用方法

由于本脚本内置了所有的绝对路径参数，您只需要直接运行即可：

```bash
bash /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/scripts/BSA_GENE_extract_SharePeak.sh
```

**环境依赖**：运行该脚本必须在含有 `numpy` 和 `scipy` 包的 Python 3.x 环境中运行，因为后台的峰值检测依赖了 `scipy.signal.find_peaks` 函数。

## 关键参数解读

您可以编辑该脚本修改以下的全局变量配置：
- `PROMINENCE=0.05`: 峰的“突出度（相对于其局部谷底的高度）”。默认从 0.05 开始，代表起伏的绝对高度不能低于此值，低于此值的全当做平地。
- `SMOOTH_WINDOW=10`: 移动平均降噪的窗口大小。设定为 10 意味着它会提取前后各 5 个点的均值，将尖刺磨平，更平滑的曲线有利于提取真正宽广有力的共有峰。
- `MIN_WIDTH=5`: 峰的最小底宽要求。默认设定下，必须连续凸起跨越至少 5 个数据点（在您的数据中约等于 1.5MB 长度）才承认它是峰，单点的突刺会被直接忽略。
- `DISTANCE_WINDOW=3000000`: 跨样本的“共有”判定距离限度。处于同一染色体且升降方向相同的峰，如果它们横坐标间的距离在此数值（默认 3Mb）以内，则会被认定为同一批群体的**共有簇（Shared Peak Cluster）**。
- `FLANK_BP=1000000`: 提取基因的边界外延宽度（默认 1Mb）。在找到共有峰簇的最左顶点和最右顶点后，在这构成的总宽度基础上向左右两翼继续扩展 1Mb，用于检索落入此范围内的候补基因。
- `MIN_SUPPORT=2`: 所需的最少共有样本数量。用于生成过滤版的高置信基因清单（大于等于该数值）。

## 输入输出概览

* **核心输入**：各个样本群体的 `plot.tsv` 文件、基因组 `GFF3` 注释文件。
* **默认生成目录**：`/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/output/gene_share_peak/`
* 关于具体的输出产物和结果文件解读详见配套脚本的文档 `gene_share_peak.md`。
