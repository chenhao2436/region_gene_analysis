# interact_html 说明

生成时间: 2026-04-13

## 脚本作用

`interact_html` 用来把 `bsa.fig.py -snpindex` 风格的 BSA 静态图转换成交互式 HTML。它读取多个 `plot.tsv/snpindex` 文件，按原图逻辑把每条染色体的位置累加成全基因组横坐标，并用 Plotly 输出可缩放、可悬停查看点位信息的 HTML 图。

本版按需求做了以下限制：

1. 不使用 ED 作图逻辑。
2. 不绘制 threshold 红线。
3. 不绘制 Share Peak 浅红色背景块。
4. 主要保留原图排版，只增加交互悬停信息和必要图例/色条。

## 输入文件

调度脚本中的 `PLOT_FILES` 数组包含 6 个输入文件：

```bash
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000.plot.tsv
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000.plot.tsv
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000.plot.tsv
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000.plot.tsv
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000.plot.tsv
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/301_new_genome_joint_F2_3000000.plot.tsv
```

脚本读取的核心列如下：

| 列号 | 含义 |
|---|---|
| 0 | 染色体名 |
| 1 | 窗口中心位置 bp |
| 5 | DELTA/SNP-index 差值，作为 y 轴 |
| 6 或末列 | SNP 数量，作为点颜色 |

## 输出文件

输出目录为：

```bash
/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/output
```

生成文件如下：

```bash
bsa_interactive_combined.html
bsa_interactive_subplots.html
bsa_interactive_301.html
bsa_interactive_168.html
bsa_interactive_184.html
bsa_interactive_36.html
bsa_interactive_70.html
bsa_interactive_F2_Canz_3MB.html
```

## 运行方法

推荐直接运行调度脚本：

```bash
bash /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/scripts/BSA_INTERACT_HTML.sh
```

如果手动运行 Python 脚本，示例为：

```bash
python /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/scripts/interact_html \
  --plot-tsv a.plot.tsv b.plot.tsv \
  --sample-labels 301 168 \
  --output-dir /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/output
```

## 图形逻辑

1. 染色体按编号排序，例如 `Chr01`、`Chr02`、`Chr10`。
2. 每条染色体使用该染色体最后一个窗口位置作为长度，再累加到全基因组 x 轴。
3. 点的 x 轴是 `染色体累计起点 + 窗口中心位置`。
4. 点的 y 轴是第 5 列 DELTA/SNP-index 差值。
5. 点颜色来自 SNP 数量，并按每个样本的 90% 分位数截断后映射到 `plasma_r` 色阶。
6. 鼠标悬停点位时显示样本名、染色体、位置、DELTA 值、SNP 数量。

## 环境依赖

需要 Python 环境中安装：

```bash
plotly
numpy
```

---

以下为旧版说明文本，以上方新版说明为准。

# interact_html 说明

生成时间: 2026-04-13

## 原理与作用

本脚本是对之前 BSA 静态出图的一个**颠覆性升级**。它直接读取原始数据的 `plot.tsv` 文件，并使用 Plotly 库将其转变为一个**高性能、全矢量的交互式 HTML 网页**。
不再需要对着静态图片发愁、无法查看具体坐标与数值了——您可以直接在浏览器双击打开生成的 HTML，鼠标点到哪里，就能看到哪里的精细信息。

它具有如下特性：
1. **纯净折线**：过滤了多余的背景参考线（HIGH/LOW等），只为您画出核心起伏线（`DELTA_SNP_INDEX`）。
2. **多维色阶映射**：继承了原有静态图最优秀的特点，利用点的颜色（从黄到蓝紫）展示该窗口蕴含的 `SNP_N` 个数，辅助判断该位点的可信度。
3. **全维基因组图谱**：将各染色体首尾相接拼成一个全基因组长轴，并在轴中间标记染色体名，在边界画辅助分割线，让染色体全貌一览无余。
4. **共有峰区间自适应高亮**：它会自动读取跑出来的 `share_peak.all_clusters.tsv`。所有被机器判定为包含共有峰的区域，会**在底层加上一层浅红色的半透明色带**。如果你觉得这块红带子挡住你视线了，只需点击右侧图例即可将它**关闭消失**！

## 使用方法

推荐直接运行已经包装好的 `BSA_INTERACT_HTML.sh`：

```bash
bash /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/scripts/BSA_INTERACT_HTML.sh
```

**环境依赖**：需要 Python 中已安装以下库：
- `plotly` （如果报错请执行 `pip install plotly`）
- `numpy`

## 输出结果解读

脚本会聪明地在目录 `output/interact_html/output` 下同时生成两个文件，满足您的对比需要：

### 1. `bsa_interactive_subplots.html`（子图版）
如同您平时作图排版一样，6 个样本（包含 `301_new`）自上而下分成 6 行。
这种做法的优点是：**每个样本线互不干扰**。并且当你在上方图按住左键框选放大某一 QTL 峰区间时，下面 5 个子图会自动跟着**同步放大对齐**在这个相同染色体区间上，比肉眼核对静态图强百倍。

### 2. `bsa_interactive_combined.html`（重叠版）
将 6 个群体的图压在同一个坐标框架内。因为颜色有深浅遮挡，因此设计了交互式图例：
**点击图例条目**即可单独关掉或打开某个样本的显示，非常适合找出那几条走势有异议的关键群体验算差异。

## 交互微操指南
- **鼠标悬停**：移到点上，光标旁出现黑框，记录了（所属群体样本、所属染色体、精准 Position（以加逗号格式，如 5,400,000）、纵坐标 Delta 值，甚至该点包含几个 SNP）。
- **左键框选**：按住左键拉小框，松开即可钻入微观区域。
- **拖拽平移**：点击图顶上工具条的四向箭头（Pan），左键摁住图可以直接左右滑动。
- **双击复原**：双击图内任何空白区，镜头自动拉回到总览尺寸。
- **共有峰区段控制**：在右侧图例点击名为 `Share Peak (>= 2)` 的条目，图上浅红色候选背景框就会消失，再点一次重现，并且悬停在背景框边缘也能直接拿到该候选区段的 Cluster 号。
