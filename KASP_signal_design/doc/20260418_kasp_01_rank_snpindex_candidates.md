# 20260418_kasp_01_rank_snpindex_candidates.py

## 运行场景
用于从一个给定峰区中, 先提取左右两个向内收的 1 Mb 区间; 如果某一侧候选不足, 再按固定步长向峰区中心扩宽, 并在 snpindex 文件里筛出 `DELTA_SNP_INDEX > 0` 的候选位点排序。

## 脚本目标
- 解析 `Chr05:248.7-255.3MB` 这类区间
- 先拆成左右两个 1 Mb 区间, 不足时继续向峰区中心扩宽
- 自动识别 `CHROM/POS/DELTA_SNP_INDEX`
- 输出结构化候选表和汇总表

## 输入文件
- `--snpindex`: 上游 snpindex 表

## 输出文件
- `--output-candidates`: `01_ranked_candidates.tsv`
- `--output-summary`: `01_rank_summary.tsv`

## 参数解释
- `--region`: 峰区区间
- `--snpindex`: snpindex 文件路径
- `--output-candidates`: 候选表输出路径
- `--output-summary`: 汇总表输出路径
- `--initial-rank-limit`: 默认 20, 仅用于汇总中记录前筛规模
- `--expansion-step-bp`: 向峰区中心扩宽的步长, 默认 1000000 bp

## 使用示例
```bash
python 20260418_kasp_01_rank_snpindex_candidates.py \
  --region Chr05:248.7-255.3MB \
  --snpindex sample.plot.tsv \
  --output-candidates 01_ranked_candidates.tsv \
  --output-summary 01_rank_summary.tsv
```

## 结果表字段说明
`01_ranked_candidates.tsv`
- `side`: `left` 或 `right`
- `chrom`: 染色体名
- `pos`: 位点坐标
- `delta_snp_index`: DELTA_SNP_INDEX 数值
- `rank`: 该侧的降序排名
- `window_start`: 当前侧窗口起点
- `window_end`: 当前侧窗口终点
- `expansion_level`: 扩宽层级, `0` 表示初始 1 Mb 窗口
- `source_snpindex`: 来源文件

`01_rank_summary.tsv`
- `side`: 左右侧
- `chrom`: 染色体
- `initial_window_start/initial_window_end`: 初始 1 Mb 区间
- `max_window_start/max_window_end`: 该侧可扩宽到的最大范围
- `total_candidates`: 候选总数
- `initial_rank_limit`: 前筛目标数量
- `initial_slice_count`: 实际能纳入前筛的数量
- `expansion_levels`: 该侧共使用了几层扩宽
- `max_delta_snp_index/min_delta_snp_index`: 当前侧最大和最小 DELTA
- `source_snpindex`: 来源文件

## 常见报错与排查
- `无法解析区间`: 区间格式不是 `Chr:248.7-255.3MB`
- `输入区间长度必须至少为 2 Mb`: 左右各 1 Mb 无法从更短区间中拆出
- `表头缺少 DELTA_SNP_INDEX 列`: 请检查列名
- `无表头 snpindex 文件列数不足`: 无表头时默认按第 1/2/6 列读取
