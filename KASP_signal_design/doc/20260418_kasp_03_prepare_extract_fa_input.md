# 20260418_kasp_03_prepare_extract_fa_input.py

## 运行场景
用于把通过过滤的位点整理成 `extract_fa.py` 能直接读取的四列表, 同时输出最终推荐标记清单。

## 脚本目标
- 左右两侧分别整理最终入选位点
- 生成无表头四列表供 `extract_fa.py` 使用
- 输出推荐清单并解释为什么不足

## 输入文件
- `--pass-sites`: `02_pass_sites.tsv`
- `--checked-sites`: `02_checked_sites.tsv`

## 输出文件
- `--output-extract-input`: `03_extract_fa_input.tsv`
- `--output-recommended`: `03_recommended_markers.tsv`

## 参数解释
- `--target-markers-per-side`: 每侧目标数量, 默认 2

## 使用示例
```bash
python 20260418_kasp_03_prepare_extract_fa_input.py \
  --pass-sites 02_pass_sites.tsv \
  --checked-sites 02_checked_sites.tsv \
  --output-extract-input 03_extract_fa_input.tsv \
  --output-recommended 03_recommended_markers.tsv
```

## 结果表字段说明
`03_extract_fa_input.tsv`
- 无表头
- 每行 4 列: `chrom pos ref alt`
- 直接给 `extract_fa.py --input` 使用

`03_recommended_markers.tsv`
- `side`: 左右侧
- `selected_count`: 当前侧最终选中的位点数
- `target_count`: 当前侧目标位点数
- `checked_candidates`: 当前侧已检查候选数
- `pass_count`: 当前侧通过过滤总数
- `shortage_reason`: 若不足目标数量, 这里写明原因统计
- `chrom/pos/ref/alt`: 最终入选位点信息
- `delta_snp_index`: 排序依据
- `checked_rank`: 当前侧排序名次
- `parent1_gt/parent2_gt/high_pool_gt/low_pool_gt`: 基因型信息

## 常见报错与排查
- `通过位点文件不存在`: 第 2 步未跑完或路径错误
- `检查明细文件不存在`: 第 2 步明细文件缺失
- `03_extract_fa_input.tsv` 为空: 说明没有任何位点通过硬过滤
- `shortage_reason` 显示多个失败原因: 说明某一侧需要扩大候选或复核区间
