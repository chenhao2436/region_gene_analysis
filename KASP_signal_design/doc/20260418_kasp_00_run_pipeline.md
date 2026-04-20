# 20260418_kasp_00_run_pipeline.sh

## 运行场景
用于在 Linux 服务器上一条命令跑完整个 KASP 候选筛选流程, 从峰区区间开始, 一直跑到 `extract_fa.py` 生成最终序列。

## 脚本目标
- 统一读取 `20260418_kasp_pipeline.ini`
- 接受命令行覆盖关键参数
- 串联执行 `01 -> 02 -> 03 -> extract_fa.py`
- 固定输出结构化结果文件和日志文件

## 输入文件
- 配置文件: `20260418_kasp_pipeline.ini`
- snpindex 文件: `snpindex_path`
- VCF 文件: `vcf_path`
- 参考基因组: `reference_fasta`
- 序列提取脚本: `scripts/extract_fa.py`

## 输出文件
- `01_ranked_candidates.tsv`
- `01_rank_summary.tsv`
- `02_checked_sites.tsv`
- `02_pass_sites.tsv`
- `03_extract_fa_input.tsv`
- `03_recommended_markers.tsv`
- `04_kasp_sequences.tsv`
- `pipeline.log`

## 参数解释
- `--config`: 指定配置文件路径
- `--region`: 峰区范围, 例如 `Chr05:248.7-255.3MB`
- `--snpindex`: 临时覆盖 `snpindex_path`
- `--high-pool`: 临时覆盖高池样本名
- `--low-pool`: 临时覆盖低池样本名
- `--outdir`: 指定输出目录
- 配置项 `expansion_step_bp`: 某一侧不足目标数时, 向峰区中心扩宽的步长, 默认 `1000000`

## 使用示例
```bash
bash 20260418_kasp_00_run_pipeline.sh doc/config/20260418_kasp_pipeline.ini
```

```bash
bash 20260418_kasp_00_run_pipeline.sh --config doc/config/20260418_kasp_pipeline.ini
```

```bash
bash 20260418_kasp_00_run_pipeline.sh \
  --config doc/config/20260418_kasp_pipeline.ini \
  --region Chr05:248.7-255.3MB \
  --snpindex /path/to/sample.plot.tsv \
  --high-pool Hunchi_H \
  --low-pool Hunchi_L \
  --outdir /data2/chenh/bsa/F2:3_WORKSPACE/bsa_HunchiData/BSA_out/kasp_signal_design/results/Chr05_248.7_255.3MB
```

## 结果表说明
- `01_ranked_candidates.tsv`: 左右两侧按 `DELTA_SNP_INDEX` 排序的候选位点
- 若某一侧初始 1 Mb 区间不足, 会按 `expansion_step_bp` 继续向峰区中心扩宽后再补选
- `02_checked_sites.tsv`: 每个已检查位点的过滤状态和原因
- `02_pass_sites.tsv`: 真正通过所有硬过滤条件的位点
- `03_extract_fa_input.tsv`: 给 `extract_fa.py` 用的四列表, 无表头
- `03_recommended_markers.tsv`: 左右两侧最终推荐位点清单, 以及不足原因
- `04_kasp_sequences.tsv`: `extract_fa.py` 生成的最终 KASP 序列表

## 常见报错与排查
- `配置文件不存在`: 配置文件路径不对, 或传入的是相对路径但当前目录不对
- `region 未设置`: 配置文件和命令行都没给 `region`
- `snpindex_path 未设置`: 需要在配置文件填写或用 `--snpindex` 覆盖
- `high_pool_sample 和 low_pool_sample 必须提供`: 高低池样本名还没填
- `VCF 缺少索引文件`: 请确认 `.tbi` 或 `.csi` 存在
- `缺少命令: bcftools/tabix/python`: 运行环境依赖不完整
- `extract_fa.py 不存在`: 请确认 `scripts/extract_fa.py` 位于预期目录
