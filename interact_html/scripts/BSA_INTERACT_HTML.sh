#!/usr/bin/env bash
set -euo pipefail

# 生成时间: 2026-04-13
# 脚本名称: BSA_INTERACT_HTML.sh
# 脚本作用:
#   调用 interact_html，把 bsa.fig.py -snpindex 风格的静态 BSA 图转换为交互式 HTML。
#   本版不绘制 ED threshold 红线，不绘制 Share Peak 背景块。
# 输入文件:
#   PLOT_FILES 数组中的 6 个 plot.tsv/snpindex 文件。
# 输出文件:
#   ${OUTPUT_DIR}/bsa_interactive_combined.html
#   ${OUTPUT_DIR}/bsa_interactive_subplots.html
#   ${OUTPUT_DIR}/bsa_interactive_301.html
#   ${OUTPUT_DIR}/bsa_interactive_168.html
#   ${OUTPUT_DIR}/bsa_interactive_184.html
#   ${OUTPUT_DIR}/bsa_interactive_36.html
#   ${OUTPUT_DIR}/bsa_interactive_70.html
#   ${OUTPUT_DIR}/bsa_interactive_F2_Canz_3MB.html
# 使用方法:
#   bash /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html/scripts/BSA_INTERACT_HTML.sh

ROOT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html
SCRIPT_DIR=${ROOT_DIR}/scripts
OUTPUT_DIR=${ROOT_DIR}/output

PLOT_FILES=(
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/301_new_genome_joint_F2_3000000.plot.tsv
)

SAMPLE_LABELS=(
  301
  168
  184
  36
  70
  F2_Canz_3MB
)

mkdir -p "${OUTPUT_DIR}"

echo "[START] Generating BSA interactive HTML plots..."
python "${SCRIPT_DIR}/interact_html" \
  --plot-tsv "${PLOT_FILES[@]}" \
  --sample-labels "${SAMPLE_LABELS[@]}" \
  --output-dir "${OUTPUT_DIR}"

echo "================================================================="
echo "[DONE] BSA_INTERACT_HTML.sh finished."
echo "Your interactive plots are located at:"
echo "1. ${OUTPUT_DIR}/bsa_interactive_combined.html"
echo "2. ${OUTPUT_DIR}/bsa_interactive_subplots.html"
echo "3. ${OUTPUT_DIR}/bsa_interactive_301.html"
echo "4. ${OUTPUT_DIR}/bsa_interactive_168.html"
echo "5. ${OUTPUT_DIR}/bsa_interactive_184.html"
echo "6. ${OUTPUT_DIR}/bsa_interactive_36.html"
echo "7. ${OUTPUT_DIR}/bsa_interactive_70.html"
echo "8. ${OUTPUT_DIR}/bsa_interactive_F2_Canz_3MB.html"
echo "================================================================="
exit 0

# 生成时间: 2026-04-13
# 脚本名称: BSA_INTERACT_HTML.sh
# 脚本作用:
# 调用 interact_html 生成基于 Plotly 的交互式网页（矢量可视化）
# 提供了全部 6 个输入 plot.tsv 样本，并调用前置生成的共有峰文件作为图背景层。

ROOT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/interact_html
SCRIPT_DIR=${ROOT_DIR}/scripts
OUTPUT_DIR=${ROOT_DIR}/output

CLUSTER_TSV=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/output/gene_share_peak_smooth/share_peak.all_clusters.tsv
MIN_CLUSTER_SUPPORT=2

PLOT_FILES=(
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak/plot_tsv/301_new_genome_joint_F2_3000000.plot.tsv
)

mkdir -p "${OUTPUT_DIR}"

echo "[START] Generating interactive HTML plots..."
python "${SCRIPT_DIR}/interact_html" \
  --plot-tsv "${PLOT_FILES[@]}" \
  --output-dir "${OUTPUT_DIR}" \
  --cluster-tsv "${CLUSTER_TSV}" \
  --min-cluster-support "${MIN_CLUSTER_SUPPORT}"

echo "================================================================="
echo "[DONE] BSA_INTERACT_HTML.sh finished."
echo "Your interactive plots are located at:"
echo "1. ${OUTPUT_DIR}/bsa_interactive_subplots.html"
echo "2. ${OUTPUT_DIR}/bsa_interactive_combined.html"
echo "================================================================="
