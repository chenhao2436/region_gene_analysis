#!/usr/bin/env bash
set -euo pipefail

# 生成时间: 2026-04-10
# 脚本名称: BSA_GENE_extract_SharePeak.sh
# 脚本作用:
# 1. 统一调度基于趋势峰值提取和位置聚类的新流程 (gene_share_peak)
# 2. 从 5 个群体的 BSA plot.tsv 中寻找趋势相同的共同峰区域
# 3. 提取支持簇中的候选基因，支持后续变异检查

ROOT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract_SharePeak
SCRIPT_DIR=${ROOT_DIR}/scripts
OUTPUT_DIR=${ROOT_DIR}/output

GFF3=/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.gff3

PROMINENCE=0.05
SMOOTH_WINDOW=10
MIN_WIDTH=5
DISTANCE_WINDOW=3000000
FLANK_BP=1000000
MIN_SUPPORT=2

PLOT_FILES=(
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_36_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_70_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_168_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint_301_3000000.plot.tsv
)


mkdir -p "${OUTPUT_DIR}/gene_share_peak"

echo "[1/1] gene_share_peak (Joint Peak Calling and Clustering)"
python "${SCRIPT_DIR}/gene_share_peak" \
  --plot-tsv "${PLOT_FILES[@]}" \
  --gff3 "${GFF3}" \
  --output-dir "${OUTPUT_DIR}/gene_share_peak" \
  --prominence "${PROMINENCE}" \
  --smooth-window "${SMOOTH_WINDOW}" \
  --min-width "${MIN_WIDTH}" \
  --distance-window "${DISTANCE_WINDOW}" \
  --flank "${FLANK_BP}" \
  --min-support "${MIN_SUPPORT}"

echo "================================================================="
echo "[DONE] BSA_GENE_extract_SharePeak.sh finished extracting shared peaks."
echo "You can check the summary of shared peaks at:"
echo "${OUTPUT_DIR}/gene_share_peak/share_peak.all_clusters.tsv"
echo "Filtered bed and genes for extraction are also inside this folder."
echo "================================================================="
