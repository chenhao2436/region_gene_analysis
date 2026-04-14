#!/usr/bin/env bash
set -euo pipefail

# 生成时间: 2026-03-26
# 脚本名称: BSA_GENE_extract.sh
# 脚本作用:
# 1. 统一调度 gene_extract / gene_common / gene_region_vcf / gene_out
# 2. 从 5 个群体的 BSA plot.tsv 中提取显著信号窗口对应基因
# 3. 生成至少 2 个群体支持和 5 个群体全部支持两套交集基因结果
# 4. 对交集基因本身范围内的所有变异做 Annovar 注释
# 5. 最终整合出基因支持情况、拟南芥功能、变异注释三类结果

ROOT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract
SCRIPT_DIR=${ROOT_DIR}/scripts
OUTPUT_DIR=${ROOT_DIR}/output

GTF_FILE=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract/output/gene_region_vcf/at_least_2/Canz.genome.gtf
GFF3=/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.gff3
REFERENCE_FASTA=/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa
JOINT_VCF=/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz
RBH_FILE=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/gene_extract/AT_Canz_reciprocal_best_hits.txt
TAIR_ANNOTATION=/data2/chenh/genome/TAIR_gene.annotate.unique_v1.2

DELTA_THRESHOLD=0.4
FLANK_BP=1000000
MIN_SUPPORT=2
ALL_SUPPORT=5

PLOT_FILES=(
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_36_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_36_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_36_1000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_70_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_70_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_70_1000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_168_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_168_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_168_1000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_184_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_1000000.plot.tsv
  /data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_301_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_301_1000000/hunchi_and_qinben.prepared_single_sample.genome.joint_301_1000000.plot.tsv
)

mkdir -p \
  "${OUTPUT_DIR}/gene_extract" \
  "${OUTPUT_DIR}/gene_common" \
  "${OUTPUT_DIR}/gene_region_vcf/at_least_2" \
  "${OUTPUT_DIR}/gene_region_vcf/all_5" \
  "${OUTPUT_DIR}/gene_out"

echo "[1/4] gene_extract"
python "${SCRIPT_DIR}/gene_extract" \
  --plot-tsv "${PLOT_FILES[@]}" \
  --gff3 "${GFF3}" \
  --output-dir "${OUTPUT_DIR}/gene_extract" \
  --delta-threshold "${DELTA_THRESHOLD}" \
  --flank "${FLANK_BP}"

echo "[2/4] gene_common"
python "${SCRIPT_DIR}/gene_common" \
  --sample-candidate-tsv "${OUTPUT_DIR}/gene_extract/gene_extract.sample_candidate_genes.tsv" \
  --output-dir "${OUTPUT_DIR}/gene_common" \
  --min-support "${MIN_SUPPORT}" \
  --all-support "${ALL_SUPPORT}"

echo "[3/4] gene_region_vcf : at_least_2"
python "${SCRIPT_DIR}/gene_region_vcf" \
  --gene-bed "${OUTPUT_DIR}/gene_common/gene_common.at_least_2.bed" \
  --gene-table "${OUTPUT_DIR}/gene_common/gene_common.at_least_2.tsv" \
  --vcf "${JOINT_VCF}" \
  --reference-fasta "${REFERENCE_FASTA}" \
  --gtf "${GTF_FILE}" \
  --label at_least_2 \
  --output-dir "${OUTPUT_DIR}/gene_region_vcf/at_least_2"

echo "[3/4] gene_region_vcf : all_5"
python "${SCRIPT_DIR}/gene_region_vcf" \
  --gene-bed "${OUTPUT_DIR}/gene_common/gene_common.all_5.bed" \
  --gene-table "${OUTPUT_DIR}/gene_common/gene_common.all_5.tsv" \
  --vcf "${JOINT_VCF}" \
  --reference-fasta "${REFERENCE_FASTA}" \
  --gtf "${GTF_FILE}" \
  --label all_5 \
  --output-dir "${OUTPUT_DIR}/gene_region_vcf/all_5"

echo "[4/4] gene_out"
python "${SCRIPT_DIR}/gene_out" \
  --atleast2-table "${OUTPUT_DIR}/gene_common/gene_common.at_least_2.tsv" \
  --all5-table "${OUTPUT_DIR}/gene_common/gene_common.all_5.tsv" \
  --rbh-file "${RBH_FILE}" \
  --tair-annotation "${TAIR_ANNOTATION}" \
  --annovar-prefix-atleast2 "${OUTPUT_DIR}/gene_region_vcf/at_least_2/gene_region_vcf.at_least_2" \
  --annovar-prefix-all5 "${OUTPUT_DIR}/gene_region_vcf/all_5/gene_region_vcf.all_5" \
  --output-dir "${OUTPUT_DIR}/gene_out"

echo "[DONE] BSA_GENE_extract.sh finished"
