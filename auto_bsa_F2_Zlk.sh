#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 自动 BSA 分析脚本 —— 使用 F2_BSA.py（@zlk）
# 根据 VCF 表头按样本名自动匹配列号，减少人工配错率
# =====================================================================

# ======= 文件路径配置 =======
VCF_FILE=/data2/chenh/bsa/CleanData/bsa_HunchiData/hunchi_vcf/hunchi_and_qinben.prepared_single_sample.cleaned.genome.joint.vcf.gz
BSA_SCRIPT=/data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/F2_BSA.py.py
FIG_SCRIPT=/data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/bsa.fig.py

# ======= 窗口大小 =======
window=1000000

# ======= 亲本名称配置 =======
# p1 (parent1) 对应 high-parent，b1 对应 high-bulk（H 混池，表现 p1 性状）
# p2 (parent2) 对应 low-parent， b2 对应 low-bulk （L 混池，表现 p2 性状）
HIGH_PARENT=C6_Parent
LOW_PARENT=C3_Parent

# ======= 样本编号（VCF 中的 F9-{sample}H / F9-{sample}L）=======
SAMPLES=(168 184 301 36 70)

# ======= 深度过滤参数 =======
MIN_PARENT_DP=5     # 亲本最小深度 (-d1)
MIN_BULK_DP=10      # 混池最小深度 (-d2)
MAX_PARENT_DP=80    # 亲本最大深度 (-d3)
MAX_BULK_DP=200     # 混池最大深度 (-d4)

# =====================================================================
# 从 VCF 表头自动解析样本列号（0-indexed）
# =====================================================================
echo "正在读取 VCF 表头..."
HEADER=$(bcftools view -h "$VCF_FILE" | grep "^#CHROM")
IFS=$'\t' read -ra COLS <<< "$HEADER"

get_col_index() {
    local name=$1
    for i in "${!COLS[@]}"; do
        if [[ "${COLS[$i]}" == "$name" ]]; then
            echo "$i"
            return
        fi
    done
    echo "ERROR: 样本 '${name}' 在 VCF 表头中未找到！请检查样本名。" >&2
    echo "VCF 中可用样本: ${COLS[*]:9}" >&2
    exit 1
}

# 解析亲本列号
P1_COL=$(get_col_index "$HIGH_PARENT")
P2_COL=$(get_col_index "$LOW_PARENT")

VCF_BASE="301_new_genome_joint"

OUTPUT_ROOT=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output

echo ""
echo "============================================"
echo "  F2_BSA 自动分析"
echo "============================================"
echo "VCF 文件 : $VCF_FILE"
echo "窗口大小 : $window"
echo "亲本深度 : ${MIN_PARENT_DP}-${MAX_PARENT_DP}x"
echo "混池深度 : ${MIN_BULK_DP}-${MAX_BULK_DP}x"
echo "============================================"
echo ""

for sample in "${SAMPLES[@]}"; do
  HIGH_BULK="F9-${sample}H"
  LOW_BULK="F9-${sample}L"

  # 按样本名自动匹配列号
  B1_COL=$(get_col_index "$HIGH_BULK")
  B2_COL=$(get_col_index "$LOW_BULK")

  # ---- 打印亲本-混池映射关系 ----
  echo "============================================"
  echo "  样本组: ${sample}"
  echo "--------------------------------------------"
  echo "  亲本1 (p1) : ${HIGH_PARENT}   列号=${P1_COL}"
  echo "  混池1 (b1) : ${HIGH_BULK}     列号=${B1_COL}"
  echo "  亲本2 (p2) : ${LOW_PARENT}    列号=${P2_COL}"
  echo "  混池2 (b2) : ${LOW_BULK}      列号=${B2_COL}"
  echo "============================================"

  PREFIX="${VCF_BASE}_${sample}_${window}"
  OUTPUT_DIR="${OUTPUT_ROOT}/${PREFIX}"

  mkdir -p "$OUTPUT_DIR"

  OUTFILE="${OUTPUT_DIR}/${PREFIX}"
  PLOT_PDF="${OUTPUT_DIR}/${PREFIX}.pdf"

  # 跳过已完成的样本
  if [[ -s "$PLOT_PDF" ]]; then
    echo "[跳过] ${PREFIX} : PDF 已存在"
    echo ""
    continue
  fi

  # 2.1 使用 F2_BSA.py 计算 snpindex
  if [[ ! -s "${OUTFILE}" ]]; then
    echo "[运行] F2_BSA.py snpindex 计算中..."
    python "$BSA_SCRIPT" \
      -snpindex \
      -vcf "$VCF_FILE" \
      -p1 "$P1_COL" \
      -p2 "$P2_COL" \
      -b1 "$B1_COL" \
      -b2 "$B2_COL" \
      -d1 "$MIN_PARENT_DP" \
      -d2 "$MIN_BULK_DP" \
      -d3 "$MAX_PARENT_DP" \
      -d4 "$MAX_BULK_DP" \
      -w "$window" \
      -o "$OUTFILE"
    echo "[完成] snpindex ok: ${PREFIX}"
  else
    echo "[跳过] ${PREFIX} : snpindex 结果已存在"
  fi

  # 2.2 绘图
  if [[ -s "${OUTFILE}" && ! -s "$PLOT_PDF" ]]; then
    echo "[运行] 绘图中..."
    python "$FIG_SCRIPT" \
      -f "${OUTFILE}" \
      -snpindex \
      -o "$PLOT_PDF"
    echo "[完成] fig ok: ${PREFIX}"
  fi

  echo "恭喜你，${PREFIX} all be ok"
  echo ""
done

echo "============================================"
echo "  所有样本处理完毕！"
echo "============================================"
