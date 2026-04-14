#自动进行vcf计算snpindex及绘图；shell脚本
#输入，vcf文件，群体名，窗口大小，后续命名前缀为vcf_sample_window
#1.变量设计
VCF_FILE=/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz
sample=36
window=1000000
PREFIX=${VCF_FILE}_${sample}_${window}
OUTPUT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/$PREFIX
mkdir -p $OUTPUT_DIR
#2.处理逻辑
#2.1 snpindex计算
python /data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/gvcf_plidy_30_BSA.py --vcf $VCF_FILE --high-bulk F9-${sample}H   --low-bulk F9-${sample}L   --high-parent C6_Parent   --low-parent C3_Parent --out-prefix ploidy_30_3.18   --window-size $window   --step-size ${window}/10   --min-parent-dp 5   --min-bulk-dp 10   --min-window-snps 15   --allow-parent-ad-rescue
echo "snoindex ok "
#2.2 绘图
python /data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/bsa.fig.py   -f ${PREFIX}.tsv   -snpindex   -o ${VCF_FILE}_${window}.pdf
echo "fig ok"
echo "恭喜你，${VCF_FILE}_${sample}_${window} all be ok"
上面是我写的原版，下面是gpt改正过的版本，都很棒
#!/usr/bin/env bash
set -euo pipefail

# 自动进行 VCF 计算 SNP-index 及绘图；shell 脚本
# 输入：VCF 文件、群体名、窗口大小
# 输出命名：vcf_sample_window
# 已处理文件不覆盖，存在则跳过

VCF_FILE=/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz
HIGH_PARENT=C6_Parent
LOW_PARENT=C3_Parent
window=1000000

# 你的 5 个群体样本号，按需改这里
SAMPLES=(36 70 71 72 73)

# 用文件名做前缀，避免路径里的 / 把命名弄坏
VCF_BASE=$(basename "$VCF_FILE")
VCF_BASE=${VCF_BASE%.vcf.gz}
VCF_BASE=${VCF_BASE%.gz}

for sample in "${SAMPLES[@]}"; do
  PREFIX="${VCF_BASE}_${sample}_${window}"
  OUTPUT_DIR=/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/$PREFIX
  BSA_DIR="$OUTPUT_DIR/$PREFIX"

  PLOT_TSV="$BSA_DIR/$PREFIX.plot.tsv"
  PLOT_PDF="$BSA_DIR/$PREFIX.pdf"

  mkdir -p "$OUTPUT_DIR"

  if [[ -s "$PLOT_PDF" ]]; then
    echo "skip ${PREFIX} : pdf already exists"
    continue
  fi

  if [[ ! -s "$PLOT_TSV" ]]; then
    step_size=$((window / 10))
    python /data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/gvcf_plidy_30_BSA.py \
      --vcf "$VCF_FILE" \
      --high-bulk "F9-${sample}H" \
      --low-bulk "F9-${sample}L" \
      --high-parent "$HIGH_PARENT" \
      --low-parent "$LOW_PARENT" \
      --out-prefix "$PREFIX" \
      --output-root "$OUTPUT_DIR" \
      --window-size "$window" \
      --step-size "$step_size" \
      --min-parent-dp 5 \
      --min-bulk-dp 10 \
      --min-window-snps 15 \
      --allow-parent-ad-rescue
    echo "snpindex ok: ${PREFIX}"
  else
    echo "skip ${PREFIX} : plot tsv already exists"
  fi

  if [[ -s "$PLOT_TSV" && ! -s "$PLOT_PDF" ]]; then
    python /data2/chenh/bsa/CleanData/gvcf_plidy_30/scripts/bsa.fig.py \
      -f "$PLOT_TSV" \
      -snpindex \
      -o "$PLOT_PDF"
    echo "fig ok: ${PREFIX}"
  fi

  echo "恭喜你，${PREFIX} all be ok"
done
