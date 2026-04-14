#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

###############################################################################
# 脚本名称:
#   prepare_single_sample_gvcf_list_gpt_3m8d.sh
#
# 脚本目标:
#   读取一个 joint-calling 用的 gVCF 列表。
#   如果列表中某个 gVCF 含有多个样本，则把它拆成多个“单样本 gVCF”，
#   并生成一个新的 prepared.gvcf.list，供 GenomicsDBImport 使用。
#
# 适用场景:
#   你当前的 Parents_Combined.g.vcf.gz 含有两个亲本样本，而 GATK 4.1.8.1 的
#   GenomicsDBImport 要求输入为单样本 gVCF，因此必须先拆。
#
# 工作流程:
#   1. 审计列表中每个 gVCF 的样本数
#   2. 单样本 gVCF: 原样写入新的 prepared list
#   3. 多样本 gVCF: 使用 GATK SelectVariants 按 sample 拆分
#   4. 对拆分后的每个单样本 gVCF 建索引并做 ValidateVariants -gvcf 检查
#   5. 生成新的 prepared.gvcf.list
#
# 默认输出目录:
#   /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d
#
# 用法:
#   bash prepare_single_sample_gvcf_list_gpt_3m8d.sh \
#       /data2/chenh/bsa/CleanData/vcf/hunchi_and_qinben.gvcf.list \
#       /data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa \
#       /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d
###############################################################################

LIST_FILE="${1:-/data2/chenh/bsa/CleanData/vcf/hunchi_and_qinben.gvcf.list}"
REFERENCE="${2:-/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa}"
OUTDIR="${3:-/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d}"

REFERENCE="${REFERENCE%/}"
OUTDIR="${OUTDIR%/}"

RUN_PREFIX="$(basename "$LIST_FILE")"
RUN_PREFIX="${RUN_PREFIX%.gvcf.list}"
RUN_PREFIX="${RUN_PREFIX%.list}"

LOG_DIR="${OUTDIR}/logs"
SPLIT_DIR="${OUTDIR}/split_gvcf"
TMP_DIR="${OUTDIR}/tmp"
CLEAN_LIST="${OUTDIR}/${RUN_PREFIX}.cleaned.gvcf.list"
AUDIT_TSV="${OUTDIR}/${RUN_PREFIX}.sample_audit.tsv"
PREPARED_LIST="${OUTDIR}/${RUN_PREFIX}.prepared_single_sample.gvcf.list"
MASTER_LOG="${LOG_DIR}/${RUN_PREFIX}.prepare_$(date '+%Y%m%d_%H%M%S').log"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "未找到命令: $1"
}

extract_samples_from_gvcf() {
    local gvcf="$1"
    local header_line

    set +o pipefail
    header_line="$(gzip -dc "$gvcf" | awk 'BEGIN{FS="\t"} /^#CHROM/ {print; exit}')"
    local status=$?
    set -o pipefail

    [[ $status -eq 0 ]] || return 1
    [[ -n "$header_line" ]] || return 1

    awk -F'\t' '
        {
            if (NF <= 9) exit 1
            for (i = 10; i <= NF; i++) print $i
        }
    ' <<<"$header_line"
}

sanitize_name() {
    printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

mkdir -p "$OUTDIR" "$LOG_DIR" "$SPLIT_DIR" "$TMP_DIR"
exec > >(tee -a "$MASTER_LOG") 2>&1

need_cmd bash
need_cmd gatk
need_cmd gzip
need_cmd awk
need_cmd sed
need_cmd tee

[[ -f "$LIST_FILE" ]] || die "列表文件不存在: $LIST_FILE"
[[ -f "$REFERENCE" ]] || die "参考基因组不存在: $REFERENCE"

log "开始准备单样本 gVCF 列表"
log "LIST_FILE   = $LIST_FILE"
log "REFERENCE   = $REFERENCE"
log "OUTDIR      = $OUTDIR"
log "MASTER_LOG  = $MASTER_LOG"

awk 'NF && $1 !~ /^#/' "$LIST_FILE" > "$CLEAN_LIST"
[[ -s "$CLEAN_LIST" ]] || die "清理后列表为空: $CLEAN_LIST"

printf "gvcf_path\tsample_count\tsample_names\taction\n" > "$AUDIT_TSV"
: > "$PREPARED_LIST"

while IFS= read -r gvcf; do
    [[ -f "$gvcf" ]] || die "gVCF 文件不存在: $gvcf"
    [[ -f "${gvcf}.tbi" || -f "${gvcf}.csi" || -f "${gvcf}.idx" ]] || die "gVCF 缺少索引: $gvcf"

    mapfile -t samples < <(extract_samples_from_gvcf "$gvcf")
    [[ ${#samples[@]} -gt 0 ]] || die "无法读取样本名: $gvcf"

    sample_csv="$(IFS=,; echo "${samples[*]}")"

    if [[ ${#samples[@]} -eq 1 ]]; then
        printf "%s\t%d\t%s\tkeep\n" "$gvcf" "${#samples[@]}" "$sample_csv" >> "$AUDIT_TSV"
        printf "%s\n" "$gvcf" >> "$PREPARED_LIST"
        log "保留单样本 gVCF: $gvcf"
        continue
    fi

    printf "%s\t%d\t%s\tsplit\n" "$gvcf" "${#samples[@]}" "$sample_csv" >> "$AUDIT_TSV"
    log "拆分多样本 gVCF: $gvcf"
    log "  样本名: $sample_csv"

    input_base="$(basename "$gvcf")"
    input_base="${input_base%.g.vcf.gz}"
    input_base="${input_base%.vcf.gz}"

    for sample_name in "${samples[@]}"; do
        safe_sample="$(sanitize_name "$sample_name")"
        out_gvcf="${SPLIT_DIR}/${input_base}.${safe_sample}.g.vcf.gz"
        sample_tmp="${TMP_DIR}/${input_base}.${safe_sample}"
        mkdir -p "$sample_tmp"

        log "  开始拆分样本: $sample_name"

        gatk \
            --java-options "-Xms8g -Xmx8g -Djava.io.tmpdir=${sample_tmp}" \
            SelectVariants \
            -R "$REFERENCE" \
            -V "$gvcf" \
            --sample-name "$sample_name" \
            -O "$out_gvcf" \
            --tmp-dir "$sample_tmp"

        gatk IndexFeatureFile -I "$out_gvcf"

        gatk \
            --java-options "-Xms4g -Xmx4g -Djava.io.tmpdir=${sample_tmp}" \
            ValidateVariants \
            -R "$REFERENCE" \
            -V "$out_gvcf" \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --tmp-dir "$sample_tmp"

        printf "%s\n" "$out_gvcf" >> "$PREPARED_LIST"
        rm -rf "$sample_tmp"
        log "  完成样本: $sample_name -> $out_gvcf"
    done
done < "$CLEAN_LIST"

log "样本审计表: $AUDIT_TSV"
log "新的单样本 gVCF 列表: $PREPARED_LIST"
log "准备完成"
