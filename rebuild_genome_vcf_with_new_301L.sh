#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

###############################################################################
# 脚本名称:
#   rebuild_genome_vcf_with_new_301L.sh
#
# 脚本目标:
#   在不改变最终文件名和正式输出路径的前提下，用新的 301L gVCF 替换旧输入，
#   重新执行全基因组、全样本 joint-calling，保证最终 VCF 的整体一致性。
#
# 设计原则:
#   1. 只替换 301L 的上游 gVCF 输入，其余 11 个输入保持不变；
#   2. joint-calling 仍按原流程执行:
#         GenomicsDBImport -> GenotypeGVCFs -> GatherVcfs
#   3. 所有中间文件都写入临时运行目录；
#   4. 只有当临时最终 VCF 通过验证后，才覆盖正式结果；
#   5. 覆盖完成后删除临时目录，只保留正式最终 VCF 和索引。
#
# 默认关键路径:
#   主脚本:
#     /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/gvcf2vcf_genome_gpt_3m6d.sh
#   基线输入 list:
#     /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.cleaned.gvcf.list
#   新 301L:
#     /data2/chenh/bsa/CleanData/301_check/gvcf/F9-301L.g.vcf.gz
#   正式最终 VCF:
#     /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz
#
# 用法:
#   bash rebuild_genome_vcf_with_new_301L.sh
#
#   bash rebuild_genome_vcf_with_new_301L.sh \
#       /data2/chenh/bsa/CleanData/301_check/gvcf/F9-301L.g.vcf.gz
#
#   bash rebuild_genome_vcf_with_new_301L.sh \
#       /data2/chenh/bsa/CleanData/301_check/gvcf/F9-301L.g.vcf.gz \
#       /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.cleaned.gvcf.list
###############################################################################

MAIN_SCRIPT="${MAIN_SCRIPT:-/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/gvcf2vcf_genome_gpt_3m6d.sh}"
BASE_LIST="${BASE_LIST:-/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/hunchi_and_qinben.prepared_single_sample.cleaned.gvcf.list}"
NEW_301L_GVCF="${1:-/data2/chenh/bsa/CleanData/301_check/gvcf/F9-301L.g.vcf.gz}"
if [[ $# -ge 2 ]]; then
    BASE_LIST="$2"
fi

REFERENCE="${REFERENCE:-/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa}"
FINAL_OUTDIR="${FINAL_OUTDIR:-/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d}"
FINAL_VCF="${FINAL_OUTDIR}/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz"

RUN_ID="rerun_301L_$(date '+%Y%m%d_%H%M%S')"
TEMP_ROOT="${FINAL_OUTDIR}/.${RUN_ID}"
TEMP_LIST="${TEMP_ROOT}/hunchi_and_qinben.prepared_single_sample.cleaned.gvcf.list"
TEMP_OUTDIR="${TEMP_ROOT}/joint_call"
TEMP_LOG="${TEMP_ROOT}/${RUN_ID}.log"
TEMP_FINAL_VCF="${TEMP_OUTDIR}/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz"
TEMP_FINAL_TBI="${TEMP_FINAL_VCF}.tbi"
TEMP_FINAL_CSI="${TEMP_FINAL_VCF}.csi"
TEMP_FINAL_IDX="${TEMP_FINAL_VCF}.idx"

EXPECTED_SAMPLE_301L="F9-301L"
USER_NAME="${USER:-$(id -un 2>/dev/null || echo user)}"

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

has_index() {
    local file="$1"
    [[ -f "${file}.tbi" || -f "${file}.csi" || -f "${file}.idx" ]]
}

extract_single_sample_name() {
    local gvcf="$1"
    local line

    set +o pipefail
    line="$(gzip -dc "$gvcf" | awk 'BEGIN{FS="\t"} /^#CHROM/ {print; exit}')"
    local status=$?
    set -o pipefail

    [[ $status -eq 0 ]] || return 1
    [[ -n "$line" ]] || return 1

    awk -F'\t' '
        {
            if (NF != 10) {
                exit 1
            }
            print $10
        }
    ' <<<"$line"
}

assert_single_sample_gvcf() {
    local gvcf="$1"
    local sample_name

    [[ -f "$gvcf" ]] || die "gVCF 不存在: $gvcf"
    has_index "$gvcf" || die "gVCF 缺少索引(.tbi/.csi/.idx): $gvcf"

    sample_name="$(extract_single_sample_name "$gvcf")" || die "无法读取 gVCF sample 名: $gvcf"
    [[ "$sample_name" == "$EXPECTED_SAMPLE_301L" ]] || die "新 301L gVCF 的 sample 名不是 ${EXPECTED_SAMPLE_301L}: ${sample_name}"
}

replace_301l_in_list() {
    local input_list="$1"
    local output_list="$2"
    local replacement="$3"

    awk -v replacement="$replacement" '
        BEGIN { changed = 0 }
        {
            if ($0 ~ /(^|\/)F9-301L\.g\.vcf\.gz$/) {
                print replacement
                changed++
            } else {
                print
            }
        }
        END {
            if (changed != 1) {
                exit 7
            }
        }
    ' "$input_list" > "$output_list" || die "替换 301L 路径失败；请确认基线 list 中恰好只有 1 个 F9-301L.g.vcf.gz 条目"
}

compare_lists_except_301l() {
    local base_list="$1"
    local temp_list="$2"
    local base_filtered="${TEMP_ROOT}/base_without_301L.list"
    local temp_filtered="${TEMP_ROOT}/temp_without_301L.list"

    grep -v '/F9-301L\.g\.vcf\.gz$' "$base_list" > "$base_filtered"
    grep -v '/F9-301L\.g\.vcf\.gz$' "$temp_list" > "$temp_filtered"

    diff -u "$base_filtered" "$temp_filtered" >/dev/null || die "临时 list 除 301L 外还发生了其他变化"
}

query_samples() {
    local vcf="$1"
    local out="$2"
    bcftools query -l "$vcf" > "$out"
    [[ -s "$out" ]] || die "无法从 VCF 读取样本列表: $vcf"
}

detect_final_index_path() {
    local vcf="$1"
    if [[ -f "${vcf}.tbi" ]]; then
        printf '%s\n' "${vcf}.tbi"
    elif [[ -f "${vcf}.csi" ]]; then
        printf '%s\n' "${vcf}.csi"
    elif [[ -f "${vcf}.idx" ]]; then
        printf '%s\n' "${vcf}.idx"
    else
        return 1
    fi
}

cleanup_temp_root() {
    local dir="$1"
    [[ -n "$dir" && -d "$dir" ]] || return 0
    rm -rf "$dir"
}

cleanup_main_script_scratch() {
    local list_path="$1"
    local run_prefix scratch_tag

    run_prefix="$(basename "$list_path")"
    run_prefix="${run_prefix%.gvcf.list}"
    run_prefix="${run_prefix%.list}"
    scratch_tag="${run_prefix}_scratch"

    rm -rf "/dev/shm/${USER_NAME}/${scratch_tag}" 2>/dev/null || true
    rm -rf "/scratch/${USER_NAME}/${scratch_tag}" 2>/dev/null || true
    rm -rf "/tmp/${USER_NAME}/${scratch_tag}" 2>/dev/null || true
}

mkdir -p "$TEMP_ROOT"
exec > >(tee -a "$TEMP_LOG") 2>&1

need_cmd bash
need_cmd awk
need_cmd grep
need_cmd diff
need_cmd gzip
need_cmd bcftools
need_cmd rm
need_cmd mv
need_cmd cp

[[ -f "$MAIN_SCRIPT" ]] || die "主脚本不存在: $MAIN_SCRIPT"
[[ -f "$BASE_LIST" ]] || die "基线输入 list 不存在: $BASE_LIST"
[[ -f "$REFERENCE" ]] || die "参考基因组不存在: $REFERENCE"
[[ -d "$FINAL_OUTDIR" ]] || die "正式输出目录不存在: $FINAL_OUTDIR"

log "开始用新 301L 重建全基因组全样本 VCF"
log "MAIN_SCRIPT   = $MAIN_SCRIPT"
log "BASE_LIST     = $BASE_LIST"
log "NEW_301L_GVCF = $NEW_301L_GVCF"
log "REFERENCE     = $REFERENCE"
log "FINAL_OUTDIR  = $FINAL_OUTDIR"
log "FINAL_VCF     = $FINAL_VCF"
log "TEMP_ROOT     = $TEMP_ROOT"

assert_single_sample_gvcf "$NEW_301L_GVCF"

replace_301l_in_list "$BASE_LIST" "$TEMP_LIST" "$NEW_301L_GVCF"
compare_lists_except_301l "$BASE_LIST" "$TEMP_LIST"

TEMP_LIST_LINES="$(wc -l < "$TEMP_LIST" | tr -d ' ')"
[[ "$TEMP_LIST_LINES" -eq 12 ]] || die "临时 list 应包含 12 个单样本 gVCF，实际为: $TEMP_LIST_LINES"

log "临时替换版 list 已生成: $TEMP_LIST"
log "开始运行主流程到临时目录"

bash "$MAIN_SCRIPT" "$TEMP_LIST" "$REFERENCE" "$TEMP_OUTDIR"

[[ -s "$TEMP_FINAL_VCF" ]] || die "临时最终 VCF 不存在或为空: $TEMP_FINAL_VCF"
TEMP_INDEX_PATH="$(detect_final_index_path "$TEMP_FINAL_VCF")" || die "临时最终 VCF 缺少索引: $TEMP_FINAL_VCF"

log "开始验证临时最终 VCF"

base_samples_file="${TEMP_ROOT}/base.samples.txt"
temp_samples_file="${TEMP_ROOT}/temp.samples.txt"
query_samples "$TEMP_FINAL_VCF" "$temp_samples_file"

if [[ -s "$FINAL_VCF" ]]; then
    query_samples "$FINAL_VCF" "$base_samples_file"
    diff -u "$base_samples_file" "$temp_samples_file" >/dev/null || die "临时最终 VCF 的样本列表/顺序与旧正式 VCF 不一致"
else
    # 正式 VCF 若不存在，则至少保证 12 个样本且包含 F9-301L。
    temp_sample_count="$(wc -l < "$temp_samples_file" | tr -d ' ')"
    [[ "$temp_sample_count" -eq 12 ]] || die "临时最终 VCF 的样本数不是 12: $temp_sample_count"
fi

grep -qx "$EXPECTED_SAMPLE_301L" "$temp_samples_file" || die "临时最终 VCF 中未找到样本 ${EXPECTED_SAMPLE_301L}"

bcftools view -h "$TEMP_FINAL_VCF" | grep -q '^##FORMAT=<ID=GT,' || die "临时最终 VCF 缺少 FORMAT/GT"
bcftools view -h "$TEMP_FINAL_VCF" | grep -q '^##FORMAT=<ID=DP,' || die "临时最终 VCF 缺少 FORMAT/DP"
bcftools view -h "$TEMP_FINAL_VCF" | grep -q '^##INFO=<ID=DP,' || die "临时最终 VCF 缺少 INFO/DP"

# 抽查 301L 至少可被正常查询。
bcftools view -H -s "$EXPECTED_SAMPLE_301L" "$TEMP_FINAL_VCF" | head -n 1 >/dev/null || die "无法在临时最终 VCF 中查询 ${EXPECTED_SAMPLE_301L}"
bcftools stats "$TEMP_FINAL_VCF" > "${TEMP_ROOT}/final.bcftools.stats.txt"

log "验证通过，开始覆盖正式最终 VCF"

rm -f "${FINAL_VCF}" "${FINAL_VCF}.tbi" "${FINAL_VCF}.csi" "${FINAL_VCF}.idx"
mv -f "$TEMP_FINAL_VCF" "$FINAL_VCF"

case "$TEMP_INDEX_PATH" in
    *.tbi) mv -f "$TEMP_INDEX_PATH" "${FINAL_VCF}.tbi" ;;
    *.csi) mv -f "$TEMP_INDEX_PATH" "${FINAL_VCF}.csi" ;;
    *.idx) mv -f "$TEMP_INDEX_PATH" "${FINAL_VCF}.idx" ;;
    *) die "未知的临时索引类型: $TEMP_INDEX_PATH" ;;
esac

[[ -s "$FINAL_VCF" ]] || die "覆盖后的正式最终 VCF 为空: $FINAL_VCF"
detect_final_index_path "$FINAL_VCF" >/dev/null || die "覆盖后的正式最终 VCF 缺少索引: $FINAL_VCF"

log "正式最终 VCF 已更新完成"
log "正式结果: $FINAL_VCF"

log "开始清理临时目录"
cleanup_main_script_scratch "$TEMP_LIST"
cleanup_temp_root "$TEMP_ROOT"

log "清理完成；仅保留正式最终 VCF 与索引"
