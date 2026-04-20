#!/usr/bin/env bash
# 生成日期: 2026-04-18
# 脚本作用: 串联执行 KASP 候选排序、VCF 过滤、extract_fa 输入整理和序列提取。
# 输入文件: 20260418_kasp_pipeline.ini、snpindex 文件、VCF 文件、参考基因组、scripts/extract_fa.py。
# 输出文件: 01_ranked_candidates.tsv、01_rank_summary.tsv、02_checked_sites.tsv、02_pass_sites.tsv、03_extract_fa_input.tsv、03_recommended_markers.tsv、04_kasp_sequences.tsv、pipeline.log。
# 原理与逻辑:
# 1. 统一从配置文件读取默认值, 并允许命令行覆盖关键参数。
# 2. 启动前严格检查依赖和输入文件, 防止流程跑到中途才报错。
# 3. 固定执行 01 -> 02 -> 03 -> extract_fa.py, 保证最终一条命令跑完整条链路。
# 使用方法:
# bash 20260418_kasp_00_run_pipeline.sh 20260418_kasp_pipeline.ini
# bash 20260418_kasp_00_run_pipeline.sh --config 20260418_kasp_pipeline.ini --region Chr05:248.7-255.3MB --snpindex input.tsv --high-pool SAMPLE_H --low-pool SAMPLE_L --outdir ./result_dir

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_FILE="${BASE_DIR}/20260418_kasp_pipeline.ini"

REGION_OVERRIDE=""
SNPINDEX_OVERRIDE=""
HIGH_POOL_OVERRIDE=""
LOW_POOL_OVERRIDE=""
OUTDIR_OVERRIDE=""

usage() {
    cat <<'EOF'
用法:
  bash 20260418_kasp_00_run_pipeline.sh [config.ini]
  bash 20260418_kasp_00_run_pipeline.sh --config config.ini
  bash 20260418_kasp_00_run_pipeline.sh --config config.ini --region Chr05:248.7-255.3MB --snpindex input.tsv --high-pool SAMPLE_H --low-pool SAMPLE_L [--outdir ./result_dir]
EOF
}

require_cmd() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[错误] 缺少命令: ${cmd}" >&2
        exit 1
    }
}

load_config() {
    if [[ ! -f "${CONFIG_FILE}" ]]; then
        echo "[错误] 配置文件不存在: ${CONFIG_FILE}" >&2
        exit 1
    fi
    set -a
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
    set +a
}

sanitize_region() {
    local raw="$1"
    raw="${raw//:/_}"
    raw="${raw//-/_}"
    raw="${raw// /}"
    raw="${raw//\//_}"
    echo "${raw}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --region)
            REGION_OVERRIDE="$2"
            shift 2
            ;;
        --snpindex)
            SNPINDEX_OVERRIDE="$2"
            shift 2
            ;;
        --high-pool)
            HIGH_POOL_OVERRIDE="$2"
            shift 2
            ;;
        --low-pool)
            LOW_POOL_OVERRIDE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR_OVERRIDE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            if [[ -z "${POSITIONAL_CONFIG:-}" ]]; then
                POSITIONAL_CONFIG="$1"
                shift
            else
                echo "[错误] 未知参数: $1" >&2
                usage
                exit 1
            fi
            ;;
    esac
done

if [[ -n "${POSITIONAL_CONFIG:-}" ]]; then
    CONFIG_FILE="${POSITIONAL_CONFIG}"
fi

load_config

region="${REGION_OVERRIDE:-${region:-}}"
snpindex_path="${SNPINDEX_OVERRIDE:-${snpindex_path:-}}"
vcf_path="${vcf_path:-}"
reference_fasta="${reference_fasta:-}"
parent1_sample="${parent1_sample:-}"
parent2_sample="${parent2_sample:-}"
high_pool_sample="${HIGH_POOL_OVERRIDE:-${high_pool_sample:-}}"
low_pool_sample="${LOW_POOL_OVERRIDE:-${low_pool_sample:-}}"
flank_no_variant_bp="${flank_no_variant_bp:-50}"
target_markers_per_side="${target_markers_per_side:-2}"
initial_rank_limit="${initial_rank_limit:-20}"
expansion_step_bp="${expansion_step_bp:-1000000}"
extract_fa_script="${SCRIPT_DIR}/extract_fa.py"

if [[ -z "${region}" ]]; then
    echo "[错误] region 未设置, 请在配置文件中填写或通过 --region 传入" >&2
    exit 1
fi
if [[ -z "${snpindex_path}" ]]; then
    echo "[错误] snpindex_path 未设置, 请在配置文件中填写或通过 --snpindex 传入" >&2
    exit 1
fi
if [[ -z "${high_pool_sample}" || -z "${low_pool_sample}" ]]; then
    echo "[错误] high_pool_sample 和 low_pool_sample 必须提供" >&2
    exit 1
fi

if [[ -n "${OUTDIR_OVERRIDE}" ]]; then
    outdir="${OUTDIR_OVERRIDE}"
elif [[ -n "${outdir:-}" ]]; then
    outdir="${outdir}"
else
    outdir="${BASE_DIR}/results/$(sanitize_region "${region}")"
fi

require_cmd python
require_cmd bcftools
require_cmd tabix

[[ -f "${snpindex_path}" ]] || { echo "[错误] snpindex 文件不存在: ${snpindex_path}" >&2; exit 1; }
[[ -f "${vcf_path}" ]] || { echo "[错误] VCF 文件不存在: ${vcf_path}" >&2; exit 1; }
[[ -f "${vcf_path}.tbi" || -f "${vcf_path}.csi" ]] || { echo "[错误] VCF 缺少索引文件: ${vcf_path}.tbi/.csi" >&2; exit 1; }
[[ -f "${reference_fasta}" ]] || { echo "[错误] 参考基因组不存在: ${reference_fasta}" >&2; exit 1; }
[[ -f "${extract_fa_script}" ]] || { echo "[错误] extract_fa.py 不存在: ${extract_fa_script}" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/20260418_kasp_01_rank_snpindex_candidates.py" ]] || { echo "[错误] 缺少步骤 1 脚本" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/20260418_kasp_02_filter_isolated_parental_snps.py" ]] || { echo "[错误] 缺少步骤 2 脚本" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/20260418_kasp_03_prepare_extract_fa_input.py" ]] || { echo "[错误] 缺少步骤 3 脚本" >&2; exit 1; }

mkdir -p "${outdir}"
log_file="${outdir}/pipeline.log"
exec > >(tee -a "${log_file}") 2>&1

echo "[开始] KASP 标记筛选流程"
echo "[参数] region=${region}"
echo "[参数] snpindex_path=${snpindex_path}"
echo "[参数] vcf_path=${vcf_path}"
echo "[参数] reference_fasta=${reference_fasta}"
echo "[参数] parent1_sample=${parent1_sample}"
echo "[参数] parent2_sample=${parent2_sample}"
echo "[参数] high_pool_sample=${high_pool_sample}"
echo "[参数] low_pool_sample=${low_pool_sample}"
echo "[参数] outdir=${outdir}"
echo "[参数] expansion_step_bp=${expansion_step_bp}"

ranked_candidates="${outdir}/01_ranked_candidates.tsv"
rank_summary="${outdir}/01_rank_summary.tsv"
checked_sites="${outdir}/02_checked_sites.tsv"
pass_sites="${outdir}/02_pass_sites.tsv"
extract_input="${outdir}/03_extract_fa_input.tsv"
recommended_markers="${outdir}/03_recommended_markers.tsv"
kasp_sequences="${outdir}/04_kasp_sequences.tsv"

echo "[步骤1] 提取并排序候选位点"
python "${SCRIPT_DIR}/20260418_kasp_01_rank_snpindex_candidates.py" \
    --region "${region}" \
    --snpindex "${snpindex_path}" \
    --output-candidates "${ranked_candidates}" \
    --output-summary "${rank_summary}" \
    --initial-rank-limit "${initial_rank_limit}" \
    --expansion-step-bp "${expansion_step_bp}"

echo "[步骤2] 过滤双亲纯合差异且无邻近变异的位点"
python "${SCRIPT_DIR}/20260418_kasp_02_filter_isolated_parental_snps.py" \
    --candidates "${ranked_candidates}" \
    --vcf "${vcf_path}" \
    --parent1 "${parent1_sample}" \
    --parent2 "${parent2_sample}" \
    --high-pool "${high_pool_sample}" \
    --low-pool "${low_pool_sample}" \
    --flank-bp "${flank_no_variant_bp}" \
    --target-markers-per-side "${target_markers_per_side}" \
    --initial-rank-limit "${initial_rank_limit}" \
    --output-checked "${checked_sites}" \
    --output-pass "${pass_sites}"

echo "[步骤3] 整理 extract_fa 输入和推荐清单"
python "${SCRIPT_DIR}/20260418_kasp_03_prepare_extract_fa_input.py" \
    --pass-sites "${pass_sites}" \
    --checked-sites "${checked_sites}" \
    --output-extract-input "${extract_input}" \
    --output-recommended "${recommended_markers}" \
    --target-markers-per-side "${target_markers_per_side}"

echo "[步骤4] 调用 extract_fa.py 提取序列"
python "${extract_fa_script}" \
    --input "${extract_input}" \
    --output "${kasp_sequences}" \
    --ref "${reference_fasta}"

echo "[完成] 全流程结束"
echo "[结果] 候选排序: ${ranked_candidates}"
echo "[结果] 排序汇总: ${rank_summary}"
echo "[结果] 检查明细: ${checked_sites}"
echo "[结果] 通过位点: ${pass_sites}"
echo "[结果] extract_fa 输入: ${extract_input}"
echo "[结果] 推荐清单: ${recommended_markers}"
echo "[结果] KASP 序列: ${kasp_sequences}"
echo "[结果] 日志: ${log_file}"
