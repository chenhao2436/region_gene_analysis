#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

###############################################################################
# 脚本名称:
#   gvcf2vcf_genome_gpt_3m6d.sh
#
# 脚本目标:
#   将 hunchi_and_qinben.gvcf.list 中列出的 11 个 gVCF 文件，使用
#   GATK4 的 "GenomicsDBImport -> GenotypeGVCFs" 流程，从头联合分型，
#   最终得到 1 个全基因组 joint VCF 文件。
#
# 为什么不用旧脚本的 CombineGVCFs:
#   1. CombineGVCFs 会直接把多个 gVCF 合成更大的 gVCF，中间文件很大；
#   2. 你之前脚本按染色体并发 12 个任务，IO、内存、临时文件压力都很高；
#   3. 染色体 9 卡住时，很难判断是某条染色体的问题，还是全局资源被打爆。
#
# 为什么这里改成 GenomicsDBImport:
#   1. GenomicsDBImport 专门为多 gVCF 联合分型设计；
#   2. 每个染色体/contig 单独建一个 GenomicsDB workspace，更容易定位失败点；
#   3. 可以断点续跑，已经成功的 interval 不必重复跑；
#   4. 最后再按参考基因组顺序把各 interval 的 joint VCF 合并成 1 个全基因组 VCF。
#
# 重要设计原则:
#   1. "全基因组" 不等于 "一次性把整个基因组导入到一个超大 GenomicsDB"。
#      对大基因组、共享存储、或者 scaffold 较多的数据，这样做更容易卡住。
#      本脚本采用更稳的方案:
#         - 按参考基因组 .fai 中的 contig 顺序逐个处理
#         - 每个 contig:
#             GenomicsDBImport  ->  GenotypeGVCFs
#         - 所有 contig 完成后再 GatherVcfs 合并
#
#   2. 注意: GenomicsDBImport 要求输入的每个 gVCF 都是“单样本 gVCF”。
#      这不是 sample-name-map 的限制，而是工具本身的输入限制。
#      因此如果列表里有 multi-sample gVCF（例如双亲本合并在一个 gVCF 中），
#      必须先拆成多个单样本 gVCF，再运行本脚本。
#
#   3. 默认输出的是“联合分型后的变异位点 VCF”。
#      如果你真正想要“全位点(all-sites) VCF”，包括非变异位点，也可以做到，
#      但文件会非常大，速度和磁盘占用都会显著上升。这个开关在下面配置区中：
#         INCLUDE_NON_VARIANT_SITES="yes"
#      默认保持为 "no"。
#
# 输出目录结构:
#   /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/
#   ├── logs/                         # 总日志 + 每个 interval 的详细日志
#   ├── per_interval/                 # 每个 interval 的 workspace 和 joint VCF
#   ├── merge_rounds/                 # GatherVcfs 分轮合并中间结果
#   ├── status/                       # 每个 interval 的成功/失败标记
#   ├── tmp/                          # GATK 临时目录
#   ├── *.cleaned.gvcf.list           # 清洗后的输入列表
#   ├── *.sample_audit.tsv            # 每个 gVCF 对应的 sample 审计结果
#   ├── *.multisample_gvcf.tsv        # 若存在多样本 gVCF，会在这里列出并提前退出
#   ├── *.intervals.tsv               # 按参考基因组 .fai 生成的 interval 清单
#   ├── *.per_interval_vcf.list       # 最终待合并的分片 VCF 列表
#   └── *.genome.joint.vcf.gz         # 最终全基因组 joint VCF
#
# 依赖:
#   - bash (建议 >= 4)
#   - gatk
#   - gzip
#   - awk
#   - tee
#   - cp
#   - samtools (仅在参考基因组缺少 .fai 且允许自动构建时需要)
#
# 用法:
#   1) 直接使用脚本中的默认路径:
#        bash gvcf2vcf_genome_gpt_3m6d.sh
#
#   2) 覆盖默认输入:
#        bash gvcf2vcf_genome_gpt_3m6d.sh \
#            /path/to/hunchi_and_qinben.gvcf.list \
#            /path/to/Canz.genome.fa \
#            /path/to/output_dir
#
# 注意:
#   - 本脚本默认优先优化共享存储场景，而不是盲目开并发。
#   - 如果输出目录在 NFS/Lustre 这类共享存储上，真正的瓶颈通常是文件系统延迟，
#     而不是内存。因此下面会优先尝试把 GenomicsDB workspace 和 tmp 放到本地盘。
#   - 即使 CPU 线程很多，也不建议一开始就把 MAX_JOBS 调很高。
#     对当前 12 个单样本 gVCF 的任务，3 个并发通常比 8~12 个并发更稳。
###############################################################################

########################################
# 一、用户可调配置区
########################################

# gVCF 列表文件。
# 如果已经存在拆分好的 prepared list，则默认优先使用它；否则退回原始列表。
DEFAULT_PREPARED_LIST="/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d/hunchi_and_qinben.prepared_single_sample.gvcf.list"
DEFAULT_ORIGINAL_LIST="/data2/chenh/bsa/CleanData/vcf/hunchi_and_qinben.gvcf.list"
if [[ -f "$DEFAULT_PREPARED_LIST" ]]; then
    LIST_FILE="$DEFAULT_PREPARED_LIST"
else
    LIST_FILE="$DEFAULT_ORIGINAL_LIST"
fi

# 参考基因组 FASTA。注意不能写成目录；即使误写了末尾斜杠，下面也会自动去掉。
REFERENCE="/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa"

# 所有结果统一输出到这个目录。
OUTDIR="/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d"

# 是否输出 all-sites VCF：
#   no  = 默认，只输出联合分型后的变异位点（通常就是大家说的 joint VCF）
#   yes = 额外保留非变异位点，文件会非常大；只有你明确需要时才打开
INCLUDE_NON_VARIANT_SITES="no"

# 每次同时处理几个 interval。
# 这里按“约 100GB 内存 + 共享存储”的场景调到 3 个并发。
MAX_JOBS=3

# GenomicsDBImport 的 Java 堆内存。
# 对 11 个 gVCF 来说，这个值不需要夸张；但考虑到某些 contig 可能很大，
# 当前默认给到 20g，兼顾速度和总并发内存占用。
JAVA_XMX_IMPORT="20g"

# GenotypeGVCFs 的 Java 堆内存。
# 联合分型通常比导入更吃内存，可略高一些。
JAVA_XMX_GENOTYPE="28g"

# GatherVcfs 最终合并时的 Java 堆内存。
JAVA_XMX_GATHER="12g"

# GenomicsDBImport 读取线程数。
# 共享文件系统上线程开太大反而容易加重抖动，因此保守开到 4。
READER_THREADS=4

# GenomicsDBImport 的 batch size。
# "auto" 表示自动使用 gVCF 文件数（你这里默认是 11），即单批导入。
# 由于你这里样本量并不大，单批通常更简单，也能避免生成很多 fragment。
GVCF_BATCH_SIZE="auto"

# 最终合并时，每轮最多把多少个分片 VCF 合并成 1 个更大的 VCF。
# 这样做是为了避免 contig/scaffold 很多时，单条命令参数过长。
GATHER_CHUNK_SIZE=100

# 如果参考基因组缺少 .fai 或 .dict，是否自动创建。
#   yes = 自动补齐
#   no  = 缺什么就报错退出，由你手工准备
AUTO_BUILD_REFERENCE_INDEX="yes"

# 后台任务达到 MAX_JOBS 后，主循环每隔多少秒检查一次空闲槽位。
SLOT_CHECK_SLEEP=5

# 是否启用 GATK 官方针对共享 POSIX 文件系统(NFS/Lustre)的优化开关。
# Broad 官方文档说明：在共享文件系统变慢时，建议设为 true。
USE_SHARED_POSIXFS_OPTIMIZATIONS="yes"

# 是否优先把 GenomicsDB workspace 和临时文件放到本地 scratch，而不是 /data2。
# 对你当前这种 /data2 已用 93% 的场景，强烈建议保持 yes。
PREFER_LOCAL_SCRATCH="yes"

# 本地 scratch 根目录。
# auto: 按顺序尝试 /dev/shm、/scratch/$USER、/tmp/$USER
# 也可以手工改成明确路径，例如 /home/chenh/tmp 或 /ssd/chenh/gatk_tmp
LOCAL_SCRATCH_ROOT="auto"

# 如果 prepared gVCF 列表里记着旧路径，而文件后来被你移动了，
# 脚本会按 basename 到下面这些目录里自动寻找并重定向。
RELOCATED_GVCF_SEARCH_DIRS=(
    "/data2/chenh/bsa/CleanData/bsa_HunchiData/parent_vcf"
    "/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d/split_gvcf"
)

########################################
# 二、命令行参数覆盖默认值（可选）
########################################

usage() {
    cat <<'EOF'
用法:
  bash gvcf2vcf_genome_gpt_3m6d.sh [gvcf_list] [reference.fa] [outdir]

示例:
  bash gvcf2vcf_genome_gpt_3m6d.sh

  bash gvcf2vcf_genome_gpt_3m6d.sh \
      /data2/chenh/bsa/CleanData/vcf/hunchi_and_qinben.gvcf.list \
      /data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa \
      /data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d
EOF
}

if [[ $# -gt 3 ]]; then
    usage
    exit 1
fi

[[ $# -ge 1 ]] && LIST_FILE="$1"
[[ $# -ge 2 ]] && REFERENCE="$2"
[[ $# -ge 3 ]] && OUTDIR="$3"

########################################
# 三、路径规范化与派生变量
########################################

# 防止把 reference 写成 "...fa/" 这种形式，自动去掉末尾斜杠。
REFERENCE="${REFERENCE%/}"
OUTDIR="${OUTDIR%/}"

RUN_PREFIX="$(basename "$LIST_FILE")"
RUN_PREFIX="${RUN_PREFIX%.gvcf.list}"
RUN_PREFIX="${RUN_PREFIX%.list}"

REF_DIR="$(dirname "$REFERENCE")"
REF_BASE="$(basename "$REFERENCE")"
REF_STEM="${REF_BASE%.*}"
REF_DICT="${REF_DIR}/${REF_STEM}.dict"

LOG_DIR="${OUTDIR}/logs"
STATUS_DIR="${OUTDIR}/status"
PER_INTERVAL_ROOT="${OUTDIR}/per_interval"
MERGE_ROOT="${OUTDIR}/merge_rounds"
SCRATCH_ROOT=""
TMP_ROOT=""
SCRATCH_DB_ROOT=""

CLEAN_LIST="${OUTDIR}/${RUN_PREFIX}.cleaned.gvcf.list"
RESOLVED_LIST="${OUTDIR}/${RUN_PREFIX}.resolved_input.gvcf.list"
SAMPLE_AUDIT="${OUTDIR}/${RUN_PREFIX}.sample_audit.tsv"
MULTISAMPLE_AUDIT="${OUTDIR}/${RUN_PREFIX}.multisample_gvcf.tsv"
INTERVAL_META="${OUTDIR}/${RUN_PREFIX}.intervals.tsv"
PER_INTERVAL_VCF_LIST="${OUTDIR}/${RUN_PREFIX}.per_interval_vcf.list"
FINAL_VCF="${OUTDIR}/${RUN_PREFIX}.genome.joint.vcf.gz"

TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
MASTER_LOG="${LOG_DIR}/${RUN_PREFIX}.run_${TIMESTAMP}.log"

########################################
# 四、基础函数
########################################

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

has_variant_index() {
    local file="$1"
    [[ -f "${file}.tbi" || -f "${file}.idx" || -f "${file}.csi" ]]
}

ensure_variant_index() {
    local file="$1"
    if has_variant_index "$file"; then
        return 0
    fi

    log "为文件建立索引: $file"
    gatk IndexFeatureFile -I "$file"

    has_variant_index "$file" || die "索引创建失败: $file"
}

# 从压缩 gVCF 的表头里提取所有 sample 名称。
# 这里故意只读取到 #CHROM 行为止，避免把整个大文件解压出来。
extract_samples_from_gvcf() {
    local gvcf="$1"
    local header_line

    # 关闭 pipefail，只取 awk 的退出状态，避免 gzip 因下游提前退出而被误判失败。
    set +o pipefail
    header_line="$(gzip -dc "$gvcf" | awk 'BEGIN{FS="\t"} /^#CHROM/ {print; exit}')"
    local status=$?
    set -o pipefail

    [[ $status -eq 0 ]] || return 1
    [[ -n "$header_line" ]] || return 1

    awk -F'\t' '
        {
            if (NF <= 9) {
                exit 1
            }
            for (i = 10; i <= NF; i++) {
                print $i
            }
        }
    ' <<<"$header_line"
}

# 把 contig 名称改成适合做目录/文件名的安全字符串。
sanitize_name() {
    printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

# 输入列表中的文件如果被用户移动过，则按 basename 到候选目录里补查。
resolve_gvcf_path() {
    local raw_path="$1"
    local base_name candidate_dir candidate_path

    if [[ -f "$raw_path" ]]; then
        printf '%s\n' "$raw_path"
        return 0
    fi

    base_name="$(basename "$raw_path")"
    for candidate_dir in "${RELOCATED_GVCF_SEARCH_DIRS[@]}"; do
        candidate_path="${candidate_dir%/}/${base_name}"
        if [[ -f "$candidate_path" ]]; then
            printf '%s\n' "$candidate_path"
            return 0
        fi
    done

    return 1
}

# 选择用于 GenomicsDB workspace 和 tmp 的 scratch 目录。
# 原则:
#   1. 优先本地盘，减少 /data2 共享存储压力；
#   2. 若找不到可写本地目录，则回退到 OUTDIR/scratch；
#   3. 只把“高频读写的中间文件”放 scratch，最终 VCF 和日志仍写回 OUTDIR。
choose_scratch_root() {
    local run_tag="$1"
    local candidate=""

    if [[ "$PREFER_LOCAL_SCRATCH" == "yes" ]]; then
        if [[ "$LOCAL_SCRATCH_ROOT" != "auto" ]]; then
            candidate="${LOCAL_SCRATCH_ROOT%/}/${run_tag}"
            mkdir -p "$candidate" 2>/dev/null || true
            if [[ -d "$candidate" && -w "$candidate" ]]; then
                printf '%s\n' "$candidate"
                return 0
            fi
        else
            for candidate in \
                "/dev/shm/${USER:-$(id -un 2>/dev/null || echo user)}/${run_tag}" \
                "/scratch/${USER:-$(id -un 2>/dev/null || echo user)}/${run_tag}" \
                "/tmp/${USER:-$(id -un 2>/dev/null || echo user)}/${run_tag}"
            do
                mkdir -p "$candidate" 2>/dev/null || true
                if [[ -d "$candidate" && -w "$candidate" ]]; then
                    printf '%s\n' "$candidate"
                    return 0
                fi
            done
        fi
    fi

    candidate="${OUTDIR}/scratch/${run_tag}"
    mkdir -p "$candidate"
    printf '%s\n' "$candidate"
}

# 按 chunk 调用 GatherVcfs，避免一次性传入过多输入文件。
gather_chunk() {
    local out_vcf="$1"
    shift

    local -a gather_args=()
    local input_vcf
    for input_vcf in "$@"; do
        gather_args+=(-I "$input_vcf")
    done

    log "开始 GatherVcfs: 输出=$(basename "$out_vcf") 输入数=$#"
    gatk \
        --java-options "-Xms${JAVA_XMX_GATHER} -Xmx${JAVA_XMX_GATHER} -Djava.io.tmpdir=${TMP_ROOT}/gather" \
        GatherVcfs \
        "${gather_args[@]}" \
        -O "$out_vcf"

    ensure_variant_index "$out_vcf"
}

########################################
# 五、前置检查与准备
########################################

mkdir -p "$OUTDIR" "$LOG_DIR" "$STATUS_DIR" "$PER_INTERVAL_ROOT" "$MERGE_ROOT"

SCRATCH_ROOT="$(choose_scratch_root "${RUN_PREFIX}_scratch")"
TMP_ROOT="${SCRATCH_ROOT}/tmp"
SCRATCH_DB_ROOT="${SCRATCH_ROOT}/genomicsdb"
mkdir -p "$SCRATCH_ROOT" "$TMP_ROOT" "${TMP_ROOT}/gather" "$SCRATCH_DB_ROOT"

exec > >(tee -a "$MASTER_LOG") 2>&1

trap 'log "脚本异常终止；请优先查看总日志: '"$MASTER_LOG"'"' ERR

need_cmd bash
need_cmd gatk
need_cmd gzip
need_cmd awk
need_cmd sed
need_cmd tee
need_cmd cp

if (( BASH_VERSINFO[0] < 4 )); then
    die "本脚本需要 bash >= 4；当前版本为: $BASH_VERSION"
fi

[[ -f "$LIST_FILE" ]] || die "gVCF 列表文件不存在: $LIST_FILE"
[[ -f "$REFERENCE" ]] || die "参考基因组不存在: $REFERENCE"

log "脚本开始运行"
log "LIST_FILE        = $LIST_FILE"
log "REFERENCE        = $REFERENCE"
log "OUTDIR           = $OUTDIR"
log "SCRATCH_ROOT     = $SCRATCH_ROOT"
log "RESOLVED_LIST    = $RESOLVED_LIST"
log "INCLUDE_NON_VAR  = $INCLUDE_NON_VARIANT_SITES"
log "MAX_JOBS         = $MAX_JOBS"
log "JAVA_XMX_IMPORT  = $JAVA_XMX_IMPORT"
log "JAVA_XMX_GENOTYPE= $JAVA_XMX_GENOTYPE"
log "JAVA_XMX_GATHER  = $JAVA_XMX_GATHER"
log "READER_THREADS   = $READER_THREADS"
log "SHARED_POSIXFS   = $USE_SHARED_POSIXFS_OPTIMIZATIONS"
log "MASTER_LOG       = $MASTER_LOG"

if command -v df >/dev/null 2>&1; then
    log "输出目录文件系统信息:"
    df -h "$OUTDIR" | sed 's/^/[df OUTDIR] /'
    log "scratch 目录文件系统信息:"
    df -h "$SCRATCH_ROOT" | sed 's/^/[df SCRATCH] /'
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    if [[ "$AUTO_BUILD_REFERENCE_INDEX" != "yes" ]]; then
        die "缺少参考基因组索引: ${REFERENCE}.fai"
    fi
    need_cmd samtools
    log "未发现 ${REFERENCE}.fai，开始执行 samtools faidx"
    samtools faidx "$REFERENCE"
fi

if [[ ! -f "$REF_DICT" ]]; then
    if [[ "$AUTO_BUILD_REFERENCE_INDEX" != "yes" ]]; then
        die "缺少参考基因组字典: $REF_DICT"
    fi
    log "未发现参考字典 $REF_DICT，开始执行 GATK CreateSequenceDictionary"
    gatk CreateSequenceDictionary -R "$REFERENCE" -O "$REF_DICT"
fi

########################################
# 六、整理 gVCF 输入，并做 sample 审计
########################################

# 去掉空行和注释行，得到真正要处理的 gVCF 列表。
awk 'NF && $1 !~ /^#/' "$LIST_FILE" > "$CLEAN_LIST"
TOTAL_GVCF="$(wc -l < "$CLEAN_LIST" | tr -d ' ')"
[[ "$TOTAL_GVCF" -gt 0 ]] || die "清理后没有可用的 gVCF 输入: $CLEAN_LIST"

EXPECTED_GVCF_COUNT=11
if [[ "$RUN_PREFIX" == *prepared_single_sample* ]]; then
    EXPECTED_GVCF_COUNT=12
fi

if [[ "${TOTAL_GVCF}" -ne "${EXPECTED_GVCF_COUNT}" ]]; then
    log "警告: 清理后的 gVCF 文件数不是 ${EXPECTED_GVCF_COUNT}, 而是 ${TOTAL_GVCF}. 脚本仍按实际数量运行."
fi

if [[ "$GVCF_BATCH_SIZE" == "auto" ]]; then
    IMPORT_BATCH_SIZE="$TOTAL_GVCF"
else
    IMPORT_BATCH_SIZE="$GVCF_BATCH_SIZE"
fi

log "清洗后的 gVCF 列表: $CLEAN_LIST"
log "实际 gVCF 文件数: ${TOTAL_GVCF}"
log "GenomicsDBImport batch size: $IMPORT_BATCH_SIZE"

: > "$RESOLVED_LIST"
: > "$SAMPLE_AUDIT"
printf "gvcf_path\tsample_count\tsample_names\n" > "$SAMPLE_AUDIT"

# 这里不生成 sample-name-map；直接检查每个 gVCF 是否满足“单样本 gVCF”要求。
GVCF_ARGS=()
declare -A SEEN_SAMPLES=()
TOTAL_SAMPLES=0
INPUT_INDEX=0
MULTISAMPLE_COUNT=0

while IFS= read -r gvcf; do
    ((INPUT_INDEX += 1))
    resolved_gvcf="$(resolve_gvcf_path "$gvcf")" || die "gVCF 文件不存在，且无法在候选目录中重定位: $gvcf"
    if [[ "$resolved_gvcf" != "$gvcf" ]]; then
        log "输入路径已重定位:"
        log "  原路径: $gvcf"
        log "  新路径: $resolved_gvcf"
    fi
    gvcf="$resolved_gvcf"

    if [[ ! -f "${gvcf}.tbi" && ! -f "${gvcf}.csi" && ! -f "${gvcf}.idx" ]]; then
        die "gVCF 缺少索引文件(.tbi/.csi/.idx): $gvcf"
    fi

    printf "%s\n" "$gvcf" >> "$RESOLVED_LIST"
    GVCF_ARGS+=(-V "$gvcf")

    current_samples=()
    while IFS= read -r sample_name; do
        [[ -n "$sample_name" ]] || continue

        if [[ -n "${SEEN_SAMPLES[$sample_name]:-}" ]]; then
            die "发现重复 sample 名: ${sample_name}
首次出现: ${SEEN_SAMPLES[$sample_name]}
再次出现: $gvcf"
        fi

        SEEN_SAMPLES["$sample_name"]="$gvcf"
        current_samples+=("$sample_name")
        ((TOTAL_SAMPLES += 1))
    done < <(extract_samples_from_gvcf "$gvcf")

    [[ ${#current_samples[@]} -gt 0 ]] || die "无法从 gVCF 表头提取 sample 名称: $gvcf"

    sample_csv="$(IFS=,; echo "${current_samples[*]}")"
    printf "%s\t%d\t%s\n" "$gvcf" "${#current_samples[@]}" "$sample_csv" >> "$SAMPLE_AUDIT"
    log "输入[${INPUT_INDEX}/${TOTAL_GVCF}] ${gvcf}"
    log "  样本数: ${#current_samples[@]}"
    log "  样本名: $sample_csv"

    if [[ ${#current_samples[@]} -gt 1 ]]; then
        if [[ "$MULTISAMPLE_COUNT" -eq 0 ]]; then
            printf "gvcf_path\tsample_count\tsample_names\n" > "$MULTISAMPLE_AUDIT"
        fi
        printf "%s\t%d\t%s\n" "$gvcf" "${#current_samples[@]}" "$sample_csv" >> "$MULTISAMPLE_AUDIT"
        ((MULTISAMPLE_COUNT += 1))
    fi
done < "$CLEAN_LIST"

log "sample 审计完成: $SAMPLE_AUDIT"
log "总样本数(跨全部 gVCF 去重后): $TOTAL_SAMPLES"

if [[ "$MULTISAMPLE_COUNT" -gt 0 ]]; then
    die "检测到 ${MULTISAMPLE_COUNT} 个 multi-sample gVCF。
GenomicsDBImport 在当前流程下要求每个输入文件都是单样本 gVCF。
请先拆分这些文件，再用拆分后的新列表重跑。
多样本清单见: $MULTISAMPLE_AUDIT"
fi

########################################
# 七、根据参考基因组 .fai 生成 interval 清单
########################################

# 格式:
#   序号(补零)   contig名   contig长度   安全文件名
awk '
    BEGIN{OFS="\t"}
    {
        safe = $1
        gsub(/[^A-Za-z0-9._-]/, "_", safe)
        printf "%05d\t%s\t%s\t%s\n", NR, $1, $2, safe
    }
' "${REFERENCE}.fai" > "$INTERVAL_META"

TOTAL_INTERVALS="$(wc -l < "$INTERVAL_META" | tr -d ' ')"
[[ "$TOTAL_INTERVALS" -gt 0 ]] || die "没有从 ${REFERENCE}.fai 读到任何 contig"

log "interval 清单生成完成: $INTERVAL_META"
log "参考基因组中的 contig/scaffold 数量: $TOTAL_INTERVALS"

########################################
# 八、逐 interval 运行 GenomicsDBImport 和 GenotypeGVCFs
########################################

run_one_interval() {
    local idx="$1"
    local interval="$2"
    local interval_len="$3"
    local safe_name="$4"

    local tag="${idx}_${safe_name}"
    local interval_dir="${PER_INTERVAL_ROOT}/${tag}"
    local db_dir="${SCRATCH_DB_ROOT}/${tag}.workspace"
    local out_vcf="${interval_dir}/${tag}.joint.vcf.gz"
    local tmp_dir="${TMP_ROOT}/${tag}"
    local interval_log="${LOG_DIR}/${tag}.log"
    local ok_flag="${STATUS_DIR}/${tag}.ok"
    local fail_flag="${STATUS_DIR}/${tag}.fail"

    rm -f "$ok_flag" "$fail_flag"

    {
        log "[$tag] 开始处理 interval=${interval} length=${interval_len}"
        mkdir -p "$interval_dir" "$tmp_dir"

        # 如果这个 interval 的最终结果已经存在且有索引，则直接跳过。
        # 这样脚本重新执行时，不会重复计算已经成功的 interval。
        if [[ -s "$out_vcf" ]] && has_variant_index "$out_vcf"; then
            log "[$tag] 已检测到完整结果，跳过重跑: $out_vcf"
            touch "$ok_flag"
            return 0
        fi

        # 如果 workspace 残留而最终 VCF 不完整，最稳妥的方式是删除残留后重建。
        # 因为不完整的 GenomicsDB workspace 经常会导致下一次导入报错或表现异常。
        if [[ -d "$db_dir" ]]; then
            log "[$tag] 发现残留 workspace，删除后重建: $db_dir"
            rm -rf "$db_dir"
        fi

        rm -f "$out_vcf" "${out_vcf}.tbi" "${out_vcf}.idx" "${out_vcf}.csi"

        log "[$tag] Step 1/2: GenomicsDBImport"
        local -a import_extra_args=()
        if [[ "$USE_SHARED_POSIXFS_OPTIMIZATIONS" == "yes" ]]; then
            import_extra_args+=(--genomicsdb-shared-posixfs-optimizations true)
        fi

        if ! gatk \
            --java-options "-Xms${JAVA_XMX_IMPORT} -Xmx${JAVA_XMX_IMPORT} -Djava.io.tmpdir=${tmp_dir}" \
            GenomicsDBImport \
            -R "$REFERENCE" \
            "${GVCF_ARGS[@]}" \
            -L "$interval" \
            --genomicsdb-workspace-path "$db_dir" \
            --batch-size "$IMPORT_BATCH_SIZE" \
            --reader-threads "$READER_THREADS" \
            "${import_extra_args[@]}" \
            --tmp-dir "$tmp_dir"; then
            log "[$tag] GenomicsDBImport 失败，停止当前 interval，避免产生后续连锁假报错"
            return 1
        fi

        log "[$tag] Step 2/2: GenotypeGVCFs"
        local -a genotype_extra_args=()
        if [[ "$INCLUDE_NON_VARIANT_SITES" == "yes" ]]; then
            genotype_extra_args+=(--include-non-variant-sites)
        fi
        if [[ "$USE_SHARED_POSIXFS_OPTIMIZATIONS" == "yes" ]]; then
            genotype_extra_args+=(--genomicsdb-shared-posixfs-optimizations true)
        fi

        if ! gatk \
            --java-options "-Xms${JAVA_XMX_GENOTYPE} -Xmx${JAVA_XMX_GENOTYPE} -Djava.io.tmpdir=${tmp_dir}" \
            GenotypeGVCFs \
            -R "$REFERENCE" \
            -V "gendb://${db_dir}" \
            -L "$interval" \
            "${genotype_extra_args[@]}" \
            -O "$out_vcf" \
            --tmp-dir "$tmp_dir"; then
            log "[$tag] GenotypeGVCFs 失败，停止当前 interval"
            return 1
        fi

        if ! ensure_variant_index "$out_vcf"; then
            log "[$tag] 输出 VCF 索引失败，停止当前 interval"
            return 1
        fi

        # interval 成功后清理临时目录和 workspace，避免 scratch 越积越大。
        rm -rf "$db_dir"
        rm -rf "$tmp_dir"

        touch "$ok_flag"
        log "[$tag] 完成: $out_vcf"
    } > >(tee -a "$interval_log") 2>&1 || {
        touch "$fail_flag"
        log "[$tag] 失败，请查看日志: $interval_log"
        return 1
    }
}

log "开始逐 interval 联合分型"

if [[ "$MAX_JOBS" -le 1 ]]; then
    while IFS=$'\t' read -r idx interval interval_len safe_name; do
        run_one_interval "$idx" "$interval" "$interval_len" "$safe_name"
    done < "$INTERVAL_META"
else
    while IFS=$'\t' read -r idx interval interval_len safe_name; do
        while [[ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAX_JOBS" ]]; do
            sleep "$SLOT_CHECK_SLEEP"
        done

        run_one_interval "$idx" "$interval" "$interval_len" "$safe_name" &
    done < "$INTERVAL_META"

    # 后台任务即便有失败，也先等全部结束，再统一汇总失败 interval。
    wait || true
fi

fail_flags=( "${STATUS_DIR}"/*.fail )
if [[ ${#fail_flags[@]} -gt 0 ]]; then
    log "以下 interval 处理失败，请按对应日志排查:"
    for fail_flag in "${fail_flags[@]}"; do
        log "  $(basename "${fail_flag%.fail}")"
    done
    die "存在失败 interval，已停止进入最终合并步骤。"
fi

ok_flags=( "${STATUS_DIR}"/*.ok )
log "interval 处理完成: 成功 ${#ok_flags[@]} / ${TOTAL_INTERVALS}"
[[ ${#ok_flags[@]} -eq "$TOTAL_INTERVALS" ]] || die "成功 interval 数与参考基因组 contig 数不一致"

########################################
# 九、生成分片 VCF 列表并递归合并为全基因组 VCF
########################################

: > "$PER_INTERVAL_VCF_LIST"

while IFS=$'\t' read -r idx interval interval_len safe_name; do
    tag="${idx}_${safe_name}"
    shard_vcf="${PER_INTERVAL_ROOT}/${tag}/${tag}.joint.vcf.gz"

    [[ -s "$shard_vcf" ]] || die "缺少 interval 输出 VCF: $shard_vcf"
    ensure_variant_index "$shard_vcf"
    printf "%s\n" "$shard_vcf" >> "$PER_INTERVAL_VCF_LIST"
done < "$INTERVAL_META"

log "分片 VCF 列表生成完成: $PER_INTERVAL_VCF_LIST"

CURRENT_LIST="$PER_INTERVAL_VCF_LIST"
ROUND=1

while true; do
    mapfile -t ROUND_INPUTS < "$CURRENT_LIST"
    [[ ${#ROUND_INPUTS[@]} -gt 0 ]] || die "合并输入列表为空: $CURRENT_LIST"

    if [[ ${#ROUND_INPUTS[@]} -eq 1 ]]; then
        log "最终只剩 1 个输入，复制为最终结果"
        cp -f "${ROUND_INPUTS[0]}" "$FINAL_VCF"
        rm -f "${FINAL_VCF}.tbi" "${FINAL_VCF}.idx" "${FINAL_VCF}.csi"
        ensure_variant_index "$FINAL_VCF"
        break
    fi

    ROUND_DIR="${MERGE_ROOT}/round_$(printf '%02d' "$ROUND")"
    mkdir -p "$ROUND_DIR"
    NEXT_LIST="${ROUND_DIR}/next_round.list"
    : > "$NEXT_LIST"

    log "开始最终合并第 $ROUND 轮: 输入数=${#ROUND_INPUTS[@]}"

    chunk_inputs=()
    chunk_id=0
    for shard_vcf in "${ROUND_INPUTS[@]}"; do
        chunk_inputs+=("$shard_vcf")

        if [[ ${#chunk_inputs[@]} -eq "$GATHER_CHUNK_SIZE" ]]; then
            chunk_out="${ROUND_DIR}/chunk_$(printf '%05d' "$chunk_id").vcf.gz"
            gather_chunk "$chunk_out" "${chunk_inputs[@]}"
            printf "%s\n" "$chunk_out" >> "$NEXT_LIST"
            chunk_inputs=()
            ((chunk_id += 1))
        fi
    done

    if [[ ${#chunk_inputs[@]} -gt 0 ]]; then
        chunk_out="${ROUND_DIR}/chunk_$(printf '%05d' "$chunk_id").vcf.gz"
        gather_chunk "$chunk_out" "${chunk_inputs[@]}"
        printf "%s\n" "$chunk_out" >> "$NEXT_LIST"
    fi

    CURRENT_LIST="$NEXT_LIST"
    ((ROUND += 1))
done

########################################
# 十、收尾与结果汇总
########################################

log "全基因组 joint VCF 生成完成"
log "最终结果: $FINAL_VCF"

if has_variant_index "$FINAL_VCF"; then
    log "最终索引存在: 是"
else
    log "最终索引存在: 否"
fi

log "主日志: $MASTER_LOG"
log "脚本执行完成"
