#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

BASE_VCF="/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/VCF_FILE/BSA_FINAL_RESULT_COMPLETE.vcf.gz"
PREPARED_GVCF_LIST="/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_prepare_single_sample_gpt_3m8d/hunchi_and_qinben.prepared_single_sample.gvcf.list"
REFERENCE="/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.fa"
REFERENCE_FAI="${REFERENCE}.fai"
REFERENCE_DICT="/data2/chenh/cal_tree/1000+vcf/compare_tree/Canz_ref/Canz.genome.dict"

OUTDIR="/data2/chenh/bsa/CleanData/vcf/fixed_Chr09Head65MB_genomo_vcf"
FINAL_VCF="${OUTDIR}/BSA_genomo_3.12.vcf.gz"

THREADS=16
READER_THREADS=16
IMPORT_BATCH_SIZE=12
JAVA_XMX_IMPORT="64g"
JAVA_XMX_GENOTYPE="96g"
BIN_SIZE=5000000

CHR09_REGION="Chr09"
PRE09_REGION="Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08"
POST09_REGION="Chr10,Chr11,Chr12"
GENOME_REGION="Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12"
STRICT_EXPR='QUAL>20 && INFO/DP>10 && F_MISSING<0.5 && INFO/MQ>30'
RELAXED_EXPR='QUAL>10 && INFO/DP>5 && F_MISSING<0.9'
REQUIRED_PARENTS=("C3_Parent" "C6_Parent")

TMP_DIR="${OUTDIR}/tmp"
QC_DIR="${OUTDIR}/qc"
LOG_DIR="${OUTDIR}/logs"
RUN_LOG="${QC_DIR}/run.log"
CHR09_QC="${QC_DIR}/head_qc.tsv"
GENOME_QC="${QC_DIR}/genome_qc.tsv"
WARNINGS_FILE="${QC_DIR}/warnings.txt"
FINAL_SUMMARY="${QC_DIR}/final_summary.txt"
JOINTCALL_SAMPLE_AUDIT="${QC_DIR}/jointcall_sample_audit.tsv"
CHR09_STATS="${QC_DIR}/chr09_full.bcftools.stats.txt"
GENOME_STATS="${QC_DIR}/genome_3.12.bcftools.stats.txt"

TMP_CLEAN_LIST="${TMP_DIR}/prepared_single_sample.cleaned.gvcf.list"
TMP_BASE_SAMPLES="${TMP_DIR}/base.samples.txt"
TMP_BASE_SAMPLES_SORTED="${TMP_DIR}/base.samples.sorted.txt"
TMP_CHR09_SAMPLES="${TMP_DIR}/chr09.samples.txt"
TMP_CHR09_SAMPLES_SORTED="${TMP_DIR}/chr09.samples.sorted.txt"
TMP_CHR09_RAW_VCF="${TMP_DIR}/Chr09.12samples.joint.raw.vcf.gz"
TMP_CHR09_MERGED="${TMP_DIR}/Chr09.12samples.joint.vcf.gz"
TMP_CHR09_ORDERED="${TMP_DIR}/Chr09.12samples.ordered.vcf.gz"
TMP_PART1="${TMP_DIR}/Part1_Before09.vcf.gz"
TMP_PART3="${TMP_DIR}/Part3_After09.vcf.gz"

CHR09_SNPINDEX_CHECK="FAIL"
GENOME_SNPINDEX_CHECK="FAIL"

SCRATCH_ROOT=""
JOINT_DB_DIR=""
JOINT_TMP_DIR=""

declare -a GVCF_ARGS=()

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

warn() {
    log "WARNING: $*"
    printf '%s\n' "$*" >> "$WARNINGS_FILE"
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"
}

assert_file() {
    local file="$1"
    [[ -s "$file" ]] || die "Missing or empty file: $file"
}

has_index() {
    local file="$1"
    [[ -f "${file}.tbi" || -f "${file}.csi" || -f "${file}.idx" ]]
}

ensure_index() {
    local file="$1"
    if has_index "$file"; then
        return 0
    fi
    log "Indexing $file"
    tabix -f -p vcf "$file"
    has_index "$file" || die "Failed to index $file"
}

choose_scratch_root() {
    local run_tag="$1"
    local user_name
    user_name="${USER:-$(id -un 2>/dev/null || echo user)}"

    if [[ -d "/dev/shm" && -w "/dev/shm" ]]; then
        printf '%s\n' "/dev/shm/${user_name}/${run_tag}"
        return 0
    fi
    if [[ -d "/scratch/${user_name}" && -w "/scratch/${user_name}" ]]; then
        printf '%s\n' "/scratch/${user_name}/${run_tag}"
        return 0
    fi
    if [[ -d "/tmp" && -w "/tmp" ]]; then
        printf '%s\n' "/tmp/${user_name}/${run_tag}"
        return 0
    fi
    printf '%s\n' "${OUTDIR}/scratch/${run_tag}"
}

extract_samples() {
    local vcf="$1"
    local out="$2"
    bcftools query -l "$vcf" > "$out"
    [[ -s "$out" ]] || die "No samples found in $vcf"
}

compare_sample_sets() {
    local left="$1"
    local right="$2"
    local left_sorted="$3"
    local right_sorted="$4"
    sort "$left" > "$left_sorted"
    sort "$right" > "$right_sorted"
    diff -u "$left_sorted" "$right_sorted" >/dev/null
}

prepare_input_list() {
    awk 'NF && $1 !~ /^#/' "$PREPARED_GVCF_LIST" > "$TMP_CLEAN_LIST"
    [[ -s "$TMP_CLEAN_LIST" ]] || die "Prepared gVCF list is empty after cleanup: $TMP_CLEAN_LIST"
}

audit_prepared_gvcfs() {
    local gvcf sample_name sample_file sample_count total_samples=0
    declare -A seen_samples=()

    printf "gvcf_path\tsample_count\tsample_name\n" > "$JOINTCALL_SAMPLE_AUDIT"
    GVCF_ARGS=()

    while IFS= read -r gvcf; do
        [[ -n "$gvcf" ]] || continue
        assert_file "$gvcf"
        has_index "$gvcf" || die "gVCF is missing index: $gvcf"

        sample_file="${TMP_DIR}/$(basename "$gvcf").samples.txt"
        extract_samples "$gvcf" "$sample_file"
        sample_count="$(wc -l < "$sample_file" | tr -d ' ')"
        [[ "$sample_count" -eq 1 ]] || die "Prepared list contains multi-sample gVCF: $gvcf"

        sample_name="$(head -n 1 "$sample_file")"
        [[ -n "$sample_name" ]] || die "Failed to extract sample name from $gvcf"
        [[ -z "${seen_samples[$sample_name]:-}" ]] || die "Duplicate sample detected in prepared list: $sample_name"

        seen_samples["$sample_name"]="$gvcf"
        GVCF_ARGS+=(-V "$gvcf")
        ((total_samples += 1))
        printf "%s\t1\t%s\n" "$gvcf" "$sample_name" >> "$JOINTCALL_SAMPLE_AUDIT"
    done < "$TMP_CLEAN_LIST"

    [[ "$total_samples" -eq 12 ]] || die "Expected 12 unique samples in prepared list, found ${total_samples}"

    for parent_sample in "${REQUIRED_PARENTS[@]}"; do
        [[ -n "${seen_samples[$parent_sample]:-}" ]] || die "Required parent sample missing from prepared list: $parent_sample"
    done
}

extract_header_defs() {
    local vcf="$1"
    local out="$2"
    bcftools view -h "$vcf" | awk '
        /^##contig=<ID=/ || /^##INFO=<ID=/ || /^##FORMAT=<ID=/ {
            line = $0
            key = line
            sub(/^##/, "", key)
            sub(/=<ID=/, "\t", key)
            sub(/,.*/, "", key)
            print key "\t" line
        }
    ' > "$out"
}

compare_header_compatibility() {
    local reference_vcf="$1"
    local query_vcf="$2"
    local label="$3"
    local ref_defs="${TMP_DIR}/${label}.ref_defs.tsv"
    local qry_defs="${TMP_DIR}/${label}.qry_defs.tsv"

    extract_header_defs "$reference_vcf" "$ref_defs"
    extract_header_defs "$query_vcf" "$qry_defs"

    awk -F'\t' '
        NR == FNR {
            ref[$1 FS $2] = $3
            next
        }
        {
            key = $1 FS $2
            if (!(key in ref)) {
                printf "Missing definition in base header for %s: %s\n", key, $3 > "/dev/stderr"
                bad = 1
            } else if (ref[key] != $3) {
                printf "Conflict for %s: %s <> %s\n", key, ref[key], $3 > "/dev/stderr"
                bad = 1
            }
        }
        END {
            exit bad
        }
    ' "$ref_defs" "$qry_defs" || die "Header compatibility check failed for $query_vcf against $reference_vcf"
}

has_header_id() {
    local vcf="$1"
    local header_type="$2"
    local header_id="$3"
    bcftools view -h "$vcf" | grep -q "^##${header_type}=<ID=${header_id},"
}

first_last_pos() {
    local vcf="$1"
    local region="$2"
    bcftools query -r "$region" -f '%POS\n' "$vcf" | awk '
        NR == 1 { first = $1 }
        { last = $1 }
        END {
            if (NR == 0) print "NA\tNA"
            else print first "\t" last
        }
    '
}

count_records() {
    local vcf="$1"
    local region="${2:-}"
    if [[ -n "$region" ]]; then
        bcftools view -H -r "$region" "$vcf" | wc -l | tr -d ' '
    else
        bcftools view -H "$vcf" | wc -l | tr -d ' '
    fi
}

count_biallelic_snps() {
    local vcf="$1"
    local region="${2:-}"
    local expr="${3:-}"
    local -a cmd=(bcftools view -H -v snps -m2 -M2)
    if [[ -n "$region" ]]; then
        cmd+=(-r "$region")
    fi
    if [[ -n "$expr" ]]; then
        cmd+=(-i "$expr")
    fi
    cmd+=("$vcf")
    "${cmd[@]}" | wc -l | tr -d ' '
}

count_bins_with_biallelic_snps() {
    local vcf="$1"
    local region="$2"
    local expr="${3:-}"
    local -a cmd=(bcftools view -H -r "$region" -v snps -m2 -M2)
    if [[ -n "$expr" ]]; then
        cmd+=(-i "$expr")
    fi
    cmd+=("$vcf")
    "${cmd[@]}" | awk -v bin_size="$BIN_SIZE" '
        {
            bin = int(($2 - 1) / bin_size) + 1
            seen[bin] = 1
        }
        END {
            for (bin in seen) count++
            print count + 0
        }
    '
}

append_qc_row() {
    local outfile="$1"
    local scope="$2"
    local tier="$3"
    local records="$4"
    local snps="$5"
    local bins="$6"
    local first_pos="$7"
    local last_pos="$8"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$scope" "$tier" "$records" "$snps" "$bins" "$first_pos" "$last_pos" "$(date '+%F %T')" >> "$outfile"
}

compute_qc_metrics() {
    local vcf="$1"
    local region="$2"
    local label="$3"
    local outfile="$4"
    local expr="${5:-}"
    local records snps bins pos_range first_pos last_pos

    records="$(count_records "$vcf" "$region")"
    snps="$(count_biallelic_snps "$vcf" "$region" "$expr")"
    bins="$(count_bins_with_biallelic_snps "$vcf" "$region" "$expr")"
    pos_range="$(first_last_pos "$vcf" "$region")"
    first_pos="${pos_range%%$'\t'*}"
    last_pos="${pos_range##*$'\t'}"

    append_qc_row "$outfile" "$region" "$label" "$records" "$snps" "$bins" "$first_pos" "$last_pos"
}

run_joint_call_for_chr09() {
    rm -rf "$JOINT_DB_DIR" "$JOINT_TMP_DIR"
    mkdir -p "$JOINT_TMP_DIR"
    [[ ! -e "$JOINT_DB_DIR" ]] || die "Failed to remove stale GenomicsDB workspace: $JOINT_DB_DIR"
    rm -f "$TMP_CHR09_RAW_VCF" "${TMP_CHR09_RAW_VCF}.tbi" "${TMP_CHR09_RAW_VCF}.csi" "${TMP_CHR09_RAW_VCF}.idx"

    log "Running GenomicsDBImport for ${CHR09_REGION}"
    gatk \
        --java-options "-Xms${JAVA_XMX_IMPORT} -Xmx${JAVA_XMX_IMPORT} -Djava.io.tmpdir=${JOINT_TMP_DIR}" \
        GenomicsDBImport \
        -R "$REFERENCE" \
        "${GVCF_ARGS[@]}" \
        -L "$CHR09_REGION" \
        --genomicsdb-workspace-path "$JOINT_DB_DIR" \
        --batch-size "$IMPORT_BATCH_SIZE" \
        --reader-threads "$READER_THREADS" \
        --tmp-dir "$JOINT_TMP_DIR"

    log "Running GenotypeGVCFs for ${CHR09_REGION}"
    gatk \
        --java-options "-Xms${JAVA_XMX_GENOTYPE} -Xmx${JAVA_XMX_GENOTYPE} -Djava.io.tmpdir=${JOINT_TMP_DIR}" \
        GenotypeGVCFs \
        -R "$REFERENCE" \
        -V "gendb://${JOINT_DB_DIR}" \
        -L "$CHR09_REGION" \
        -O "$TMP_CHR09_RAW_VCF" \
        --tmp-dir "$JOINT_TMP_DIR"

    ensure_index "$TMP_CHR09_RAW_VCF"
}

normalize_chr09_to_base_order() {
    extract_samples "$BASE_VCF" "$TMP_BASE_SAMPLES"
    sort "$TMP_BASE_SAMPLES" > "$TMP_BASE_SAMPLES_SORTED"
    extract_samples "$TMP_CHR09_RAW_VCF" "$TMP_CHR09_SAMPLES"

    if ! compare_sample_sets "$TMP_BASE_SAMPLES" "$TMP_CHR09_SAMPLES" "$TMP_BASE_SAMPLES_SORTED" "$TMP_CHR09_SAMPLES_SORTED"; then
        die "Sample set mismatch between new 12-sample Chr09 VCF and $BASE_VCF"
    fi

    compare_header_compatibility "$BASE_VCF" "$TMP_CHR09_RAW_VCF" "base_vs_new_chr09"

    if diff -u "$TMP_BASE_SAMPLES" "$TMP_CHR09_SAMPLES" >/dev/null; then
        cp -f "$TMP_CHR09_RAW_VCF" "$TMP_CHR09_MERGED"
        rm -f "${TMP_CHR09_MERGED}.tbi" "${TMP_CHR09_MERGED}.csi" "${TMP_CHR09_MERGED}.idx"
        ensure_index "$TMP_CHR09_MERGED"
        return 0
    fi

    log "Reordering new 12-sample Chr09 VCF to match base sample order"
    bcftools view -S "$TMP_BASE_SAMPLES" -Oz --threads "$THREADS" -o "$TMP_CHR09_ORDERED" "$TMP_CHR09_RAW_VCF"
    ensure_index "$TMP_CHR09_ORDERED"
    cp -f "$TMP_CHR09_ORDERED" "$TMP_CHR09_MERGED"
    rm -f "${TMP_CHR09_MERGED}.tbi" "${TMP_CHR09_MERGED}.csi" "${TMP_CHR09_MERGED}.idx"
    ensure_index "$TMP_CHR09_MERGED"
}

run_bcftools_stats() {
    local vcf="$1"
    local out="$2"
    bcftools stats "$vcf" > "$out"
}

write_summary() {
    local chr09_gt="$1"
    local chr09_fmt_dp="$2"
    local chr09_ad="$3"
    local chr09_info_dp="$4"
    local chr09_info_mq="$5"
    local genome_gt="$6"
    local genome_fmt_dp="$7"
    local genome_ad="$8"
    local genome_info_dp="$9"
    local genome_info_mq="${10}"
    local base_pre09="${11}"
    local final_pre09="${12}"
    local base_post09="${13}"
    local final_post09="${14}"
    local final_chr09_raw="${15}"
    local final_chr09_relaxed="${16}"
    local final_chr09_strict="${17}"

    cat > "$FINAL_SUMMARY" <<EOF
CHR09_SNPINDEX_CHECK=${CHR09_SNPINDEX_CHECK}
GENOME_SNPINDEX_CHECK=${GENOME_SNPINDEX_CHECK}
FINAL_VCF=${FINAL_VCF}
PREPARED_GVCF_LIST=${PREPARED_GVCF_LIST}
REFERENCE=${REFERENCE}

[chr09_fields]
GT=${chr09_gt}
FORMAT_DP=${chr09_fmt_dp}
FORMAT_AD=${chr09_ad}
INFO_DP=${chr09_info_dp}
INFO_MQ=${chr09_info_mq}

[genome_fields]
GT=${genome_gt}
FORMAT_DP=${genome_fmt_dp}
FORMAT_AD=${genome_ad}
INFO_DP=${genome_info_dp}
INFO_MQ=${genome_info_mq}

[non_chr09_region_counts]
BASE_PRE09=${base_pre09}
FINAL_PRE09=${final_pre09}
BASE_POST09=${base_post09}
FINAL_POST09=${final_post09}

[chr09_counts]
FINAL_RAW=${final_chr09_raw}
FINAL_RELAXED=${final_chr09_relaxed}
FINAL_STRICT=${final_chr09_strict}

[reports]
JOINTCALL_SAMPLE_AUDIT=${JOINTCALL_SAMPLE_AUDIT}
CHR09_QC=${CHR09_QC}
GENOME_QC=${GENOME_QC}
CHR09_STATS=${CHR09_STATS}
GENOME_STATS=${GENOME_STATS}
WARNINGS=${WARNINGS_FILE}
RUN_LOG=${RUN_LOG}
EOF
}

cleanup() {
    rm -rf "$TMP_DIR"
    rm -rf "$SCRATCH_ROOT"
}

main() {
    mkdir -p "$OUTDIR" "$TMP_DIR" "$QC_DIR" "$LOG_DIR"
    : > "$WARNINGS_FILE"
    printf "scope\ttier\trecords\tbiallelic_snps\tbins_with_biallelic_snps\tfirst_pos\tlast_pos\ttimestamp\n" > "$CHR09_QC"
    printf "scope\ttier\trecords\tbiallelic_snps\tbins_with_biallelic_snps\tfirst_pos\tlast_pos\ttimestamp\n" > "$GENOME_QC"

    exec > >(tee -a "$RUN_LOG") 2>&1

    log "Starting full-Chr09 rerun + genome repair workflow"

    need_cmd bcftools
    need_cmd tabix
    need_cmd bgzip
    need_cmd awk
    need_cmd sort
    need_cmd diff
    need_cmd tee
    need_cmd gatk
    need_cmd cp
    need_cmd grep

    assert_file "$BASE_VCF"
    assert_file "$REFERENCE"
    assert_file "$REFERENCE_FAI"
    assert_file "$REFERENCE_DICT"
    assert_file "$PREPARED_GVCF_LIST"

    ensure_index "$BASE_VCF"
    prepare_input_list
    audit_prepared_gvcfs

    SCRATCH_ROOT="$(choose_scratch_root "chr09_full_12samples_repair")"
    JOINT_DB_DIR="${SCRATCH_ROOT}/Chr09.12samples.workspace"
    JOINT_TMP_DIR="${SCRATCH_ROOT}/Chr09.12samples.tmp"
    mkdir -p "$SCRATCH_ROOT"

    run_joint_call_for_chr09
    normalize_chr09_to_base_order

    local chr09_gt="no"
    local chr09_fmt_dp="no"
    local chr09_ad="no"
    local chr09_info_dp="no"
    local chr09_info_mq="no"

    has_header_id "$TMP_CHR09_MERGED" "FORMAT" "GT" && chr09_gt="yes"
    has_header_id "$TMP_CHR09_MERGED" "FORMAT" "DP" && chr09_fmt_dp="yes"
    has_header_id "$TMP_CHR09_MERGED" "FORMAT" "AD" && chr09_ad="yes"
    has_header_id "$TMP_CHR09_MERGED" "INFO" "DP" && chr09_info_dp="yes"
    has_header_id "$TMP_CHR09_MERGED" "INFO" "MQ" && chr09_info_mq="yes"

    [[ "$chr09_gt" == "yes" ]] || die "New 12-sample Chr09 VCF is missing FORMAT/GT"
    [[ "$chr09_fmt_dp" == "yes" ]] || die "New 12-sample Chr09 VCF is missing FORMAT/DP"
    [[ "$chr09_info_dp" == "yes" ]] || die "New 12-sample Chr09 VCF is missing INFO/DP"
    if [[ "$chr09_ad" != "yes" ]]; then
        warn "New 12-sample Chr09 VCF does not contain FORMAT/AD; check downstream SNP-index script if it depends on AD"
    fi
    if [[ "$chr09_info_mq" != "yes" ]]; then
        warn "New 12-sample Chr09 VCF does not contain INFO/MQ; strict SNP-index filter will be skipped"
    fi

    compute_qc_metrics "$TMP_CHR09_MERGED" "$CHR09_REGION" "raw" "$CHR09_QC"
    run_bcftools_stats "$TMP_CHR09_MERGED" "$CHR09_STATS"

    local chr09_strict_snps="0"
    local chr09_relaxed_snps="0"
    local chr09_strict_bins="0"
    local chr09_relaxed_bins="0"
    local chr09_raw_snps

    if [[ "$chr09_info_mq" == "yes" ]]; then
        compute_qc_metrics "$TMP_CHR09_MERGED" "$CHR09_REGION" "strict" "$CHR09_QC" "$STRICT_EXPR"
        chr09_strict_snps="$(count_biallelic_snps "$TMP_CHR09_MERGED" "$CHR09_REGION" "$STRICT_EXPR")"
        chr09_strict_bins="$(count_bins_with_biallelic_snps "$TMP_CHR09_MERGED" "$CHR09_REGION" "$STRICT_EXPR")"
    else
        append_qc_row "$CHR09_QC" "$CHR09_REGION" "strict_skipped_no_mq" "NA" "NA" "NA" "NA" "NA"
    fi

    compute_qc_metrics "$TMP_CHR09_MERGED" "$CHR09_REGION" "relaxed" "$CHR09_QC" "$RELAXED_EXPR"
    chr09_relaxed_snps="$(count_biallelic_snps "$TMP_CHR09_MERGED" "$CHR09_REGION" "$RELAXED_EXPR")"
    chr09_relaxed_bins="$(count_bins_with_biallelic_snps "$TMP_CHR09_MERGED" "$CHR09_REGION" "$RELAXED_EXPR")"
    chr09_raw_snps="$(count_biallelic_snps "$TMP_CHR09_MERGED" "$CHR09_REGION")"

    if [[ "$chr09_raw_snps" -gt 0 ]] && {
        [[ "$chr09_relaxed_snps" -gt 0 && "$chr09_relaxed_bins" -gt 0 ]] ||
        [[ "$chr09_strict_snps" -gt 0 && "$chr09_strict_bins" -gt 0 ]]
    }; then
        CHR09_SNPINDEX_CHECK="PASS"
    else
        die "New 12-sample Chr09 VCF does not pass SNP-index usability checks"
    fi

    log "Extracting non-Chr09 regions from base VCF"
    bcftools view -r "$PRE09_REGION" -Oz --threads "$THREADS" -o "$TMP_PART1" "$BASE_VCF"
    bcftools view -r "$POST09_REGION" -Oz --threads "$THREADS" -o "$TMP_PART3" "$BASE_VCF"
    ensure_index "$TMP_PART1"
    ensure_index "$TMP_PART3"

    printf "%s\n%s\n%s\n" "$TMP_PART1" "$TMP_CHR09_MERGED" "$TMP_PART3" > "${TMP_DIR}/genome_inputs.list"
    bcftools concat -f "${TMP_DIR}/genome_inputs.list" -Oz --threads "$THREADS" -o "$FINAL_VCF"
    ensure_index "$FINAL_VCF"

    local final_samples="${TMP_DIR}/final.samples.txt"
    extract_samples "$FINAL_VCF" "$final_samples"
    diff -u "$TMP_BASE_SAMPLES" "$final_samples" >/dev/null || die "Final genome sample order differs from base VCF"

    local genome_gt="no"
    local genome_fmt_dp="no"
    local genome_ad="no"
    local genome_info_dp="no"
    local genome_info_mq="no"

    has_header_id "$FINAL_VCF" "FORMAT" "GT" && genome_gt="yes"
    has_header_id "$FINAL_VCF" "FORMAT" "DP" && genome_fmt_dp="yes"
    has_header_id "$FINAL_VCF" "FORMAT" "AD" && genome_ad="yes"
    has_header_id "$FINAL_VCF" "INFO" "DP" && genome_info_dp="yes"
    has_header_id "$FINAL_VCF" "INFO" "MQ" && genome_info_mq="yes"

    [[ "$genome_gt" == "yes" ]] || die "Final genome VCF is missing FORMAT/GT"
    [[ "$genome_fmt_dp" == "yes" ]] || die "Final genome VCF is missing FORMAT/DP"
    [[ "$genome_info_dp" == "yes" ]] || die "Final genome VCF is missing INFO/DP"
    if [[ "$genome_ad" != "yes" ]]; then
        warn "Final genome VCF does not contain FORMAT/AD; check downstream SNP-index script if it depends on AD"
    fi
    if [[ "$genome_info_mq" != "yes" ]]; then
        warn "Final genome VCF does not contain INFO/MQ; strict genome SNP-index filter will be skipped"
    fi

    compute_qc_metrics "$FINAL_VCF" "$CHR09_REGION" "chr09_raw" "$GENOME_QC"
    compute_qc_metrics "$FINAL_VCF" "$GENOME_REGION" "genome_raw" "$GENOME_QC"
    run_bcftools_stats "$FINAL_VCF" "$GENOME_STATS"

    local genome_strict_snps="0"
    local genome_relaxed_snps="0"
    local final_chr09_strict="0"
    local final_chr09_relaxed="0"
    local final_chr09_raw
    local base_pre09
    local final_pre09
    local base_post09
    local final_post09

    if [[ "$genome_info_mq" == "yes" ]]; then
        compute_qc_metrics "$FINAL_VCF" "$GENOME_REGION" "genome_strict" "$GENOME_QC" "$STRICT_EXPR"
        compute_qc_metrics "$FINAL_VCF" "$CHR09_REGION" "chr09_strict" "$GENOME_QC" "$STRICT_EXPR"
        genome_strict_snps="$(count_biallelic_snps "$FINAL_VCF" "$GENOME_REGION" "$STRICT_EXPR")"
        final_chr09_strict="$(count_biallelic_snps "$FINAL_VCF" "$CHR09_REGION" "$STRICT_EXPR")"
    else
        append_qc_row "$GENOME_QC" "$GENOME_REGION" "genome_strict_skipped_no_mq" "NA" "NA" "NA" "NA" "NA"
        append_qc_row "$GENOME_QC" "$CHR09_REGION" "chr09_strict_skipped_no_mq" "NA" "NA" "NA" "NA" "NA"
    fi

    compute_qc_metrics "$FINAL_VCF" "$GENOME_REGION" "genome_relaxed" "$GENOME_QC" "$RELAXED_EXPR"
    compute_qc_metrics "$FINAL_VCF" "$CHR09_REGION" "chr09_relaxed" "$GENOME_QC" "$RELAXED_EXPR"
    genome_relaxed_snps="$(count_biallelic_snps "$FINAL_VCF" "$GENOME_REGION" "$RELAXED_EXPR")"
    final_chr09_relaxed="$(count_biallelic_snps "$FINAL_VCF" "$CHR09_REGION" "$RELAXED_EXPR")"
    final_chr09_raw="$(count_biallelic_snps "$FINAL_VCF" "$CHR09_REGION")"

    base_pre09="$(count_records "$BASE_VCF" "$PRE09_REGION")"
    final_pre09="$(count_records "$FINAL_VCF" "$PRE09_REGION")"
    base_post09="$(count_records "$BASE_VCF" "$POST09_REGION")"
    final_post09="$(count_records "$FINAL_VCF" "$POST09_REGION")"

    [[ "$base_pre09" == "$final_pre09" ]] || die "Chr01-08 record counts changed after repair"
    [[ "$base_post09" == "$final_post09" ]] || die "Chr10-12 record counts changed after repair"
    [[ "$final_chr09_raw" -gt 0 ]] || die "Final genome still has zero biallelic SNPs in ${CHR09_REGION}"

    if {
        [[ "$genome_relaxed_snps" -gt 0 && "$final_chr09_relaxed" -gt 0 ]] ||
        [[ "$genome_strict_snps" -gt 0 && "$final_chr09_strict" -gt 0 ]]
    }; then
        GENOME_SNPINDEX_CHECK="PASS"
    else
        die "Final genome VCF does not pass SNP-index usability checks"
    fi

    bcftools view -H -r "Chr09:1-1000000" "$FINAL_VCF" | awk 'NR == 1 { found = 1 } END { exit(found ? 0 : 1) }' \
        || die "Final genome query check failed for Chr09:1-1000000"

    write_summary \
        "$chr09_gt" "$chr09_fmt_dp" "$chr09_ad" "$chr09_info_dp" "$chr09_info_mq" \
        "$genome_gt" "$genome_fmt_dp" "$genome_ad" "$genome_info_dp" "$genome_info_mq" \
        "$base_pre09" "$final_pre09" "$base_post09" "$final_post09" \
        "$final_chr09_raw" "$final_chr09_relaxed" "$final_chr09_strict"

    cleanup
    log "Workflow completed successfully"
    log "Final VCF: $FINAL_VCF"
    log "Summary: $FINAL_SUMMARY"
}

main "$@"
