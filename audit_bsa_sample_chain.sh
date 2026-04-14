#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

usage() {
    cat <<'EOF'
Usage:
  bash audit_bsa_sample_chain.sh [base_dir] [expected_parents_csv] [outdir]

Defaults:
  base_dir             /data2/chenh/bsa/CleanData
  expected_parents_csv C3_Parent,C6_Parent
  outdir               <base_dir>/sample_chain_audit_<timestamp>

What it checks:
  1. fq pair completeness under fq_data/
  2. BAM file names vs BAM read-group SM tags
  3. gVCF file names vs header sample names
  4. gVCF list contents and sample names
  5. parent VCF / final BSA VCF sample names
  6. whether expected parent sample names are present or missing

Examples:
  bash audit_bsa_sample_chain.sh
  bash audit_bsa_sample_chain.sh /data2/chenh/bsa/CleanData C3_Parent,C6_Parent
  bash audit_bsa_sample_chain.sh /data2/chenh/bsa/CleanData C3_Parent,C6_Parent /data2/chenh/bsa/CleanData/audit_parent_check
EOF
}

if [[ $# -gt 3 ]]; then
    usage
    exit 1
fi

BASE_DIR="${1:-/data2/chenh/bsa/CleanData}"
EXPECTED_PARENTS_CSV="${2:-C3_Parent,C6_Parent}"
OUTDIR="${3:-${BASE_DIR%/}/sample_chain_audit_$(date '+%Y%m%d_%H%M%S')}"

FQ_DIR="${BASE_DIR%/}/fq_data"
BAM_DIR="${BASE_DIR%/}/bam"
GVCF_DIR="${BASE_DIR%/}/gvcf"
VCF_DIR="${BASE_DIR%/}/vcf"

LOG_FILE="${OUTDIR}/audit.log"
FQ_TSV="${OUTDIR}/fq_pairs.tsv"
BAM_TSV="${OUTDIR}/bam_samples.tsv"
GVCF_TSV="${OUTDIR}/gvcf_samples.tsv"
LIST_TSV="${OUTDIR}/gvcf_list_audit.tsv"
VCF_TSV="${OUTDIR}/vcf_samples.tsv"
SUMMARY_TXT="${OUTDIR}/summary.txt"
MATRIX_TSV="${OUTDIR}/sample_presence_matrix.tsv"
TMP_FQ_SAMPLES="${OUTDIR}/.tmp_fq_samples.txt"
TMP_BAM_SAMPLES="${OUTDIR}/.tmp_bam_samples.txt"
TMP_GVCF_SAMPLES="${OUTDIR}/.tmp_gvcf_samples.txt"
TMP_LIST_SAMPLES="${OUTDIR}/.tmp_list_samples.txt"
TMP_PARENT_VCF_SAMPLES="${OUTDIR}/.tmp_parent_vcf_samples.txt"
TMP_FINAL_VCF_SAMPLES="${OUTDIR}/.tmp_final_vcf_samples.txt"
TMP_JOINT_VCF_SAMPLES="${OUTDIR}/.tmp_joint_vcf_samples.txt"

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        printf 'ERROR: required command not found: %s\n' "$1" >&2
        exit 1
    }
}

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

join_by() {
    local sep="$1"
    shift || true
    local first=1
    local item
    for item in "$@"; do
        if [[ $first -eq 1 ]]; then
            printf '%s' "$item"
            first=0
        else
            printf '%s%s' "$sep" "$item"
        fi
    done
}

extract_vcf_samples() {
    local path="$1"
    bcftools query -l "$path" 2>/dev/null || true
}

extract_bam_rg_sm() {
    local path="$1"
    samtools view -H "$path" 2>/dev/null \
        | awk -F'\t' '/^@RG/ {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }' \
        | sort -u
}

file_basename_without_ext() {
    local name
    name="$(basename "$1")"
    name="${name%.clean.fq.gz}"
    name="${name%.fq.gz}"
    name="${name%.bam}"
    name="${name%.g.vcf.gz}"
    name="${name%.vcf.gz}"
    name="${name%.vcf}"
    printf '%s\n' "$name"
}

mkdir -p "$OUTDIR"
exec > >(tee -a "$LOG_FILE") 2>&1

need_cmd awk
need_cmd sort
need_cmd find
need_cmd grep
need_cmd samtools
need_cmd bcftools

log "BASE_DIR              = $BASE_DIR"
log "EXPECTED_PARENTS_CSV  = $EXPECTED_PARENTS_CSV"
log "OUTDIR                = $OUTDIR"

printf "sample\tr1_path\tr2_path\tstatus\n" > "$FQ_TSV"
printf "bam_path\tbam_basename\trg_sm_count\trg_sm_names\tstatus\n" > "$BAM_TSV"
printf "gvcf_path\tgvcf_basename\tsample_count\tsample_names\tstatus\n" > "$GVCF_TSV"
printf "list_file\tlisted_path\tsample_count\tsample_names\tstatus\n" > "$LIST_TSV"
printf "vcf_path\tcategory\tsample_count\tsample_names\tparent_check\n" > "$VCF_TSV"
printf "sample\tin_fq\tin_bam\tin_gvcf\tin_gvcf_list\tin_parent_vcf\tin_final_bsa_vcf\tin_joint_vcf\n" > "$MATRIX_TSV"

: > "$TMP_FQ_SAMPLES"
: > "$TMP_BAM_SAMPLES"
: > "$TMP_GVCF_SAMPLES"
: > "$TMP_LIST_SAMPLES"
: > "$TMP_PARENT_VCF_SAMPLES"
: > "$TMP_FINAL_VCF_SAMPLES"
: > "$TMP_JOINT_VCF_SAMPLES"

EXPECTED_PARENT_LIST=()
IFS=',' read -r -a EXPECTED_PARENT_LIST <<< "$EXPECTED_PARENTS_CSV"

log "Step 1/5: audit fq pairs"
if [[ -d "$FQ_DIR" ]]; then
    declare -A fq_r1=()
    declare -A fq_r2=()
    for fq in "$FQ_DIR"/*_1.clean.fq.gz; do
        [[ -e "$fq" ]] || continue
        sample="$(basename "$fq")"
        sample="${sample%_1.clean.fq.gz}"
        fq_r1["$sample"]="$fq"
    done
    for fq in "$FQ_DIR"/*_2.clean.fq.gz; do
        [[ -e "$fq" ]] || continue
        sample="$(basename "$fq")"
        sample="${sample%_2.clean.fq.gz}"
        fq_r2["$sample"]="$fq"
    done

    samples="$(
        {
            for sample in "${!fq_r1[@]}"; do
                printf '%s\n' "$sample"
            done
            for sample in "${!fq_r2[@]}"; do
                printf '%s\n' "$sample"
            done
        } | awk 'NF' | sort -u
    )"
    while IFS= read -r sample; do
        [[ -n "$sample" ]] || continue
        r1="${fq_r1[$sample]-}"
        r2="${fq_r2[$sample]-}"
        status="ok"
        [[ -n "$r1" ]] || status="missing_r1"
        [[ -n "$r2" ]] || status="missing_r2"
        if [[ -z "$r1" && -z "$r2" ]]; then
            status="missing_both"
        fi
        printf "%s\t%s\t%s\t%s\n" "$sample" "$r1" "$r2" "$status" >> "$FQ_TSV"
        printf "%s\n" "$sample" >> "$TMP_FQ_SAMPLES"
    done <<< "$samples"
else
    log "fq directory not found, skip: $FQ_DIR"
fi

log "Step 2/5: audit bam sample names"
if [[ -d "$BAM_DIR" ]]; then
    for bam in "$BAM_DIR"/*.bam; do
        [[ -e "$bam" ]] || continue
        bam_base="$(file_basename_without_ext "$bam")"
        mapfile -t sm_names < <(extract_bam_rg_sm "$bam")
        sm_count="${#sm_names[@]}"
        sm_csv="$(join_by "," "${sm_names[@]}")"
        status="ok"
        if [[ "$sm_count" -eq 0 ]]; then
            status="missing_rg_sm"
        elif [[ "$sm_count" -gt 1 ]]; then
            status="multiple_rg_sm"
        elif [[ "${sm_names[0]}" != "$bam_base" ]]; then
            status="basename_sm_mismatch"
        fi
        printf "%s\t%s\t%s\t%s\t%s\n" "$bam" "$bam_base" "$sm_count" "$sm_csv" "$status" >> "$BAM_TSV"
        printf "%s\n" "${sm_names[@]}" >> "$TMP_BAM_SAMPLES"
    done
else
    log "bam directory not found, skip: $BAM_DIR"
fi

log "Step 3/5: audit gVCF sample names"
if [[ -d "$GVCF_DIR" ]]; then
    for gvcf in "$GVCF_DIR"/*.g.vcf.gz; do
        [[ -e "$gvcf" ]] || continue
        gvcf_base="$(file_basename_without_ext "$gvcf")"
        mapfile -t g_samples < <(extract_vcf_samples "$gvcf")
        sample_count="${#g_samples[@]}"
        sample_csv="$(join_by "," "${g_samples[@]}")"
        status="ok"
        if [[ "$sample_count" -eq 0 ]]; then
            status="missing_header_sample"
        elif [[ "$sample_count" -gt 1 ]]; then
            status="multi_sample_gvcf"
        elif [[ "${g_samples[0]}" != "$gvcf_base" ]]; then
            status="basename_sample_mismatch"
        fi
        printf "%s\t%s\t%s\t%s\t%s\n" "$gvcf" "$gvcf_base" "$sample_count" "$sample_csv" "$status" >> "$GVCF_TSV"
        printf "%s\n" "${g_samples[@]}" >> "$TMP_GVCF_SAMPLES"
    done
else
    log "gvcf directory not found, skip: $GVCF_DIR"
fi

log "Step 4/5: audit gVCF lists"
list_candidates=()
[[ -f "${VCF_DIR}/gvcf.list" ]] && list_candidates+=("${VCF_DIR}/gvcf.list")
[[ -f "${VCF_DIR}/hunchi_and_qinben.gvcf.list" ]] && list_candidates+=("${VCF_DIR}/hunchi_and_qinben.gvcf.list")
[[ -f "${BASE_DIR}/dbimport_vcf/gvcf.list" ]] && list_candidates+=("${BASE_DIR}/dbimport_vcf/gvcf.list")

for list_file in "${list_candidates[@]}"; do
    while IFS= read -r listed_path; do
        [[ -n "$listed_path" ]] || continue
        [[ "$listed_path" =~ ^# ]] && continue
        if [[ ! -f "$listed_path" ]]; then
            printf "%s\t%s\t0\t\tmissing_file\n" "$list_file" "$listed_path" >> "$LIST_TSV"
            continue
        fi
        mapfile -t list_samples < <(extract_vcf_samples "$listed_path")
        sample_count="${#list_samples[@]}"
        sample_csv="$(join_by "," "${list_samples[@]}")"
        status="ok"
        [[ "$sample_count" -gt 0 ]] || status="missing_header_sample"
        if [[ "$sample_count" -gt 1 ]]; then
            status="multi_sample_input"
        fi
        printf "%s\t%s\t%s\t%s\t%s\n" "$list_file" "$listed_path" "$sample_count" "$sample_csv" "$status" >> "$LIST_TSV"
        printf "%s\n" "${list_samples[@]}" >> "$TMP_LIST_SAMPLES"
    done < "$list_file"
done

log "Step 5/5: audit parent/final VCF sample names"
candidate_vcfs=()
while IFS= read -r path; do
    candidate_vcfs+=("$path")
done < <(
    find "$BASE_DIR" -type f \
        \( -name 'Parents_Combined.g.vcf.gz' \
        -o -name 'Parents_Final_Raw.vcf.gz' \
        -o -name 'Parents_Filtered.vcf.gz' \
        -o -name 'BSA_FINAL_RESULT_COMPLETE.vcf.gz' \
        -o -name 'Pepper_All_Chr_Whole_Genome.vcf.gz' \
        -o -name 'Pepper_All_Chr_Whole_Genome.vcf' \
        -o -name 'BSA_FINAL_ULTIMATE_FIXED.vcf.gz' \
        -o -name 'BSA_READY_FOR_PYTHON.vcf.gz' \
        -o -name 'BSA_READY_V2_RELAXED.vcf.gz' \
        -o -name '*.genome.joint.vcf.gz' \) \
        | sort -u
)

for vcf in "${candidate_vcfs[@]}"; do
    category="other"
    base="$(basename "$vcf")"
    case "$base" in
        Parents_Combined.g.vcf.gz) category="parent_gvcf" ;;
        Parents_Final_Raw.vcf.gz|Parents_Filtered.vcf.gz) category="parent_vcf" ;;
        BSA_FINAL_RESULT_COMPLETE.vcf.gz|Pepper_All_Chr_Whole_Genome.vcf.gz|Pepper_All_Chr_Whole_Genome.vcf|BSA_FINAL_ULTIMATE_FIXED.vcf.gz|BSA_READY_FOR_PYTHON.vcf.gz|BSA_READY_V2_RELAXED.vcf.gz) category="final_bsa_vcf" ;;
        *.genome.joint.vcf.gz) category="joint_vcf" ;;
    esac

    mapfile -t vcf_samples < <(extract_vcf_samples "$vcf")
    sample_count="${#vcf_samples[@]}"
    sample_csv="$(join_by "," "${vcf_samples[@]}")"
    parent_check="na"

    if [[ "$category" == "parent_gvcf" || "$category" == "parent_vcf" || "$category" == "final_bsa_vcf" || "$category" == "joint_vcf" ]]; then
        missing=()
        for parent in "${EXPECTED_PARENT_LIST[@]}"; do
            [[ -n "$parent" ]] || continue
            found=0
            for sample in "${vcf_samples[@]}"; do
                if [[ "$sample" == "$parent" ]]; then
                    found=1
                    break
                fi
            done
            [[ "$found" -eq 1 ]] || missing+=("$parent")
        done
        if [[ "${#missing[@]}" -eq 0 ]]; then
            parent_check="expected_parents_present"
        else
            parent_check="missing:$(join_by "," "${missing[@]}")"
        fi
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" "$vcf" "$category" "$sample_count" "$sample_csv" "$parent_check" >> "$VCF_TSV"
    case "$category" in
        parent_gvcf|parent_vcf)
            printf "%s\n" "${vcf_samples[@]}" >> "$TMP_PARENT_VCF_SAMPLES"
            ;;
        final_bsa_vcf)
            printf "%s\n" "${vcf_samples[@]}" >> "$TMP_FINAL_VCF_SAMPLES"
            ;;
        joint_vcf)
            printf "%s\n" "${vcf_samples[@]}" >> "$TMP_JOINT_VCF_SAMPLES"
            ;;
    esac
done

all_samples="$(
    cat "$TMP_FQ_SAMPLES" "$TMP_BAM_SAMPLES" "$TMP_GVCF_SAMPLES" "$TMP_LIST_SAMPLES" \
        "$TMP_PARENT_VCF_SAMPLES" "$TMP_FINAL_VCF_SAMPLES" "$TMP_JOINT_VCF_SAMPLES" \
        | awk 'NF' | sort -u
)"

while IFS= read -r sample; do
    [[ -n "$sample" ]] || continue
    in_fq=0
    in_bam=0
    in_gvcf=0
    in_list=0
    in_parent=0
    in_final=0
    in_joint=0

    grep -Fxq "$sample" "$TMP_FQ_SAMPLES" && in_fq=1 || true
    grep -Fxq "$sample" "$TMP_BAM_SAMPLES" && in_bam=1 || true
    grep -Fxq "$sample" "$TMP_GVCF_SAMPLES" && in_gvcf=1 || true
    grep -Fxq "$sample" "$TMP_LIST_SAMPLES" && in_list=1 || true
    grep -Fxq "$sample" "$TMP_PARENT_VCF_SAMPLES" && in_parent=1 || true
    grep -Fxq "$sample" "$TMP_FINAL_VCF_SAMPLES" && in_final=1 || true
    grep -Fxq "$sample" "$TMP_JOINT_VCF_SAMPLES" && in_joint=1 || true

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$in_fq" "$in_bam" "$in_gvcf" "$in_list" "$in_parent" "$in_final" "$in_joint" \
        >> "$MATRIX_TSV"
done <<< "$all_samples"

{
    echo "BSA sample chain audit summary"
    echo
    echo "Base dir: $BASE_DIR"
    echo "Expected parents: $EXPECTED_PARENTS_CSV"
    echo
    echo "[fq pair issues]"
    awk -F'\t' 'NR==1 || $4 != "ok"' "$FQ_TSV"
    echo
    echo "[bam issues]"
    awk -F'\t' 'NR==1 || $5 != "ok"' "$BAM_TSV"
    echo
    echo "[gvcf issues]"
    awk -F'\t' 'NR==1 || $5 != "ok"' "$GVCF_TSV"
    echo
    echo "[gvcf list issues]"
    awk -F'\t' 'NR==1 || $5 != "ok"' "$LIST_TSV"
    echo
    echo "[vcf parent check]"
    awk -F'\t' 'NR==1 || $5 != "expected_parents_present"' "$VCF_TSV"
    echo
    echo "[sample presence matrix: suspicious rows]"
    awk -F'\t' 'NR==1 || $2+$3+$4+$5+$6+$7+$8 < 4' "$MATRIX_TSV"
} > "$SUMMARY_TXT"

log "Audit completed"
log "fq report      = $FQ_TSV"
log "bam report     = $BAM_TSV"
log "gvcf report    = $GVCF_TSV"
log "list report    = $LIST_TSV"
log "vcf report     = $VCF_TSV"
log "matrix         = $MATRIX_TSV"
log "summary        = $SUMMARY_TXT"
