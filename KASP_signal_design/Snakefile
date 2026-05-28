configfile: "data/config.yaml"

from pathlib import Path
import csv


SCRIPT_DIR = Path("scripts")
TEMP_DIR = Path("temp")
RESULTS_DIR = Path("results")
RULE_ENV = "data/kasp_rule_env.yaml"
REGIONS_FILE = config.get("regions_file", "data/regions.tsv")


def load_regions(path: str) -> dict[str, str]:
    regions: dict[str, str] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            region_id = row["region_id"].strip()
            region = row["region"].strip()
            if not region_id or not region:
                continue
            regions[region_id] = region
    if not regions:
        raise ValueError("data/regions.tsv 中没有可用区间")
    return regions


REGIONS = load_regions(REGIONS_FILE)
REGION_IDS = list(REGIONS.keys())


def region_value(wildcards):
    return REGIONS[wildcards.region_id]


rule all:
    input:
        expand("results/{region_id}/04_kasp_sequences.tsv", region_id=REGION_IDS)


rule rank_candidates:
    input:
        snpindex=lambda wildcards: config["snpindex_path"]
    output:
        candidates="results/{region_id}/01_ranked_candidates.tsv",
        summary="results/{region_id}/01_rank_summary.tsv"
    params:
        region=region_value,
        rank_limit=lambda wildcards: config["initial_rank_limit"],
        expansion_step=lambda wildcards: config["expansion_step_bp"],
        selection_mode=lambda wildcards: config.get("selection_mode", "flank"),
        bin_size=lambda wildcards: config.get("bin_size_bp", 1000000),
        log_dir=lambda wildcards: f"temp/{wildcards.region_id}"
    log:
        "temp/{region_id}/01_rank_candidates.log"
    conda:
        RULE_ENV
    shell:
        r"""
        mkdir -p "{params.log_dir}" "$(dirname "{output.candidates}")"
        python "{SCRIPT_DIR}/20260418_kasp_01_rank_snpindex_candidates.py" \
          --region "{params.region}" \
          --snpindex "{input.snpindex}" \
          --output-candidates "{output.candidates}" \
          --output-summary "{output.summary}" \
          --initial-rank-limit "{params.rank_limit}" \
          --expansion-step-bp "{params.expansion_step}" \
          --selection-mode "{params.selection_mode}" \
          --bin-size-bp "{params.bin_size}" \
          > "{log}" 2>&1
        """


rule filter_sites:
    input:
        candidates="results/{region_id}/01_ranked_candidates.tsv",
        vcf=lambda wildcards: config["vcf_path"]
    output:
        checked="results/{region_id}/02_checked_sites.tsv",
        passed="results/{region_id}/02_pass_sites.tsv"
    params:
        parent1=lambda wildcards: config["parent1_sample"],
        parent2=lambda wildcards: config["parent2_sample"],
        high=lambda wildcards: config["high_pool_sample"],
        low=lambda wildcards: config["low_pool_sample"],
        flank=lambda wildcards: config["flank_no_variant_bp"],
        target=lambda wildcards: config["target_markers_per_side"],
        rank_limit=lambda wildcards: config["initial_rank_limit"],
        log_dir=lambda wildcards: f"temp/{wildcards.region_id}"
    log:
        "temp/{region_id}/02_filter_sites.log"
    conda:
        RULE_ENV
    shell:
        r"""
        mkdir -p "{params.log_dir}" "$(dirname "{output.checked}")"
        python "{SCRIPT_DIR}/20260418_kasp_02_filter_isolated_parental_snps.py" \
          --candidates "{input.candidates}" \
          --vcf "{input.vcf}" \
          --parent1 "{params.parent1}" \
          --parent2 "{params.parent2}" \
          --high-pool "{params.high}" \
          --low-pool "{params.low}" \
          --flank-bp "{params.flank}" \
          --target-markers-per-side "{params.target}" \
          --initial-rank-limit "{params.rank_limit}" \
          --output-checked "{output.checked}" \
          --output-pass "{output.passed}" \
          > "{log}" 2>&1
        """


rule prepare_extract_input:
    input:
        passed="results/{region_id}/02_pass_sites.tsv",
        checked="results/{region_id}/02_checked_sites.tsv"
    output:
        extract_input="results/{region_id}/03_extract_fa_input.tsv",
        recommended="results/{region_id}/03_recommended_markers.tsv"
    params:
        target=lambda wildcards: config["target_markers_per_side"],
        log_dir=lambda wildcards: f"temp/{wildcards.region_id}"
    log:
        "temp/{region_id}/03_prepare_extract_input.log"
    conda:
        RULE_ENV
    shell:
        r"""
        mkdir -p "{params.log_dir}" "$(dirname "{output.extract_input}")"
        python "{SCRIPT_DIR}/20260418_kasp_03_prepare_extract_fa_input.py" \
          --pass-sites "{input.passed}" \
          --checked-sites "{input.checked}" \
          --output-extract-input "{output.extract_input}" \
          --output-recommended "{output.recommended}" \
          --target-markers-per-side "{params.target}" \
          > "{log}" 2>&1
        """


rule extract_fa:
    input:
        extract_input="results/{region_id}/03_extract_fa_input.tsv"
    output:
        kasp="results/{region_id}/04_kasp_sequences.tsv"
    params:
        ref=lambda wildcards: config["reference_fasta"],
        log_dir=lambda wildcards: f"temp/{wildcards.region_id}"
    log:
        "temp/{region_id}/04_extract_fa.log"
    conda:
        RULE_ENV
    shell:
        r"""
        mkdir -p "{params.log_dir}" "$(dirname "{output.kasp}")"
        python "{SCRIPT_DIR}/extract_fa.py" \
          --input "{input.extract_input}" \
          --output "{output.kasp}" \
          --ref "{params.ref}" \
          > "{log}" 2>&1
        """
