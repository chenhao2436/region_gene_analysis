#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gvcf_plidy_30_BSA.py

目标
1. 读取 joint VCF / VCF.GZ
2. 用高池、低池、双亲 4 个样本做 pooled BSA 的 SNP-index 计算
3. 输出 site 表、window 表、bsa.fig.py 兼容绘图表、Markdown 摘要

核心原理
- 只使用双亲方向明确的位点: high-parent 和 low-parent 必须是相反的纯合状态
- SNP-index 的方向统一定义为 high-parent 来源等位基因频率
- delta SNP-index = high_pool_snp_index - low_pool_snp_index
- window 阈值使用双总体比例差的正态近似标准误

兼容 bsa.fig.py 的输出
- plot.tsv 固定输出 7 列:
  CHROM POS THRESHOLD HIGH_SNP_INDEX LOW_SNP_INDEX DELTA_SNP_INDEX SNP_N
- 第一行以 # 开头, 原始 bsa.fig.py 会自动跳过
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import math
import os
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_OUTPUT_ROOT = "/data2/chenh/bsa/CleanData/gvcf_plidy_30/bsa"
Z_95 = 1.959963984540054
Z_99 = 2.5758293035489004
Z_999 = 3.2905267314919255

GT_HOM_REF = {"0/0", "0|0"}
GT_HOM_ALT = {"1/1", "1|1"}
GT_HET = {"0/1", "1/0", "0|1", "1|0"}
CHR_NUM_RE = re.compile(r"^(?:chr|Chr|CHR)?(\d+)$")


@dataclass
class SampleCall:
    gt: str
    dp: int
    ref_count: int
    alt_count: int


@dataclass
class SiteRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    parent1_gt: str
    parent2_gt: str
    high_gt: str
    low_gt: str
    parent1_state: int
    parent2_state: int
    parent2_allele_in_refalt: str
    high_dp: int
    low_dp: int
    high_parent2_count: int
    low_parent2_count: int
    high_snp_index: float
    low_snp_index: float
    delta_snp_index: float
    ed_value: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute pooled BSA SNP-index from a joint VCF generated after ploidy=30 calling."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("--high-bulk", required=True, help="High bulk sample name in VCF header")
    parser.add_argument("--low-bulk", required=True, help="Low bulk sample name in VCF header")
    parser.add_argument("--high-parent", dest="high_parent", default=None, help="Sample name of the parent that matches the high bulk allele direction")
    parser.add_argument("--low-parent", dest="low_parent", default=None, help="Sample name of the parent that matches the low bulk allele direction")
    parser.add_argument("--parent1", default=None, help="Legacy alias for --low-parent")
    parser.add_argument("--parent2", default=None, help="Legacy alias for --high-parent")
    parser.add_argument("--out-prefix", default="ploidy_30_3.18", help="Output prefix")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT, help="Structured output root")
    parser.add_argument("--window-size", type=int, default=1000000, help="Window size in bp")
    parser.add_argument("--step-size", type=int, default=100000, help="Window step in bp")
    parser.add_argument("--min-parent-dp", type=int, default=5, help="Minimum DP for each parent")
    parser.add_argument("--min-bulk-dp", type=int, default=10, help="Minimum DP for each bulk")
    parser.add_argument("--min-window-snps", type=int, default=15, help="Minimum SNP count for a usable window")
    parser.add_argument(
        "--parent-min-hom-af",
        type=float,
        default=0.9,
        help="Infer parent homozygous state from AD only if AF >= this cutoff",
    )
    parser.add_argument("--allow-parent-ad-rescue", action="store_true", help="Allow parent homozygous state rescue from AD")
    parser.add_argument("--allow-gt-count-fallback", action="store_true", help="Approximate counts from GT+DP if AD is missing")
    parser.add_argument(
        "--diagnostic-report",
        default=None,
        help="Write a chromosome-level diagnostic report TSV/MD prefix instead of failing when no usable sites are found",
    )
    return parser.parse_args()


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def safe_int(value: Optional[str], default: int = 0) -> int:
    if value in (None, "", "."):
        return default
    try:
        return int(value)
    except Exception:
        try:
            return int(float(value))
        except Exception:
            return default


def format_float(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.6f}"


def natural_chr_key(chrom: str):
    text = str(chrom).strip()
    match = CHR_NUM_RE.match(text)
    if match:
        return (0, int(match.group(1)), text)
    return (1, text)


def resolve_parent_names(args: argparse.Namespace) -> Tuple[str, str]:
    """
    Return (high_parent_name, low_parent_name).

    Priority:
    1) explicit --high-parent / --low-parent
    2) legacy --parent2 / --parent1
    """
    high_parent = args.high_parent or args.parent2
    low_parent = args.low_parent or args.parent1
    if not high_parent or not low_parent:
        raise SystemExit(
            "Please provide either --high-parent/--low-parent or the legacy --parent2/--parent1 pair."
        )
    return high_parent, low_parent


def parse_sample_call(sample_text: str, fmt_keys: Sequence[str], allow_gt_count_fallback: bool) -> SampleCall:
    values = sample_text.split(":")
    data = {key: values[idx] if idx < len(values) else "." for idx, key in enumerate(fmt_keys)}
    gt = data.get("GT", ".")
    dp = safe_int(data.get("DP"), 0)

    ref_count = 0
    alt_count = 0
    ad = data.get("AD", ".")
    if ad not in (None, "", "."):
        parts = ad.split(",")
        if len(parts) >= 2:
            ref_count = safe_int(parts[0], 0)
            alt_count = safe_int(parts[1], 0)
            if dp <= 0:
                dp = ref_count + alt_count

    if dp <= 0 and (ref_count > 0 or alt_count > 0):
        dp = ref_count + alt_count

    if dp > 0 and ref_count == 0 and alt_count == 0 and allow_gt_count_fallback:
        if gt in GT_HOM_REF:
            ref_count = dp
        elif gt in GT_HOM_ALT:
            alt_count = dp
        elif gt in GT_HET:
            ref_count = dp // 2
            alt_count = dp - ref_count

    return SampleCall(gt=gt, dp=dp, ref_count=ref_count, alt_count=alt_count)


def infer_parent_state(call: SampleCall, min_dp: int, min_hom_af: float, allow_ad_rescue: bool) -> Optional[int]:
    if call.dp < min_dp:
        return None
    if call.gt in GT_HOM_REF:
        return 0
    if call.gt in GT_HOM_ALT:
        return 1
    if call.gt in GT_HET:
        return None
    if not allow_ad_rescue:
        return None

    total = call.ref_count + call.alt_count
    if total < min_dp:
        return None
    ref_af = call.ref_count / total if total else 0.0
    alt_af = call.alt_count / total if total else 0.0
    if ref_af >= min_hom_af:
        return 0
    if alt_af >= min_hom_af:
        return 1
    return None


def compute_thresholds(target1: int, dp1: int, target2: int, dp2: int) -> Dict[str, float]:
    if dp1 <= 0 or dp2 <= 0:
        return {"thr95": math.nan, "thr99": math.nan, "thr999": math.nan}
    pooled = (target1 + target2) / (dp1 + dp2)
    variance = pooled * (1.0 - pooled) * (1.0 / dp1 + 1.0 / dp2)
    variance = max(variance, 1e-15)
    se = math.sqrt(variance)
    return {"thr95": Z_95 * se, "thr99": Z_99 * se, "thr999": Z_999 * se}


def process_vcf(args: argparse.Namespace) -> Tuple[Dict[str, List[SiteRecord]], Dict[str, int]]:
    counters: Counter = Counter()
    chrom_counters: Dict[str, Counter] = defaultdict(Counter)
    chr_sites: Dict[str, List[SiteRecord]] = defaultdict(list)
    sample_map: Dict[str, int] = {}
    high_parent_name, low_parent_name = resolve_parent_names(args)

    with open_text(args.vcf) as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                header = raw.rstrip("\n").split("\t")
                for idx, sample_name in enumerate(header[9:]):
                    sample_map[sample_name] = 9 + idx
                required = [args.high_bulk, args.low_bulk, high_parent_name, low_parent_name]
                missing = [name for name in required if name not in sample_map]
                if missing:
                    raise SystemExit("Samples not found in VCF header: " + ", ".join(missing))
                continue

            counters["total_records"] += 1
            row = raw.rstrip("\n").split("\t")
            if len(row) < 10:
                counters["short_rows"] += 1
                continue

            chrom, pos_s, _, ref, alt = row[0], row[1], row[2], row[3], row[4]
            chrom_counters[chrom]["total_records"] += 1
            if "," in alt:
                counters["multiallelic_rows"] += 1
                chrom_counters[chrom]["multiallelic_rows"] += 1
                continue
            if len(ref) != 1 or len(alt) != 1:
                counters["non_snp_rows"] += 1
                chrom_counters[chrom]["non_snp_rows"] += 1
                continue

            fmt_keys = row[8].split(":")
            pos = safe_int(pos_s, 0)
            high = parse_sample_call(row[sample_map[args.high_bulk]], fmt_keys, args.allow_gt_count_fallback)
            low = parse_sample_call(row[sample_map[args.low_bulk]], fmt_keys, args.allow_gt_count_fallback)
            high_parent = parse_sample_call(row[sample_map[high_parent_name]], fmt_keys, args.allow_gt_count_fallback)
            low_parent = parse_sample_call(row[sample_map[low_parent_name]], fmt_keys, args.allow_gt_count_fallback)

            high_parent_state = infer_parent_state(high_parent, args.min_parent_dp, args.parent_min_hom_af, args.allow_parent_ad_rescue)
            low_parent_state = infer_parent_state(low_parent, args.min_parent_dp, args.parent_min_hom_af, args.allow_parent_ad_rescue)
            if high_parent_state is None or low_parent_state is None:
                counters["parent_not_callable"] += 1
                chrom_counters[chrom]["parent_not_callable"] += 1
                continue
            if high_parent_state == low_parent_state:
                counters["parent_not_informative"] += 1
                chrom_counters[chrom]["parent_not_informative"] += 1
                continue
            if high.dp < args.min_bulk_dp or low.dp < args.min_bulk_dp:
                counters["bulk_low_dp"] += 1
                chrom_counters[chrom]["bulk_low_dp"] += 1
                continue

            if high_parent_state == 0:
                high_target = high.ref_count
                low_target = low.ref_count
                high_parent_allele = "REF"
            else:
                high_target = high.alt_count
                low_target = low.alt_count
                high_parent_allele = "ALT"

            high_index = high_target / high.dp
            low_index = low_target / low.dp
            delta = high_index - low_index
            ed_value = math.sqrt((high_index - low_index) ** 2 + ((1.0 - high_index) - (1.0 - low_index)) ** 2)

            chr_sites[chrom].append(
                SiteRecord(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    parent1_gt=low_parent.gt,
                    parent2_gt=high_parent.gt,
                    high_gt=high.gt,
                    low_gt=low.gt,
                    parent1_state=low_parent_state,
                    parent2_state=high_parent_state,
                    parent2_allele_in_refalt=high_parent_allele,
                    high_dp=high.dp,
                    low_dp=low.dp,
                    high_parent2_count=high_target,
                    low_parent2_count=low_target,
                    high_snp_index=high_index,
                    low_snp_index=low_index,
                    delta_snp_index=delta,
                    ed_value=ed_value,
                )
            )
            counters["usable_sites"] += 1
            chrom_counters[chrom]["usable_sites"] += 1

    if not sample_map:
        raise SystemExit("No #CHROM header found. Please provide a valid VCF/VCF.GZ.")

    for chrom in chr_sites:
        chr_sites[chrom].sort(key=lambda record: record.pos)
    return chr_sites, counters, chrom_counters


def write_diagnostic_report(path: str, chrom_counters: Dict[str, Counter], chroms: Sequence[str]):
    total = Counter()
    for chrom in chroms:
        total.update(chrom_counters.get(chrom, Counter()))

    with open(path, "w", encoding="utf-8") as out:
        out.write("CHROM\ttotal_records\tmultiallelic_rows\tnon_snp_rows\tparent_not_callable\tparent_not_informative\tbulk_low_dp\tusable_sites\n")
        for chrom in chroms:
            c = chrom_counters.get(chrom, Counter())
            out.write(
                f"{chrom}\t{c.get('total_records', 0)}\t{c.get('multiallelic_rows', 0)}\t{c.get('non_snp_rows', 0)}\t"
                f"{c.get('parent_not_callable', 0)}\t{c.get('parent_not_informative', 0)}\t"
                f"{c.get('bulk_low_dp', 0)}\t{c.get('usable_sites', 0)}\n"
            )
        out.write(
            f"TOTAL\t{total.get('total_records', 0)}\t{total.get('multiallelic_rows', 0)}\t{total.get('non_snp_rows', 0)}\t"
            f"{total.get('parent_not_callable', 0)}\t{total.get('parent_not_informative', 0)}\t"
            f"{total.get('bulk_low_dp', 0)}\t{total.get('usable_sites', 0)}\n"
        )


def summarise_window(records: Sequence[SiteRecord]) -> Dict[str, float]:
    if not records:
        return {
            "SNP_N": 0,
            "HIGH_TOTAL_DP": 0,
            "LOW_TOTAL_DP": 0,
            "HIGH_TOTAL_PARENT2_COUNT": 0,
            "LOW_TOTAL_PARENT2_COUNT": 0,
            "HIGH_SNP_INDEX": math.nan,
            "LOW_SNP_INDEX": math.nan,
            "DELTA_SNP_INDEX": math.nan,
            "MEAN_DELTA_SNP_INDEX": math.nan,
            "ED_MEAN": math.nan,
            "ABS_DELTA_THRESHOLD_95": math.nan,
            "ABS_DELTA_THRESHOLD_99": math.nan,
            "ABS_DELTA_THRESHOLD_999": math.nan,
        }

    high_target = sum(record.high_parent2_count for record in records)
    low_target = sum(record.low_parent2_count for record in records)
    high_dp = sum(record.high_dp for record in records)
    low_dp = sum(record.low_dp for record in records)
    thresholds = compute_thresholds(high_target, high_dp, low_target, low_dp)
    high_index = high_target / high_dp if high_dp else math.nan
    low_index = low_target / low_dp if low_dp else math.nan
    delta = high_index - low_index if high_dp and low_dp else math.nan

    return {
        "SNP_N": len(records),
        "HIGH_TOTAL_DP": high_dp,
        "LOW_TOTAL_DP": low_dp,
        "HIGH_TOTAL_PARENT2_COUNT": high_target,
        "LOW_TOTAL_PARENT2_COUNT": low_target,
        "HIGH_SNP_INDEX": high_index,
        "LOW_SNP_INDEX": low_index,
        "DELTA_SNP_INDEX": delta,
        "MEAN_DELTA_SNP_INDEX": statistics.fmean(record.delta_snp_index for record in records),
        "ED_MEAN": statistics.fmean(record.ed_value for record in records),
        "ABS_DELTA_THRESHOLD_95": thresholds["thr95"],
        "ABS_DELTA_THRESHOLD_99": thresholds["thr99"],
        "ABS_DELTA_THRESHOLD_999": thresholds["thr999"],
    }


def build_windows(chr_sites: Dict[str, List[SiteRecord]], args: argparse.Namespace) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for chrom in sorted(chr_sites, key=natural_chr_key):
        records = chr_sites[chrom]
        positions = [record.pos for record in records]
        max_pos = positions[-1]
        window_start = 1
        while window_start <= max_pos:
            window_end = window_start + args.window_size - 1
            window_center = window_start + (args.window_size // 2)
            left = bisect.bisect_left(positions, window_start)
            right = bisect.bisect_right(positions, window_end)
            chunk = records[left:right]
            summary = summarise_window(chunk)
            row: Dict[str, object] = {
                "CHROM": chrom,
                "WINDOW_START": window_start,
                "WINDOW_END": window_end,
                "WINDOW_CENTER": window_center,
                "USABLE_WINDOW": 1 if summary["SNP_N"] >= args.min_window_snps else 0,
            }
            row.update(summary)
            rows.append(row)
            window_start += args.step_size
    return rows


def write_sites(path: str, chr_sites: Dict[str, List[SiteRecord]]) -> None:
    header = [
        "CHROM", "POS", "REF", "ALT", "LOW_PARENT_GT", "HIGH_PARENT_GT", "HIGH_GT", "LOW_GT",
        "LOW_PARENT_STATE", "HIGH_PARENT_STATE", "HIGH_PARENT_ALLELE", "HIGH_DP", "LOW_DP",
        "HIGH_PARENT_COUNT", "LOW_PARENT_COUNT", "HIGH_SNP_INDEX", "LOW_SNP_INDEX",
        "DELTA_SNP_INDEX", "ED",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for chrom in sorted(chr_sites, key=natural_chr_key):
            for record in chr_sites[chrom]:
                row = [
                    record.chrom, str(record.pos), record.ref, record.alt, record.parent1_gt, record.parent2_gt,
                    record.high_gt, record.low_gt, str(record.parent1_state), str(record.parent2_state),
                    record.parent2_allele_in_refalt, str(record.high_dp), str(record.low_dp),
                    str(record.high_parent2_count), str(record.low_parent2_count),
                    format_float(record.high_snp_index), format_float(record.low_snp_index),
                    format_float(record.delta_snp_index), format_float(record.ed_value),
                ]
                handle.write("\t".join(row) + "\n")


def write_windows(path: str, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("CHROM\tWINDOW_START\tWINDOW_END\tWINDOW_CENTER\tUSABLE_WINDOW\n")
        return
    columns = list(rows[0].keys())
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            fields: List[str] = []
            for column in columns:
                value = row[column]
                fields.append(format_float(value) if isinstance(value, float) else str(value))
            handle.write("\t".join(fields) + "\n")


def write_plot_table(path: str, rows: Sequence[Dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("#CHROM\tPOS\tTHRESHOLD\tHIGH_SNP_INDEX\tLOW_SNP_INDEX\tDELTA_SNP_INDEX\tSNP_N\n")
        for row in rows:
            if int(row.get("USABLE_WINDOW", 0)) != 1:
                continue
            handle.write(
                "\t".join(
                    [
                        str(row["CHROM"]),
                        str(row["WINDOW_CENTER"]),
                        format_float(float(row["ABS_DELTA_THRESHOLD_99"])),
                        format_float(float(row["HIGH_SNP_INDEX"])),
                        format_float(float(row["LOW_SNP_INDEX"])),
                        format_float(float(row["DELTA_SNP_INDEX"])),
                        str(row["SNP_N"]),
                    ]
                ) + "\n"
            )


def top_windows(rows: Sequence[Dict[str, object]], top_n: int = 20) -> List[Dict[str, object]]:
    usable = [row for row in rows if int(row.get("USABLE_WINDOW", 0)) == 1]
    usable = [row for row in usable if not math.isnan(float(row.get("DELTA_SNP_INDEX", math.nan)))]
    usable.sort(key=lambda row: abs(float(row["DELTA_SNP_INDEX"])), reverse=True)
    return usable[:top_n]


def write_summary(path: str, args: argparse.Namespace, counters: Dict[str, int], rows: Sequence[Dict[str, object]]) -> None:
    top = top_windows(rows, top_n=20)
    usable_windows = sum(1 for row in rows if int(row.get("USABLE_WINDOW", 0)) == 1)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# ploidy=30 BSA summary\n\n")
        handle.write("## Input\n\n")
        handle.write(f"- VCF: `{args.vcf}`\n")
        handle.write(f"- High bulk: `{args.high_bulk}`\n")
        handle.write(f"- Low bulk: `{args.low_bulk}`\n")
        high_parent_name, low_parent_name = resolve_parent_names(args)
        handle.write(f"- High parent: `{high_parent_name}`\n")
        handle.write(f"- Low parent: `{low_parent_name}`\n")
        handle.write(f"- Window size: `{args.window_size}`\n")
        handle.write(f"- Step size: `{args.step_size}`\n")
        handle.write(f"- Min parent DP: `{args.min_parent_dp}`\n")
        handle.write(f"- Min bulk DP: `{args.min_bulk_dp}`\n")
        handle.write(f"- Min window SNPs: `{args.min_window_snps}`\n")
        handle.write("\n## QC counts\n\n")
        for key in ["total_records", "short_rows", "multiallelic_rows", "non_snp_rows", "parent_not_callable", "parent_not_informative", "bulk_low_dp", "usable_sites"]:
            handle.write(f"- {key}: `{counters.get(key, 0)}`\n")

        handle.write("\n## Output interpretation\n\n")
        handle.write("- `HIGH_SNP_INDEX` and `LOW_SNP_INDEX` are both high-parent-oriented allele frequencies.\n")
        handle.write("- `DELTA_SNP_INDEX = HIGH_SNP_INDEX - LOW_SNP_INDEX`.\n")
        handle.write("- `ABS_DELTA_THRESHOLD_99` is the per-window reference threshold used in the plot table.\n")
        handle.write(f"- Usable windows: `{usable_windows}`\n")

        handle.write("\n## Top windows by |delta SNP-index|\n\n")
        if not top:
            handle.write("No usable windows passed the minimum SNP count threshold.\n")
        else:
            handle.write("| CHROM | START | END | SNP_N | HIGH_SNP_INDEX | LOW_SNP_INDEX | DELTA_SNP_INDEX |\n")
            handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
            for row in top:
                handle.write(
                    "| {chrom} | {start} | {end} | {snp_n} | {high} | {low} | {delta} |\n".format(
                        chrom=row["CHROM"],
                        start=row["WINDOW_START"],
                        end=row["WINDOW_END"],
                        snp_n=row["SNP_N"],
                        high=format_float(float(row["HIGH_SNP_INDEX"])),
                        low=format_float(float(row["LOW_SNP_INDEX"])),
                        delta=format_float(float(row["DELTA_SNP_INDEX"])),
                    )
                )

        handle.write("\n## Plot command\n\n")
        handle.write("```bash\n")
        handle.write(
            "python bsa.fig.py "
            f"-f {os.path.basename(args.out_prefix)}.plot.tsv "
            f"-snpindex -o {os.path.basename(args.out_prefix)}.snpindex.pdf\n"
        )
        handle.write("```\n")


def main() -> None:
    args = parse_args()
    if args.window_size <= 0 or args.step_size <= 0:
        raise SystemExit("--window-size and --step-size must be positive integers")
    if args.step_size > args.window_size:
        raise SystemExit("--step-size should usually be <= --window-size")

    run_dir = os.path.join(args.output_root, args.out_prefix)
    os.makedirs(run_dir, exist_ok=True)

    result = process_vcf(args)
    if len(result) == 3:
        chr_sites, counters, chrom_counters = result
    else:
        chr_sites, counters = result  # pragma: no cover
        chrom_counters = defaultdict(Counter)
    rows = build_windows(chr_sites, args)

    sites_path = os.path.join(run_dir, f"{args.out_prefix}.sites.tsv")
    windows_path = os.path.join(run_dir, f"{args.out_prefix}.windows.tsv")
    plot_path = os.path.join(run_dir, f"{args.out_prefix}.plot.tsv")
    summary_path = os.path.join(run_dir, f"{args.out_prefix}.summary.md")

    write_sites(sites_path, chr_sites)
    write_windows(windows_path, rows)
    write_plot_table(plot_path, rows)
    write_summary(summary_path, args, counters, rows)

    if args.diagnostic_report:
        chroms = ["Chr10", "Chr11", "Chr12"]
        diag_tsv = args.diagnostic_report + ".tsv"
        diag_md = args.diagnostic_report + ".md"
        write_diagnostic_report(diag_tsv, chrom_counters, chroms)
        with open(diag_md, "w", encoding="utf-8") as out:
            out.write("# Chr10-12 diagnostic report\n\n")
            out.write(f"- VCF: `{args.vcf}`\n")
            out.write(f"- High bulk: `{args.high_bulk}`\n")
            out.write(f"- Low bulk: `{args.low_bulk}`\n")
            high_parent_name, low_parent_name = resolve_parent_names(args)
            out.write(f"- High parent: `{high_parent_name}`\n")
            out.write(f"- Low parent: `{low_parent_name}`\n")
            out.write(f"- Min parent DP: `{args.min_parent_dp}`\n")
            out.write(f"- Min bulk DP: `{args.min_bulk_dp}`\n")
            out.write(f"- Min window SNPs: `{args.min_window_snps}`\n\n")
            out.write("## Filter summary\n\n")
            out.write("```text\n")
            with open(diag_tsv, "r", encoding="utf-8") as diag_in:
                for line in diag_in:
                    out.write(line)
            out.write("```\n")
        print(f"[OK] diagnostic report: {diag_tsv}")
        print(f"[OK] diagnostic markdown: {diag_md}")

    print(f"[OK] site table: {sites_path}")
    print(f"[OK] window table: {windows_path}")
    print(f"[OK] plot table: {plot_path}")
    print(f"[OK] summary: {summary_path}")


if __name__ == "__main__":
    main()
