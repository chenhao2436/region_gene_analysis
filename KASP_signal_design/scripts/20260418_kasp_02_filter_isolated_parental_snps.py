#!/usr/bin/env python3
# 生成日期: 2026-04-18
# 脚本作用: 逐个检查高优先级候选位点, 保留双亲纯合且差异、双等位 SNP、前后 50 bp 无其他变异的孤立位点。
# 输入文件: 01_ranked_candidates.tsv, bgzip+tabix 索引的 VCF 文件。
# 输出文件: 02_checked_sites.tsv 和 02_pass_sites.tsv。
# 原理与逻辑:
# 1. 候选位点按左右侧分别依 DELTA_SNP_INDEX 降序检查。
# 2. 目标位点必须是双等位 SNP, 且两个亲本都纯合、基因型不同。
# 3. 再查询目标位点前后 50 bp 范围, 如果存在其他变异位点则剔除。
# 4. 每侧优先检查前 20 个, 不足目标数时继续向后补位, 直到够用或候选耗尽。
# 使用方法:
# python 20260418_kasp_02_filter_isolated_parental_snps.py --candidates 01_ranked_candidates.tsv --vcf input.vcf.gz --parent1 C3_Parent --parent2 C6_Parent --high-pool HIGH --low-pool LOW --output-checked 02_checked_sites.tsv --output-pass 02_pass_sites.tsv

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


CHECKED_FIELDS = [
    "side",
    "chrom",
    "pos",
    "ref",
    "alt",
    "delta_snp_index",
    "parent1_gt",
    "parent2_gt",
    "high_pool_gt",
    "low_pool_gt",
    "filter_status",
    "filter_reason",
    "checked_rank",
]


@dataclass(frozen=True)
class Candidate:
    side: str
    chrom: str
    pos: int
    delta_snp_index: float
    rank: int


@dataclass(frozen=True)
class SiteRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    sample_gts: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="过滤亲本纯合差异且前后无其他变异的 KASP 候选位点")
    parser.add_argument("--candidates", required=True, help="01_ranked_candidates.tsv 路径")
    parser.add_argument("--vcf", required=True, help="bgzip 压缩并建立索引的 VCF 文件")
    parser.add_argument("--parent1", required=True, help="亲本 1 样本名")
    parser.add_argument("--parent2", required=True, help="亲本 2 样本名")
    parser.add_argument("--high-pool", required=True, help="高池样本名")
    parser.add_argument("--low-pool", required=True, help="低池样本名")
    parser.add_argument("--flank-bp", type=int, default=50, help="两侧不能有其他变异的范围, 默认 50")
    parser.add_argument("--target-markers-per-side", type=int, default=2, help="每侧目标标记数量, 默认 2")
    parser.add_argument("--initial-rank-limit", type=int, default=20, help="优先尝试的前筛数量, 默认 20")
    parser.add_argument("--output-checked", required=True, help="检查明细输出路径")
    parser.add_argument("--output-pass", required=True, help="通过位点输出路径")
    return parser.parse_args()


def run_command(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=False, capture_output=True, text=True, encoding="utf-8")


def load_candidates(path: Path) -> dict[str, list[Candidate]]:
    grouped: dict[str, list[Candidate]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"side", "chrom", "pos", "delta_snp_index", "rank"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("候选表缺少必要列: side/chrom/pos/delta_snp_index/rank")
        for row in reader:
            side = row["side"].strip()
            grouped.setdefault(side, [])
            grouped[side].append(
                Candidate(
                    side=side,
                    chrom=row["chrom"].strip(),
                    pos=int(row["pos"]),
                    delta_snp_index=float(row["delta_snp_index"]),
                    rank=int(row["rank"]),
                )
            )
    for side in grouped:
        grouped[side].sort(key=lambda item: item.rank)
    return grouped


def parse_gt(raw_gt: str) -> tuple[bool, str]:
    gt = raw_gt.split(":", 1)[0].strip()
    if gt in {".", "./.", ".|."}:
        return False, gt
    normalized = gt.replace("|", "/")
    alleles = normalized.split("/")
    if len(alleles) != 2 or "." in alleles:
        return False, gt
    return alleles[0] == alleles[1], gt


def query_region(vcf_path: Path, region: str, sample_names: Sequence[str]) -> list[SiteRecord]:
    command = [
        "bcftools",
        "query",
        "-r",
        region,
        "-s",
        ",".join(sample_names),
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
        str(vcf_path),
    ]
    result = run_command(command)
    if result.returncode != 0:
        stderr = result.stderr.strip() or "未知错误"
        raise RuntimeError(f"bcftools query 失败: {stderr}")

    records: list[SiteRecord] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4 + len(sample_names):
            continue
        records.append(
            SiteRecord(
                chrom=parts[0],
                pos=int(parts[1]),
                ref=parts[2],
                alt=parts[3],
                sample_gts=parts[4 : 4 + len(sample_names)],
            )
        )
    return records


def evaluate_candidate(candidate: Candidate, vcf_path: Path, sample_names: Sequence[str], flank_bp: int) -> dict[str, str]:
    region = f"{candidate.chrom}:{candidate.pos}-{candidate.pos}"
    default_row = {
        "side": candidate.side,
        "chrom": candidate.chrom,
        "pos": str(candidate.pos),
        "ref": "",
        "alt": "",
        "delta_snp_index": f"{candidate.delta_snp_index:.6f}",
        "parent1_gt": "",
        "parent2_gt": "",
        "high_pool_gt": "",
        "low_pool_gt": "",
        "filter_status": "FAIL",
        "filter_reason": "",
        "checked_rank": str(candidate.rank),
    }

    try:
        target_records = query_region(vcf_path, region, sample_names)
    except RuntimeError as exc:
        default_row["filter_reason"] = str(exc)
        return default_row

    if not target_records:
        default_row["filter_reason"] = "target_site_not_found"
        return default_row
    if len(target_records) != 1:
        default_row["filter_reason"] = "multiple_records_at_target_site"
        return default_row

    target = target_records[0]
    default_row["ref"] = target.ref
    default_row["alt"] = target.alt
    default_row["parent1_gt"] = target.sample_gts[0]
    default_row["parent2_gt"] = target.sample_gts[1]
    default_row["high_pool_gt"] = target.sample_gts[2]
    default_row["low_pool_gt"] = target.sample_gts[3]

    if "," in target.alt:
        default_row["filter_reason"] = "multiallelic_site"
        return default_row
    if len(target.ref) != 1 or len(target.alt) != 1:
        default_row["filter_reason"] = "indel_or_non_snp_site"
        return default_row

    parent1_hom, _ = parse_gt(target.sample_gts[0])
    parent2_hom, _ = parse_gt(target.sample_gts[1])
    if not parent1_hom:
        default_row["filter_reason"] = "parent1_not_homozygous"
        return default_row
    if not parent2_hom:
        default_row["filter_reason"] = "parent2_not_homozygous"
        return default_row

    parent1_gt = target.sample_gts[0].split(":", 1)[0].replace("|", "/")
    parent2_gt = target.sample_gts[1].split(":", 1)[0].replace("|", "/")
    if parent1_gt == parent2_gt:
        default_row["filter_reason"] = "parents_share_same_genotype"
        return default_row

    flank_start = max(1, candidate.pos - flank_bp)
    flank_end = candidate.pos + flank_bp
    flank_region = f"{candidate.chrom}:{flank_start}-{flank_end}"
    try:
        flank_records = query_region(vcf_path, flank_region, sample_names)
    except RuntimeError as exc:
        default_row["filter_reason"] = str(exc)
        return default_row

    nearby_positions = {record.pos for record in flank_records if record.pos != candidate.pos}
    if nearby_positions:
        default_row["filter_reason"] = "nearby_variant_present"
        return default_row

    default_row["filter_status"] = "PASS"
    default_row["filter_reason"] = "pass"
    return default_row


def filter_candidates(
    grouped_candidates: dict[str, list[Candidate]],
    vcf_path: Path,
    sample_names: Sequence[str],
    flank_bp: int,
    target_markers_per_side: int,
    initial_rank_limit: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    checked_rows: list[dict[str, str]] = []
    pass_rows: list[dict[str, str]] = []

    for side in grouped_candidates:
        side_candidates = grouped_candidates.get(side, [])
        side_pass_count = 0
        for index, candidate in enumerate(side_candidates, start=1):
            if side_pass_count >= target_markers_per_side and index > initial_rank_limit:
                break
            result_row = evaluate_candidate(candidate, vcf_path, sample_names, flank_bp)
            checked_rows.append(result_row)
            if result_row["filter_status"] == "PASS":
                pass_rows.append(result_row)
                side_pass_count += 1
                if side_pass_count >= target_markers_per_side:
                    break
    return checked_rows, pass_rows


def write_rows(path: Path, rows: Sequence[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHECKED_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    vcf_path = Path(args.vcf)
    candidates_path = Path(args.candidates)
    if not vcf_path.is_file():
        raise FileNotFoundError(f"VCF 文件不存在: {vcf_path}")
    if not candidates_path.is_file():
        raise FileNotFoundError(f"候选表不存在: {candidates_path}")

    sample_names = [args.parent1, args.parent2, args.high_pool, args.low_pool]
    grouped_candidates = load_candidates(candidates_path)
    checked_rows, pass_rows = filter_candidates(
        grouped_candidates=grouped_candidates,
        vcf_path=vcf_path,
        sample_names=sample_names,
        flank_bp=args.flank_bp,
        target_markers_per_side=args.target_markers_per_side,
        initial_rank_limit=args.initial_rank_limit,
    )

    output_checked = Path(args.output_checked)
    output_pass = Path(args.output_pass)
    output_checked.parent.mkdir(parents=True, exist_ok=True)
    output_pass.parent.mkdir(parents=True, exist_ok=True)
    write_rows(output_checked, checked_rows)
    write_rows(output_pass, pass_rows)
    print(f"检查明细已输出: {output_checked}")
    print(f"通过位点已输出: {output_pass}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
