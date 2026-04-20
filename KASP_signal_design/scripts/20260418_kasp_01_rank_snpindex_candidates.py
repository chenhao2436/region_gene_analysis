#!/usr/bin/env python3
# 生成日期: 2026-04-18
# 脚本作用: 从 snpindex 文件中提取目标峰区左右两侧候选位点, 初始使用向内收的 1 Mb 区间, 不足时按 1 Mb 一档向峰区中心扩宽。
# 输入文件: snpindex 表格文件, 必须至少包含 CHROM/POS/DELTA_SNP_INDEX 列; 无表头时默认按第 1/2/6 列处理。
# 输出文件: 01_ranked_candidates.tsv 和 01_rank_summary.tsv。
# 原理与逻辑:
# 1. 将类似 Chr05:248.7-255.3MB 的区间统一换算成 bp。
# 2. 初始拆分成左侧 [start, start+1Mb] 与右侧 [end-1Mb, end] 两个区间。
# 3. 若某侧候选不足, 则该侧按 1 Mb 一档继续向峰区中心扩宽, 直到覆盖该侧半区。
# 4. 候选排序优先级为: 扩宽层级低的优先, 同层内 DELTA_SNP_INDEX 高的优先。
# 5. 结果输出为结构化表格, 供后续 VCF 过滤脚本继续使用。
# 使用方法:
# python 20260418_kasp_01_rank_snpindex_candidates.py --region Chr05:248.7-255.3MB --snpindex input.tsv --output-candidates 01_ranked_candidates.tsv --output-summary 01_rank_summary.tsv

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import Iterable, Sequence


ONE_MB = 1_000_000
SUMMARY_FIELDS = [
    "side",
    "chrom",
    "initial_window_start",
    "initial_window_end",
    "max_window_start",
    "max_window_end",
    "total_candidates",
    "initial_rank_limit",
    "initial_slice_count",
    "expansion_levels",
    "max_delta_snp_index",
    "min_delta_snp_index",
    "source_snpindex",
]
CANDIDATE_FIELDS = [
    "side",
    "chrom",
    "pos",
    "delta_snp_index",
    "rank",
    "window_start",
    "window_end",
    "expansion_level",
    "source_snpindex",
]


@dataclass(frozen=True)
class Window:
    side: str
    chrom: str
    start: int
    end: int
    initial_start: int
    initial_end: int


@dataclass(frozen=True)
class Candidate:
    side: str
    chrom: str
    pos: int
    delta_snp_index: float
    rank: int
    window_start: int
    window_end: int
    expansion_level: int
    source_snpindex: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="提取左右两个向内收 1 Mb 区间的高 DELTA_SNP_INDEX 候选位点")
    parser.add_argument("--region", required=True, help="输入峰区, 例如 Chr05:248.7-255.3MB")
    parser.add_argument("--snpindex", required=True, help="snpindex 文件路径")
    parser.add_argument("--output-candidates", required=True, help="候选位点输出路径")
    parser.add_argument("--output-summary", required=True, help="汇总表输出路径")
    parser.add_argument("--initial-rank-limit", type=int, default=20, help="前筛候选数, 默认 20")
    parser.add_argument("--expansion-step-bp", type=int, default=ONE_MB, help="每次向峰区中心扩宽的步长, 默认 1000000 bp")
    return parser.parse_args()


def parse_region(region_text: str) -> tuple[str, int, int]:
    pattern = re.compile(
        r"^\s*(?P<chrom>[^:]+)\s*:\s*(?P<start>\d+(?:\.\d+)?)\s*-\s*(?P<end>\d+(?:\.\d+)?)\s*(?P<unit>MB|M|BP)?\s*$",
        re.IGNORECASE,
    )
    match = pattern.match(region_text)
    if not match:
        raise ValueError(f"无法解析区间: {region_text}")

    chrom = match.group("chrom").strip()
    start_value = float(match.group("start"))
    end_value = float(match.group("end"))
    unit = (match.group("unit") or "MB").upper()
    multiplier = ONE_MB if unit in {"MB", "M"} else 1
    start_bp = int(round(start_value * multiplier))
    end_bp = int(round(end_value * multiplier))
    if start_bp >= end_bp:
        raise ValueError(f"区间起点必须小于终点: {region_text}")
    if end_bp - start_bp < 2 * ONE_MB:
        raise ValueError("输入区间长度必须至少为 2 Mb, 才能拆分出左右各 1 Mb 的向内区间")
    return chrom, start_bp, end_bp


def build_windows(chrom: str, start_bp: int, end_bp: int) -> list[Window]:
    midpoint = start_bp + ((end_bp - start_bp) // 2)
    left_initial_end = min(start_bp + ONE_MB, midpoint)
    right_initial_start = max(end_bp - ONE_MB, midpoint + 1)
    return [
        Window(
            side="left",
            chrom=chrom,
            start=start_bp,
            end=midpoint,
            initial_start=start_bp,
            initial_end=left_initial_end,
        ),
        Window(
            side="right",
            chrom=chrom,
            start=midpoint + 1,
            end=end_bp,
            initial_start=right_initial_start,
            initial_end=end_bp,
        ),
    ]


def normalize_header_name(value: str) -> str:
    return value.strip().lstrip("#").upper().replace("-", "_")


def detect_indices(first_row: Sequence[str]) -> tuple[bool, int, int, int]:
    normalized = [normalize_header_name(cell) for cell in first_row]
    if "CHROM" in normalized and "POS" in normalized:
        for delta_name in ("DELTA_SNP_INDEX", "DELTA_SNPINDEX"):
            if delta_name in normalized:
                return True, normalized.index("CHROM"), normalized.index("POS"), normalized.index(delta_name)
        raise ValueError("表头缺少 DELTA_SNP_INDEX 列")
    if len(first_row) < 6:
        raise ValueError("无表头 snpindex 文件列数不足, 无法按默认第 1/2/6 列解析")
    return False, 0, 1, 5


def iter_rows(path: Path) -> Iterable[list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row and any(cell.strip() for cell in row):
                yield row


def calculate_expansion_level(window: Window, pos: int, expansion_step_bp: int) -> int:
    if window.side == "left":
        if pos <= window.initial_end:
            return 0
        extra = pos - window.initial_end
    else:
        if pos >= window.initial_start:
            return 0
        extra = window.initial_start - pos
    return (extra - 1) // expansion_step_bp + 1


def expansion_window_end(window: Window, expansion_level: int, expansion_step_bp: int) -> tuple[int, int]:
    if window.side == "left":
        return window.start, min(window.end, window.initial_end + expansion_level * expansion_step_bp)
    return max(window.start, window.initial_start - expansion_level * expansion_step_bp), window.end


def collect_candidates(snpindex_path: Path, windows: Sequence[Window], expansion_step_bp: int) -> list[Candidate]:
    rows = iter_rows(snpindex_path)
    try:
        first_row = next(rows)
    except StopIteration as exc:
        raise ValueError("snpindex 文件为空") from exc

    has_header, chrom_idx, pos_idx, delta_idx = detect_indices(first_row)
    data_rows = rows if has_header else chain([first_row], rows)

    by_side: dict[str, list[tuple[int, float, int]]] = {window.side: [] for window in windows}
    left_window = next(item for item in windows if item.side == "left")
    right_window = next(item for item in windows if item.side == "right")
    target_chrom = left_window.chrom
    midpoint = left_window.end

    for row in data_rows:
        if max(chrom_idx, pos_idx, delta_idx) >= len(row):
            continue
        chrom = row[chrom_idx].strip()
        if chrom != target_chrom:
            continue
        try:
            pos = int(row[pos_idx].strip())
            delta = float(row[delta_idx].strip())
        except ValueError:
            continue
        if math.isnan(delta) or delta <= 0:
            continue
        if pos <= midpoint and left_window.start <= pos <= left_window.end:
            by_side["left"].append((pos, delta, calculate_expansion_level(left_window, pos, expansion_step_bp)))
        elif pos >= midpoint + 1 and right_window.start <= pos <= right_window.end:
            by_side["right"].append((pos, delta, calculate_expansion_level(right_window, pos, expansion_step_bp)))

    candidates: list[Candidate] = []
    for side in ("left", "right"):
        window = next(item for item in windows if item.side == side)
        sorted_rows = sorted(by_side[side], key=lambda item: (item[2], -item[1], item[0]))
        for index, (pos, delta, expansion_level) in enumerate(sorted_rows, start=1):
            current_start, current_end = expansion_window_end(window, expansion_level, expansion_step_bp)
            candidates.append(
                Candidate(
                    side=side,
                    chrom=window.chrom,
                    pos=pos,
                    delta_snp_index=delta,
                    rank=index,
                    window_start=current_start,
                    window_end=current_end,
                    expansion_level=expansion_level,
                    source_snpindex=str(snpindex_path),
                )
            )
    return candidates


def write_candidates(path: Path, candidates: Sequence[Candidate]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_FIELDS, delimiter="\t")
        writer.writeheader()
        for item in candidates:
            writer.writerow(
                {
                    "side": item.side,
                    "chrom": item.chrom,
                    "pos": item.pos,
                    "delta_snp_index": f"{item.delta_snp_index:.6f}",
                    "rank": item.rank,
                    "window_start": item.window_start,
                    "window_end": item.window_end,
                    "expansion_level": item.expansion_level,
                    "source_snpindex": item.source_snpindex,
                }
            )


def write_summary(path: Path, windows: Sequence[Window], candidates: Sequence[Candidate], initial_rank_limit: int, source_path: Path) -> None:
    grouped: dict[str, list[Candidate]] = {window.side: [] for window in windows}
    for item in candidates:
        grouped[item.side].append(item)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for window in windows:
            side_candidates = grouped[window.side]
            writer.writerow(
                {
                    "side": window.side,
                    "chrom": window.chrom,
                    "initial_window_start": window.initial_start,
                    "initial_window_end": window.initial_end,
                    "max_window_start": window.start,
                    "max_window_end": window.end,
                    "total_candidates": len(side_candidates),
                    "initial_rank_limit": initial_rank_limit,
                    "initial_slice_count": min(initial_rank_limit, len(side_candidates)),
                    "expansion_levels": side_candidates[-1].expansion_level + 1 if side_candidates else 0,
                    "max_delta_snp_index": f"{side_candidates[0].delta_snp_index:.6f}" if side_candidates else "",
                    "min_delta_snp_index": f"{side_candidates[-1].delta_snp_index:.6f}" if side_candidates else "",
                    "source_snpindex": str(source_path),
                }
            )


def main() -> int:
    args = parse_args()
    snpindex_path = Path(args.snpindex)
    if not snpindex_path.is_file():
        raise FileNotFoundError(f"snpindex 文件不存在: {snpindex_path}")

    chrom, start_bp, end_bp = parse_region(args.region)
    windows = build_windows(chrom, start_bp, end_bp)
    candidates = collect_candidates(snpindex_path, windows, args.expansion_step_bp)

    output_candidates = Path(args.output_candidates)
    output_summary = Path(args.output_summary)
    output_candidates.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    write_candidates(output_candidates, candidates)
    write_summary(output_summary, windows, candidates, args.initial_rank_limit, snpindex_path)
    print(f"候选位点已输出: {output_candidates}")
    print(f"汇总表已输出: {output_summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
