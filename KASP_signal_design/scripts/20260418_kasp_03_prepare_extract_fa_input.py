#!/usr/bin/env python3
# 生成日期: 2026-04-18
# 脚本作用: 将通过过滤的候选位点整理成 extract_fa.py 需要的四列表, 并输出推荐标记清单。
# 输入文件: 02_pass_sites.tsv 和 02_checked_sites.tsv。
# 输出文件: 03_extract_fa_input.tsv 和 03_recommended_markers.tsv。
# 原理与逻辑:
# 1. 左右两侧分别取最终通过的前 N 个标记。
# 2. extract_fa 输入文件严格保持四列无表头, 避免和 extract_fa.py 的接口冲突。
# 3. 推荐清单同时写出已检查候选数、通过数和不足原因, 方便判断是否需要放宽策略。
# 使用方法:
# python 20260418_kasp_03_prepare_extract_fa_input.py --pass-sites 02_pass_sites.tsv --checked-sites 02_checked_sites.tsv --output-extract-input 03_extract_fa_input.tsv --output-recommended 03_recommended_markers.tsv

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Sequence


RECOMMENDED_FIELDS = [
    "side",
    "selected_count",
    "target_count",
    "checked_candidates",
    "pass_count",
    "shortage_reason",
    "chrom",
    "pos",
    "ref",
    "alt",
    "delta_snp_index",
    "checked_rank",
    "parent1_gt",
    "parent2_gt",
    "high_pool_gt",
    "low_pool_gt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="整理 extract_fa.py 输入并输出推荐标记清单")
    parser.add_argument("--pass-sites", required=True, help="02_pass_sites.tsv 路径")
    parser.add_argument("--checked-sites", required=True, help="02_checked_sites.tsv 路径")
    parser.add_argument("--output-extract-input", required=True, help="extract_fa.py 输入文件路径")
    parser.add_argument("--output-recommended", required=True, help="推荐标记清单输出路径")
    parser.add_argument("--target-markers-per-side", type=int, default=2, help="每侧目标标记数量, 默认 2")
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader]


def sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: int(row["checked_rank"]))


def summarize_shortage(side_checked: list[dict[str, str]], side_pass: list[dict[str, str]], target_count: int) -> str:
    if len(side_pass) >= target_count:
        return "quota_met"
    if not side_checked:
        return "no_candidates_checked"
    reasons = Counter(
        row["filter_reason"] for row in side_checked if row["filter_status"] != "PASS" and row["filter_reason"]
    )
    if not reasons:
        return "pass_sites_insufficient"
    parts = [f"{reason}={count}" for reason, count in reasons.most_common()]
    return "pass_sites_insufficient_after_filters:" + ",".join(parts)


def write_extract_input(path: Path, selected_rows: Sequence[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in selected_rows:
            handle.write(f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\n")


def write_recommended(
    path: Path,
    checked_rows: list[dict[str, str]],
    pass_rows: list[dict[str, str]],
    target_count: int,
) -> list[dict[str, str]]:
    grouped_checked: dict[str, list[dict[str, str]]] = {}
    grouped_pass: dict[str, list[dict[str, str]]] = {}
    side_order: list[str] = []
    for row in checked_rows:
        side = row["side"]
        if side not in grouped_checked:
            grouped_checked[side] = []
            side_order.append(side)
        grouped_checked[side].append(row)
    for row in pass_rows:
        side = row["side"]
        if side not in grouped_pass:
            grouped_pass[side] = []
        if side not in side_order:
            side_order.append(side)
        grouped_pass[side].append(row)

    selected_rows: list[dict[str, str]] = []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECOMMENDED_FIELDS, delimiter="\t")
        writer.writeheader()
        for side in side_order:
            side_checked = sort_rows(grouped_checked.get(side, []))
            side_pass = sort_rows(grouped_pass.get(side, []))[:target_count]
            shortage_reason = summarize_shortage(side_checked, side_pass, target_count)
            if side_pass:
                for row in side_pass:
                    output_row = {
                        "side": side,
                        "selected_count": len(side_pass),
                        "target_count": target_count,
                        "checked_candidates": len(side_checked),
                        "pass_count": len(grouped_pass.get(side, [])),
                        "shortage_reason": shortage_reason,
                        "chrom": row["chrom"],
                        "pos": row["pos"],
                        "ref": row["ref"],
                        "alt": row["alt"],
                        "delta_snp_index": row["delta_snp_index"],
                        "checked_rank": row["checked_rank"],
                        "parent1_gt": row["parent1_gt"],
                        "parent2_gt": row["parent2_gt"],
                        "high_pool_gt": row["high_pool_gt"],
                        "low_pool_gt": row["low_pool_gt"],
                    }
                    writer.writerow(output_row)
                    selected_rows.append(row)
            else:
                writer.writerow(
                    {
                        "side": side,
                        "selected_count": 0,
                        "target_count": target_count,
                        "checked_candidates": len(side_checked),
                        "pass_count": len(grouped_pass.get(side, [])),
                        "shortage_reason": shortage_reason,
                        "chrom": "",
                        "pos": "",
                        "ref": "",
                        "alt": "",
                        "delta_snp_index": "",
                        "checked_rank": "",
                        "parent1_gt": "",
                        "parent2_gt": "",
                        "high_pool_gt": "",
                        "low_pool_gt": "",
                    }
                )
    return selected_rows


def main() -> int:
    args = parse_args()
    pass_sites = Path(args.pass_sites)
    checked_sites = Path(args.checked_sites)
    if not pass_sites.is_file():
        raise FileNotFoundError(f"通过位点文件不存在: {pass_sites}")
    if not checked_sites.is_file():
        raise FileNotFoundError(f"检查明细文件不存在: {checked_sites}")

    pass_rows = read_table(pass_sites)
    checked_rows = read_table(checked_sites)
    output_extract_input = Path(args.output_extract_input)
    output_recommended = Path(args.output_recommended)
    output_extract_input.parent.mkdir(parents=True, exist_ok=True)
    output_recommended.parent.mkdir(parents=True, exist_ok=True)
    selected_rows = write_recommended(output_recommended, checked_rows, pass_rows, args.target_markers_per_side)
    write_extract_input(output_extract_input, selected_rows)
    print(f"extract_fa 输入已输出: {output_extract_input}")
    print(f"推荐标记清单已输出: {output_recommended}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
