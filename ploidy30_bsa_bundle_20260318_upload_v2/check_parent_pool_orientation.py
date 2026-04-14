#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
验证 BSA 亲本与混池等位基因同向性的诊断脚本
用途：从 sites.tsv 提取极端信号位点，并从原始 VCF 追踪 AD（Allele Depth）值，
      通过明确的公式展示 C6、C3 与 H、L 混池的等位基因频率对应关系。
"""

import os
import sys
import random
import subprocess
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract and verify parent-pool orientation.")
    parser.add_argument("--sites", default="/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000/hunchi_and_qinben.prepared_single_sample.genome.joint_184_3000000.sites.tsv")
    parser.add_argument("--vcf", default="/data2/chenh/bsa/CleanData/vcf/gvcf2vcf_genome_gpt_3m6d/joint_vcf_bak/hunchi_and_qinben.prepared_single_sample.genome.joint.vcf.gz")
    parser.add_argument("--outdir", default="/data2/chenh/bsa/CleanData/bsa_HunchiData/BSA_out/output/qinebn_hunchi_check")
    parser.add_argument("--sample-c6", default="C6_Parent", help="Sample name for high parent")
    parser.add_argument("--sample-c3", default="C3_Parent", help="Sample name for low parent")
    parser.add_argument("--sample-h", default="F9-184H", help="Sample name for high bulk")
    parser.add_argument("--sample-l", default="F9-184L", help="Sample name for low bulk")
    return parser.parse_args()

def parse_ad(ad_str):
    if ad_str in ('.', './.', '.|.'): return 0, 0
    parts = ad_str.split(',')
    if len(parts) >= 2:
        try:
            return int(parts[0]), int(parts[1])
        except ValueError:
            return 0, 0
    return 0, 0

def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error running command: {cmd}\n{result.stderr}")
        sys.exit(1)
    return result.stdout.strip()

def main():
    args = parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Read sites.tsv and filter
    print("Parsing sites.tsv...")
    chr07_cands = []
    chr12_cands = []
    
    try:
        with open(args.sites, 'r', encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            try:
                idx_chrom = header.index('CHROM')
                idx_pos = header.index('POS')
                idx_delta = header.index('DELTA_SNP_INDEX')
                idx_h_snp = header.index('HIGH_SNP_INDEX')
                idx_l_snp = header.index('LOW_SNP_INDEX')
            except ValueError as e:
                print(f"Missing column in sites.tsv: {e}")
                sys.exit(1)
                
            for line in f:
                if not line.strip(): continue
                parts = line.strip().split('\t')
                max_idx = max(idx_delta, idx_h_snp, idx_l_snp)
                if len(parts) <= max_idx: continue
                chrom, pos, delta_str = parts[idx_chrom], parts[idx_pos], parts[idx_delta]
                
                try:
                    delta = float(delta_str)
                except ValueError:
                    continue
                    
                if chrom == 'Chr07' and delta < -0.4:
                    chr07_cands.append((chrom, pos, delta_str, parts[idx_h_snp], parts[idx_l_snp]))
                elif chrom == 'Chr12' and delta > 0.38:
                    chr12_cands.append((chrom, pos, delta_str, parts[idx_h_snp], parts[idx_l_snp]))
    except FileNotFoundError:
        print(f"File not found: {args.sites}")
        sys.exit(1)
        
    print(f"Found {len(chr07_cands)} candidates for Chr07, {len(chr12_cands)} candidates for Chr12.")
    
    random.seed(42) # For reproducibility
    chr07_selected = random.sample(chr07_cands, min(10, len(chr07_cands)))
    chr12_selected = random.sample(chr12_cands, min(10, len(chr12_cands)))
    
    selected_sites = chr07_selected + chr12_selected
    if not selected_sites:
        print("No sites met the filtering criteria.")
        sys.exit(0)
        
    # Sort logically
    selected_sites.sort(key=lambda x: (x[0], int(x[1])))
    
    # Extract positions for bcftools
    regions = [f"{s[0]}:{s[1]}-{s[1]}" for s in selected_sites]
    regions_str = ",".join(regions)
    
    # 2. Extract specific samples and fields from VCF using bcftools
    print(f"Extracting VCF records for {len(selected_sites)} sites...")
    cmd_header = f"bcftools view -h {args.vcf} | grep '^#CHROM'"
    vcf_header_line = run_cmd(cmd_header).split('\t')
    
    try:
        idx_c6 = vcf_header_line.index(args.sample_c6)
        idx_c3 = vcf_header_line.index(args.sample_c3)
        idx_h = vcf_header_line.index(args.sample_h)
        idx_l = vcf_header_line.index(args.sample_l)
    except ValueError as e:
        print(f"Sample not found in VCF header: {e}")
        print("Available samples:", vcf_header_line[9:])
        sys.exit(1)
        
    cmd_vcf = f"bcftools view -H -r {regions_str} {args.vcf}"
    vcf_data = run_cmd(cmd_vcf).split('\n')
    
    vcf_dict = {}
    for line in vcf_data:
        if not line.strip(): continue
        parts = line.split('\t')
        chrom, pos = parts[0], parts[1]
        vcf_dict[(chrom, pos)] = parts
        
    # 3. Process and output tables
    out_headers = [
        "CHROM", "POS", "REF>ALT", 
        f"C6({args.sample_c6})_AD", f"C3({args.sample_c3})_AD", 
        f"H({args.sample_h})_AD", f"L({args.sample_l})_AD",
        "C6_STATE", "C3_STATE", "HP_ALLELE(C6)",
        "H_SNP_INDEX", "L_SNP_INDEX", "DELTA_SNP_INDEX",
        "SITES_HIGH_SNP", "SITES_LOW_SNP", "SITES_DELTA"
    ]
    
    def process_chrom(chrom_target, out_file, sites):
        print(f"Writing {out_file}...")
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write("# 列说明注释：\n")
            f.write("# CHROM, POS: 染色体和物理位置\n")
            f.write("# REF>ALT: 参考基因组与替代等位基因\n")
            f.write(f"# C6({args.sample_c6})_AD: 亲本C6(高表型)的 REF,ALT reads 深度\n")
            f.write(f"# C3({args.sample_c3})_AD: 亲本C3(低表型)的 REF,ALT reads 深度\n")
            f.write(f"# H({args.sample_h})_AD: 混池H(高表型)的 REF,ALT reads 深度\n")
            f.write(f"# L({args.sample_l})_AD: 混池L(低表型)的 REF,ALT reads 深度\n")
            f.write("# C6_STATE / C3_STATE: 亲本纯合状态 (0=REF纯合, 1=ALT纯合)\n")
            f.write("# HP_ALLELE: 高表型亲本(C6)含有的等位基因方向 (REF或ALT)\n")
            f.write("# H_SNP_INDEX, L_SNP_INDEX: 混池中C6等位基因的频率 (C6_allele_count / Total_DP)\n")
            f.write("# DELTA_SNP_INDEX: H_SNP_INDEX - L_SNP_INDEX\n")
            f.write("# SITES_...: sites.tsv 中记录的原始计算值（用于精确比对）\n")
            f.write("\t".join(out_headers) + "\n")
            
            for site in sites:
                chrom, pos, sites_delta, sites_h, sites_l = site
                if chrom != chrom_target: continue
                
                parts = vcf_dict.get((chrom, pos))
                if not parts:
                    f.write(f"{chrom}\t{pos}\tN/A\t(Not Found in VCF...)\n")
                    continue
                    
                ref, alt = parts[3], parts[4]
                fmt = parts[8].split(':')
                
                def get_ad(sample_parts):
                    if sample_parts in ('.', './.', '.|.'): return "0,0"
                    sp = sample_parts.split(':')
                    if 'AD' in fmt:
                        return sp[fmt.index('AD')]
                    return "0,0"
                    
                c6_ad_str = get_ad(parts[idx_c6])
                c3_ad_str = get_ad(parts[idx_c3])
                h_ad_str  = get_ad(parts[idx_h])
                l_ad_str  = get_ad(parts[idx_l])
                
                c6_ref, c6_alt = parse_ad(c6_ad_str)
                c3_ref, c3_alt = parse_ad(c3_ad_str)
                
                c6_state = "0" if c6_ref > c6_alt else "1" if c6_alt > c6_ref else "NA"
                c3_state = "0" if c3_ref > c3_alt else "1" if c3_alt > c3_ref else "NA"
                
                hp_allele = "REF" if c6_state == "0" else "ALT" if c6_state == "1" else "NA"
                
                h_ref, h_alt = parse_ad(h_ad_str)
                l_ref, l_alt = parse_ad(l_ad_str)
                
                h_dp = h_ref + h_alt
                l_dp = l_ref + l_alt
                
                if hp_allele == "REF":
                    h_idx = h_ref / h_dp if h_dp > 0 else 0.0
                    l_idx = l_ref / l_dp if l_dp > 0 else 0.0
                elif hp_allele == "ALT":
                    h_idx = h_alt / h_dp if h_dp > 0 else 0.0
                    l_idx = l_alt / l_dp if l_dp > 0 else 0.0
                else:
                    h_idx = 0.0
                    l_idx = 0.0
                    
                calc_delta = h_idx - l_idx
                
                row = [
                    chrom, pos, f"{ref}>{alt}",
                    c6_ad_str, c3_ad_str, h_ad_str, l_ad_str,
                    c6_state, c3_state, hp_allele,
                    f"{h_idx:.4f}", f"{l_idx:.4f}", f"{calc_delta:.4f}",
                    sites_h, sites_l, sites_delta
                ]
                f.write("\t".join(row) + "\n")

    process_chrom("Chr07", os.path.join(args.outdir, "check_Chr07_orientation.tsv"), selected_sites)
    process_chrom("Chr12", os.path.join(args.outdir, "check_Chr12_orientation.tsv"), selected_sites)
    
    # 4. Generate Markdown report
    md_file = os.path.join(args.outdir, "check_orientation_report.md")
    print(f"Generating MD report {md_file}...")
    with open(md_file, "w", encoding="utf-8") as f:
        f.write("# BSA 亲本与混池方向性验证报告\n\n")
        f.write("此报告通过对极端信号区段位点的等位基因深度（AD）进行追溯，验证双亲(C6/C3)与表型混池(H/L)之间等位基因的同向性。\n\n")
        
        f.write("## 1. 原理与逻辑说明\n\n")
        f.write("- **目标**: 确认由 C6 (高表型亲本) 纯合的等位基因，是否正确在 H (高表型混池) 中富集。\n")
        f.write("- **算法**:\n")
        f.write("  - $H\\_SNP\\_index = \\frac{C6\\_Allele\\_Count\\_in\\_H}{H\\_Total\\_DP}$\n")
        f.write("  - $L\\_SNP\\_index = \\frac{C6\\_Allele\\_Count\\_in\\_L}{L\\_Total\\_DP}$\n")
        f.write("  - $DELTA = H\\_SNP\\_index - L\\_SNP\\_index$\n")
        f.write("- **判断依据**: \n")
        f.write("  - 若 DELTA 整体呈现明显正峰群，代表 H 中 C6 等位基因变多、L 中 C6 等位基因变少，此时 C6-H 对应**正确**。\n")
        f.write("  - 如果发现颠倒（预期正峰结果算出来是负峰），则说明 high-parent 参数指定错误。\n\n")
        
        f.write("## 2. 输入输出说明\n\n")
        f.write(f"- **高亲本(C6) / 低亲本(C3)**: `{args.sample_c6}` / `{args.sample_c3}`\n")
        f.write(f"- **高混池(H) / 低混池(L)**: `{args.sample_h}` / `{args.sample_l}`\n")
        f.write("- **Sites 参考文件**: `...sites.tsv` (已过滤提取位点)\n")
        f.write("- **VCF 原始文件**: `...joint.vcf.gz` (提取原始 AD)\n")
        f.write(f"- **输出目录**: `{args.outdir}`\n")
        f.write("  - `check_Chr07_orientation.tsv`: Chr07 负峰检验表格\n")
        f.write("  - `check_Chr12_orientation.tsv`: Chr12 正峰检验表格\n\n")
        
        f.write("## 3. 结果解读示例\n\n")
        f.write("请对照生成的 TSV 表格查看 `HP_ALLELE` 和 `H_AD`、`L_AD` 的数值占比。\n")
        f.write("1. **如果 C6_STATE 为 0 (REF)**: 那么 C6 的等位基因就是 REF。\n")
        f.write("   - 观察 `H_AD=REF,ALT` -> H 中的 REF reads 是否显著多于 L 中的 REF reads。\n")
        f.write("   - 计算得出的 H_SNP_INDEX 应该接近 sites.tsv 里的 SITES_HIGH_SNP。\n")
        f.write("2. **如果 C6_STATE 为 1 (ALT)**: 那么 C6 的等位基因就是 ALT。\n")
        f.write("   - 观察 `H_AD=REF,ALT` -> H 中的 ALT reads 是否显著多于 L 中的 ALT reads。\n")
        f.write("3. **综合对比**: `DELTA_SNP_INDEX`（本机重新通过AD公式计算）与 `SITES_DELTA`（分析流程原始结果）如果不差毫厘，且正负符合 Chr07(-峰) 和 Chr12(+峰) 的特点，则说明：**你的分析流程中高低亲本方向没有搞反，计算逻辑百分之百坚实可靠！**\n\n")
        
    print("Done!")

if __name__ == "__main__":
    main()
