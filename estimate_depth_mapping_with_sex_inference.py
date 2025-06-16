#!/usr/bin/env python3
"""
Infer chromosomal sex from sequencing coverage ratios from read alignments,
calculate mean depth per sample, extract duplication rate, and parse mapping stats.
"""

__author__ = "Taylor Hains"

import math
import argparse
import pandas as pd
import numpy as np
import os
from numpy import isnan
import glob

def extract_sample_name(file_path):
    base = os.path.basename(file_path)
    return base.split('.')[0]

def import_mosdepth_summary(summary_path):
    return pd.read_csv(summary_path, sep='\t')

def coverage_ratio(df, chrom1, chrom2):
    row1 = df[df.chrom == chrom1]
    row2 = df[df.chrom == chrom2]
    if row2.empty or math.isclose(0.0, row2['mean'].values[0], abs_tol=1e-9):
        return np.nan
    if row1.empty or math.isclose(0.0, row1['mean'].values[0], abs_tol=1e-9):
        return 0.0
    return row1['mean'].values[0] / row2['mean'].values[0]

def infer_sex(ratio, chrom1, chrom2):
    threshold = 0.7  # Set threshold to 0.7 as requested

    if isnan(ratio):
        return 'UNCERTAIN'
    elif ratio < threshold:
        return 'FEMALE'
    elif ratio > threshold:
        return 'MALE'
    else:
        return 'UNCERTAIN'

def extract_duplication_rate(metrics_file):
    if not os.path.exists(metrics_file):
        return np.nan
    with open(metrics_file) as f:
        for line in f:
            if line.startswith("LIBRARY"):
                headers = line.strip().split('\t')
                values = next(f).strip().split('\t')
                try:
                    idx = headers.index("PERCENT_DUPLICATION")
                    return float(values[idx]) * 100  # Convert to percentage
                except (ValueError, IndexError):
                    return np.nan
    return np.nan

def parse_flagstat(flagstat_file):
    if not os.path.exists(flagstat_file):
        return np.nan, np.nan, np.nan
    with open(flagstat_file) as f:
        lines = f.readlines()
    try:
        total_reads = int(lines[0].split()[0])
        mapped_reads = int(lines[6].split()[0])
        mapped_percent = (mapped_reads / total_reads * 100) if total_reads > 0 else np.nan
    except (IndexError, ValueError):
        total_reads = np.nan
        mapped_reads = np.nan
        mapped_percent = np.nan
    return total_reads, mapped_reads, mapped_percent

def load_file_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]

def process_summary_files(summary_files, metrics_dir, mapping_dir, sex_chrom, ref_chrom, summarize=False):
    rows = []
    for summary_path in summary_files:
        sample = extract_sample_name(summary_path)
        df = import_mosdepth_summary(summary_path)

        total_cov = df[df.chrom == "total"]['mean'].values[0] if not df[df.chrom == "total"].empty else np.nan
        ratio = coverage_ratio(df, sex_chrom, ref_chrom)
        inferred_sex = infer_sex(ratio, sex_chrom, ref_chrom)

        min_depth = math.floor(total_cov / 3) if not isnan(total_cov) else np.nan
        mean_depth = math.ceil(total_cov) if not isnan(total_cov) else np.nan
        max_depth = math.ceil(total_cov * 2) if not isnan(total_cov) else np.nan

        metrics_file = os.path.join(metrics_dir, f"{sample}.metrics.txt")
        duplication_rate = extract_duplication_rate(metrics_file)

        flagstat_file = os.path.join(mapping_dir, f"{sample}.flagstat.txt")
        total_reads, mapped_reads, mapped_percent = parse_flagstat(flagstat_file)

        rows.append({
            'sample': sample,
            'ratio': ratio,
            'infer': inferred_sex,
            'total_coverage': total_cov,
            'mean_depth': mean_depth,
            'min_depth': min_depth,
            'max_depth': max_depth,
            'duplication_percent': duplication_rate,
            'total_reads': total_reads,
            'mapped_reads': mapped_reads,
            'mapped_percent': mapped_percent
        })

    df = pd.DataFrame(rows)

    if summarize and not df.empty:
        mean_cov = df['total_coverage'].mean()
        df = pd.concat([df, pd.DataFrame([{
            'sample': 'SUMMARY',
            'ratio': '',
            'infer': '',
            'total_coverage': mean_cov,
            'mean_depth': math.ceil(mean_cov),
            'min_depth': math.floor(mean_cov / 3),
            'max_depth': math.ceil(mean_cov * 2),
            'duplication_percent': '',
            'total_reads': '',
            'mapped_reads': '',
            'mapped_percent': ''
        }])], ignore_index=True)

    return df

def main(args):
    output_dir = args.outdir
    stats_dir = os.path.join(output_dir, "stats")
    depth_dir = os.path.join(stats_dir, "depth")
    mapping_dir = os.path.join(stats_dir, "mapping")
    metrics_dir = os.path.join(stats_dir, "metrics")

    os.makedirs(depth_dir, exist_ok=True)
    os.makedirs(mapping_dir, exist_ok=True)
    os.makedirs(metrics_dir, exist_ok=True)

    ref_chrom = "total"  # default ref chrom

    if args.summary_only:
        summary_files = sorted(glob.glob(os.path.join(depth_dir, "*.mosdepth.summary.txt")))
        if not summary_files:
            raise FileNotFoundError(f"No mosdepth summary files found in {depth_dir}")
        df = process_summary_files(summary_files, metrics_dir, mapping_dir, args.sex_line, ref_chrom, summarize=args.summarize)

    else:
        import subprocess
        if not args.bamlist:
            raise ValueError("--bamlist is required when not using --summary-only")

        bam_files = load_file_list(args.bamlist)
        summary_files = []

        for bam in bam_files:
            sample = extract_sample_name(bam)
            prefix = os.path.join(depth_dir, sample)

            print(f"Running mosdepth on {bam}...")
            subprocess.run(["mosdepth", "-x", "-t", str(args.threads), prefix, bam], check=True)

            # Remove unwanted mosdepth intermediate files
            for ext in [".per-base.bed.gz", ".per-base.bed.gz.csi", ".mosdepth.global.dist.txt"]:
                try:
                    os.remove(prefix + ext)
                except FileNotFoundError:
                    pass

            summary_file = prefix + ".mosdepth.summary.txt"
            summary_files.append(summary_file)

            print(f"Running samtools flagstat on {bam}...")
            flagstat_out = os.path.join(mapping_dir, f"{sample}.flagstat.txt")
            with open(flagstat_out, "w") as f:
                subprocess.run(["samtools", "flagstat", bam], stdout=f, check=True)

        df = process_summary_files(summary_files, metrics_dir, mapping_dir, args.sex_line, ref_chrom, summarize=args.summarize)

    columns = ['sample', 'ratio', 'infer', 'total_coverage', 'mean_depth', 'min_depth', 'max_depth',
               'duplication_percent', 'total_reads', 'mapped_reads', 'mapped_percent']

    output_df = df[columns].copy()
    output_df = output_df.round({
        "ratio": 4,
        "total_coverage": 4,
        "duplication_percent": 2,
        "mapped_percent": 2
    })

    if args.output:
        output_df.to_csv(args.output, sep='\t', index=False)
        print(f"Summary table written to {args.output}")
    else:
        print(output_df.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b", "--bamlist", help="List of BAM files to process")
    parser.add_argument("--sex-line", required=True, help="Chromosome for sex ratio numerator")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Threads for mosdepth")
    parser.add_argument("-s", "--summarize", action="store_true", help="Add summary line with means")
    parser.add_argument("-o", "--output", help="Output file path")
    parser.add_argument("-d", "--outdir", help="Output directory")
    parser.add_argument("--summary-only", action="store_true",
                        help="Use existing mosdepth, flagstat and metrics files; do not run new analyses")
    args = parser.parse_args()

    main(args)
