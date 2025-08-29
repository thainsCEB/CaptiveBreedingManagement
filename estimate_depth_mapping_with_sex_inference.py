#!/usr/bin/env python3
"""
Infer chromosomal sex from sequencing coverage ratios from read alignments,
calculate mean depth per sample, extract duplication rate, and parse mapping stats.

Features
--------
- Computes per-sample mean coverage from mosdepth summaries.
- Calculates coverage ratio of a sex chromosome scaffold vs "total".
- Sex inference logic:
    * ZW (birds): FEMALE has lower Z coverage (~0.5x); MALE higher (~1x).
    * XY (mammals): FEMALE has higher X coverage (~1x); MALE lower (~0.5x).
- Flexible triggers (ONLY ONE NEEDED):
    * Provide **-S/--sex-system** alone → auto-uses default chromosome (Z for ZW, X for XY) with alias search.
    * Provide **-l/--sex-line** alone → compute ratio to that scaffold; sex reported as **UNCERTAIN** (system unknown).
    * Provide both → use both (system rules + provided scaffold, with aliasing for Z/X).
    * Provide neither → skip sex determination.
- Parses duplication rate from Picard *.metrics.txt (PERCENT_DUPLICATION).
- Parses mapping stats from `samtools flagstat` output.
- Runs mosdepth & flagstat from a BAM list, or summarizes existing outputs with `--summary-only`.

Alias search
------------
- For ZW (Z): try Z, chrZ, SUPER_Z
- For XY (X): try X, chrX, SUPER_X

Short flags
-----------
  -b/--bamlist          List of BAM files (one per line)
  -l/--sex-line         Sex-chromosome scaffold (e.g., Z, X, chrZ, chrX, SUPER_Z, SUPER_X, CMxxxxx)
  -S/--sex-system       Sex system: ZW or XY
  -t/--threads          Threads for mosdepth
  -s/--summarize        Add a summary mean row
  -o/--output           Output TSV path
  -d/--outdir           Output directory (required)
  -y/--summary-only     Use existing outputs; do not run mosdepth/flagstat

Notes
-----
- Ratio denominator is 'total' from mosdepth summary.
- Threshold for sex inference is fixed at 0.7 (empirical, user-requested).
"""

__author__ = "Taylor Hains"

import math
import argparse
import pandas as pd
import numpy as np
import os
import sys
import glob
from numpy import isnan


# ----------------------------
# I/O helpers
# ----------------------------
def extract_sample_name(file_path: str) -> str:
    """Use the basename up to first '.' as the sample name."""
    base = os.path.basename(file_path)
    return base.split('.')[0]


def import_mosdepth_summary(summary_path: str) -> pd.DataFrame:
    """Read a mosdepth summary file."""
    return pd.read_csv(summary_path, sep='\t')


def load_file_list(path: str):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


# ----------------------------
# Metrics / stats parsing
# ----------------------------
def coverage_ratio(df: pd.DataFrame, chrom1: str, chrom2: str) -> float:
    """Return mean(chrom1) / mean(chrom2). NaN if invalid."""
    row1 = df[df.chrom == chrom1]
    row2 = df[df.chrom == chrom2]
    if row2.empty:
        return np.nan
    try:
        denom = float(row2['mean'].values[0])
    except Exception:
        return np.nan
    if math.isclose(0.0, denom, abs_tol=1e-9):
        return np.nan
    if row1.empty:
        return np.nan
    try:
        numer = float(row1['mean'].values[0])
    except Exception:
        return np.nan
    if math.isclose(0.0, numer, abs_tol=1e-9):
        return 0.0
    return numer / denom


def infer_sex(ratio: float, sex_system: str, threshold: float = 0.7) -> str:
    """
    ZW: FEMALE if ratio < threshold; MALE if ratio > threshold
    XY: FEMALE if ratio > threshold; MALE if ratio < threshold
    """
    if isnan(ratio):
        return 'UNCERTAIN'
    if sex_system == 'ZW':
        if ratio < threshold:
            return 'FEMALE'
        elif ratio > threshold:
            return 'MALE'
        else:
            return 'UNCERTAIN'
    elif sex_system == 'XY':
        if ratio > threshold:
            return 'FEMALE'
        elif ratio < threshold:
            return 'MALE'
        else:
            return 'UNCERTAIN'
    return 'UNCERTAIN'


def extract_duplication_rate(metrics_file: str) -> float:
    """Parse PERCENT_DUPLICATION from a Picard-style metrics file (as percent)."""
    if not os.path.exists(metrics_file):
        return np.nan
    with open(metrics_file) as f:
        for line in f:
            if line.startswith("LIBRARY"):
                headers = line.strip().split('\t')
                values = next(f).strip().split('\t')
                try:
                    idx = headers.index("PERCENT_DUPLICATION")
                    return float(values[idx]) * 100.0
                except (ValueError, IndexError):
                    return np.nan
    return np.nan


def parse_flagstat(flagstat_file: str):
    """Parse total/mapped reads and mapped percent from samtools flagstat output."""
    if not os.path.exists(flagstat_file):
        return np.nan, np.nan, np.nan
    with open(flagstat_file) as f:
        lines = f.readlines()
    try:
        total_reads = int(lines[0].split()[0])
        mapped_reads = int(lines[6].split()[0])
        mapped_percent = (mapped_reads / total_reads * 100.0) if total_reads > 0 else np.nan
    except (IndexError, ValueError):
        total_reads = np.nan
        mapped_reads = np.nan
        mapped_percent = np.nan
    return total_reads, mapped_reads, mapped_percent


# ----------------------------
# Sex chromosome resolution
# ----------------------------
def resolve_sex_chrom_name(df: pd.DataFrame,
                           sex_system: str | None,
                           sex_line: str | None) -> str | None:
    """
    Decide which chromosome to use as the ratio numerator.

    Rules:
      - If sex_system provided but sex_line not provided:
          default to 'Z' (ZW) or 'X' (XY), then search aliases in df.chrom.
      - If both sex_system and sex_line provided:
          If sex_line is 'Z'/'X' (case-insensitive), try aliases; else require exact match.
      - If only sex_line provided (no system):
          Use EXACT name; no aliasing (system unknown).
    """
    chrom_list = df['chrom'].astype(str).tolist()

    def alias_candidates(system: str | None, line: str | None):
        if system and not line:
            return ['Z', 'chrZ', 'SUPER_Z'] if system == 'ZW' else ['X', 'chrX', 'SUPER_X']
        if system and line:
            if system == 'ZW' and line.upper() == 'Z':
                return ['Z', 'chrZ', 'SUPER_Z']
            if system == 'XY' and line.upper() == 'X':
                return ['X', 'chrX', 'SUPER_X']
            return [line]
        if line and not system:
            return [line]
        return []

    for cand in alias_candidates(sex_system, sex_line):
        if cand in chrom_list:
            return cand

    # If user supplied a non-Z/X name with system, try exact
    if sex_line and sex_line in chrom_list:
        return sex_line

    return None


# ----------------------------
# Core processing
# ----------------------------
def process_summary_files(summary_files,
                          metrics_dir,
                          mapping_dir,
                          sex_system,
                          sex_line,
                          summarize=False) -> pd.DataFrame:
    rows = []
    do_sex_det = bool(sex_system) or bool(sex_line)
    ref_chrom = "total"  # denominator

    for summary_path in summary_files:
        sample = extract_sample_name(summary_path)
        df = import_mosdepth_summary(summary_path)

        # coverage stats
        total_cov = df[df.chrom == "total"]['mean'].values[0] if not df[df.chrom == "total"].empty else np.nan
        min_depth = math.floor(total_cov / 3) if not isnan(total_cov) else np.nan
        mean_depth = math.ceil(total_cov) if not isnan(total_cov) else np.nan
        max_depth = math.ceil(total_cov * 2) if not isnan(total_cov) else np.nan

        # duplication + mapping
        metrics_file = os.path.join(metrics_dir, f"{sample}.metrics.txt")
        duplication_rate = extract_duplication_rate(metrics_file)
        flagstat_file = os.path.join(mapping_dir, f"{sample}.flagstat.txt")
        total_reads, mapped_reads, mapped_percent = parse_flagstat(flagstat_file)

        # sex determination
        ratio = np.nan
        inferred_sex = 'SKIPPED'
        sex_chr_used = ''

        if do_sex_det:
            resolved = resolve_sex_chrom_name(df, sex_system, sex_line)
            # For reporting: if nothing resolved, show what was attempted
            if resolved:
                sex_chr_used = resolved
            else:
                if sex_system and not sex_line:
                    sex_chr_used = 'Z' if sex_system == 'ZW' else 'X'
                elif sex_line:
                    sex_chr_used = f"{sex_line} (NOT FOUND)"
                else:
                    sex_chr_used = ''

            ratio = coverage_ratio(df, resolved, ref_chrom) if resolved else np.nan

            if sex_system:
                inferred_sex = infer_sex(ratio, sex_system) if resolved else 'UNCERTAIN'
                if resolved is None:
                    print(f"[WARN] Could not resolve sex chromosome in sample '{sample}' "
                          f"(system={sex_system}, line={sex_line}). Ratio set to NaN; sex=UNCERTAIN.",
                          file=sys.stderr)
            else:
                # sex_line only (no system): cannot apply MALE/FEMALE rules
                inferred_sex = 'UNCERTAIN'
                if sex_line and resolved is None:
                    print(f"[WARN] Sex scaffold '{sex_line}' not found in {sample} mosdepth summary; "
                          f"ratio set to NaN; sex=UNCERTAIN.", file=sys.stderr)
                elif sex_line:
                    print(f"[INFO] Only --sex-line provided for {sample}; reporting ratio, sex=UNCERTAIN (no system).",
                          file=sys.stderr)

        rows.append({
            'sample': sample,
            'sex_system': sex_system if sex_system else '',
            'sex_chrom_used': sex_chr_used,
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
            'sex_system': '',
            'sex_chrom_used': '',
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


def run_pipeline(args) -> pd.DataFrame:
    """
    Run mosdepth & flagstat if not summary-only, then aggregate all outputs.
    """
    output_dir = args.outdir
    stats_dir = os.path.join(output_dir, "stats")
    depth_dir = os.path.join(stats_dir, "depth")
    mapping_dir = os.path.join(stats_dir, "mapping")
    metrics_dir = os.path.join(stats_dir, "metrics")

    os.makedirs(depth_dir, exist_ok=True)
    os.makedirs(mapping_dir, exist_ok=True)
    os.makedirs(metrics_dir, exist_ok=True)

    if args.summary_only:
        summary_files = sorted(glob.glob(os.path.join(depth_dir, "*.mosdepth.summary.txt")))
        if not summary_files:
            raise FileNotFoundError(f"No mosdepth summary files found in {depth_dir}")
        return process_summary_files(
            summary_files, metrics_dir, mapping_dir,
            args.sex_system, args.sex_line, summarize=args.summarize
        )

    # else: run mosdepth & flagstat
    import subprocess
    if not args.bamlist:
        raise ValueError("--bamlist is required when not using --summary-only")

    bam_files = load_file_list(args.bamlist)
    summary_files = []

    for bam in bam_files:
        sample = extract_sample_name(bam)
        prefix = os.path.join(depth_dir, sample)

        print(f"[RUN] mosdepth on {bam}")
        subprocess.run(["mosdepth", "-x", "-t", str(args.threads), prefix, bam], check=True)

        # Remove unwanted mosdepth intermediate files
        for ext in [".per-base.bed.gz", ".per-base.bed.gz.csi", ".mosdepth.global.dist.txt"]:
            try:
                os.remove(prefix + ext)
            except FileNotFoundError:
                pass

        summary_file = prefix + ".mosdepth.summary.txt"
        summary_files.append(summary_file)

        print(f"[RUN] samtools flagstat on {bam}")
        flagstat_out = os.path.join(mapping_dir, f"{sample}.flagstat.txt")
        with open(flagstat_out, "w") as f:
            subprocess.run(["samtools", "flagstat", bam], stdout=f, check=True)

    return process_summary_files(
        summary_files, metrics_dir, mapping_dir,
        args.sex_system, args.sex_line, summarize=args.summarize
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b", "--bamlist", help="List of BAM files to process (one path per line)")
    parser.add_argument("-l", "--sex-line",
                        help=("Chromosome/scaffold used for sex determination (e.g., Z or X, or a specific contig). "
                              "When provided alone (without --sex-system), only ratio is computed and sex is UNCERTAIN."))
    parser.add_argument("-S", "--sex-system", choices=["ZW", "XY"],
                        help=("Sex chromosome system. If provided alone, uses default sex chromosome "
                              "(Z for ZW, X for XY) with alias search."))
    parser.add_argument("-t", "--threads", type=int, default=1, help="Threads for mosdepth")
    parser.add_argument("-s", "--summarize", action="store_true", help="Add summary line with means")
    parser.add_argument("-o", "--output", help="Output TSV file path")
    parser.add_argument("-d", "--outdir", required=True, help="Output directory (will create stats/ subdirs)")
    parser.add_argument("-y", "--summary-only", action="store_true",
                        help="Use existing mosdepth, flagstat, and metrics files; do not run new analyses")

    args = parser.parse_args()

    df = run_pipeline(args)

    columns = [
        'sample', 'sex_system', 'sex_chrom_used', 'ratio', 'infer',
        'total_coverage', 'mean_depth', 'min_depth', 'max_depth',
        'duplication_percent', 'total_reads', 'mapped_reads', 'mapped_percent'
    ]

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
    main()
