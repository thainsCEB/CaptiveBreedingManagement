#!/usr/bin/env python3
"""
Merge SRA accession FASTQs by sample name.

Input:
  - A two-column text file mapping SRA accessions to sample names.
      * 1st column: SRA accession (e.g., SRR12345678, ERR..., DRR...)
      * 2nd column: sample name (free text, no tabs recommended)
    The file may be tab-, space-, or comma-delimited.
  - An input directory containing FASTQ/FASTQ.GZ files whose filenames contain
    the SRA accession (common from fasterq-dump / SRA downloads, or after renaming).

Output:
  - A directory (default: <indir>/merged_raw_fastq) containing concatenated
    FASTQs per sample:
      * Paired-end: <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
      * Single-end: <sample>.fastq.gz

Notes:
  - Concatenation is performed in *compressed* byte-stream when the inputs are .gz.
    Files with mixed compression are supported; output is always gzipped.
  - Files are merged in natural-sorted order for reproducibility.
  - If a sample mixes paired and single-end files, the script will try to infer
    and merge pairs where possible and will place leftover files into R1 (single-end)
    with a warning.
"""

import argparse
import csv
import gzip
import io
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# --- Utilities ---

ACC_RE = re.compile(r'\b(?:SRR|ERR|DRR)\d+\b', re.IGNORECASE)

def naturalsort_key(s: str):
    """Natural sort: split into ints and strings."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def detect_read_pair_side(name: str) -> str:
    """
    Return 'R1', 'R2', or '' (unknown/single) based on filename patterns.
    """
    base = name.lower()
    # Common patterns: _R1, .R1, _1, .1, _read1, _unpaired, etc.
    patterns_r1 = [r'_r?1\b', r'[^a-z0-9]1[^0-9]', r'\bread1\b', r'_r1_', r'\.r1\.', r'_1_', r'_r1$', r'_1$']
    patterns_r2 = [r'_r?2\b', r'[^a-z0-9]2[^0-9]', r'\bread2\b', r'_r2_', r'\.r2\.', r'_2_', r'_r2$', r'_2$']
    for pat in patterns_r1:
        if re.search(pat, base):
            return 'R1'
    for pat in patterns_r2:
        if re.search(pat, base):
            return 'R2'
    return ''

def open_out_gz(path: Path):
    # Create parent directories and open a gzipped output stream
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, 'wb', compresslevel=6)

def stream_concat_to_gz(inputs: List[Path], outpath: Path):
    """
    Concatenate any mix of .fastq and .fastq.gz into a gzipped output file.
    We stream bytes: if input is gz -> read as binary and gunzip on the fly;
    if input is plain fastq -> read as text/bytes and compress.
    """
    with open_out_gz(outpath) as gzout:
        for f in inputs:
            if f.suffix == '.gz':
                # Decompress while streaming into gzout
                with gzip.open(f, 'rb') as fin:
                    for chunk in iter(lambda: fin.read(1024 * 1024), b''):
                        gzout.write(chunk)
            else:
                # Plain text -> compress into gz
                with open(f, 'rb') as fin:
                    for chunk in iter(lambda: fin.read(1024 * 1024), b''):
                        gzout.write(chunk)

def parse_mapping_file(path: Path) -> Dict[str, str]:
    """
    Read two-column file (delimiter auto-detected among tab/space/comma).
    Returns dict accession -> sample.
    """
    text = path.read_text(encoding='utf-8', errors='replace')
    # Try TSV first; fall back to CSV and whitespace split.
    delimiters = ['\t', ',', None]  # None = whitespace split
    mapping: Dict[str, str] = {}
    for delim in delimiters:
        try:
            reader = csv.reader(io.StringIO(text), delimiter=delim) if delim else None
            if reader is None:
                rows = [re.split(r'\s+', line.strip(), maxsplit=1) for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')]
            else:
                rows = [row for row in reader if row and not (row[0].lstrip().startswith('#'))]
            ok = True
            tmp: Dict[str, str] = {}
            for row in rows:
                if len(row) < 2:
                    ok = False
                    break
                acc = row[0].strip()
                smp = row[1].strip()
                if not ACC_RE.search(acc):
                    # Accept as-is; warn later during collection if not found.
                    pass
                tmp[acc] = smp
            if ok and tmp:
                mapping = tmp
                break
        except Exception:
            continue
    if not mapping:
        raise ValueError("Failed to parse mapping file. Ensure it has two columns: <SRA_accession> <sample_name>")
    return mapping

def find_fastqs_for_accession(indir: Path, accession: str) -> List[Path]:
    """
    Find FASTQ files under indir that contain the accession in the filename.
    """
    found: List[Path] = []
    acc_l = accession.lower()
    for root, _, files in os.walk(indir):
        for fn in files:
            lower = fn.lower()
            if (lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))
                and acc_l in lower):
                found.append(Path(root) / fn)
    return sorted(found, key=lambda p: naturalsort_key(str(p)))

def group_by_sample(mapping: Dict[str, str]) -> Dict[str, List[str]]:
    by_sample: Dict[str, List[str]] = {}
    for acc, smp in mapping.items():
        by_sample.setdefault(smp, []).append(acc)
    # natural sort accessions per sample
    for smp in by_sample:
        by_sample[smp].sort(key=naturalsort_key)
    return by_sample

def merge_one_sample(sample: str, accessions: List[str], indir: Path, outdir: Path, force: bool, log_lines: List[str]) -> Tuple[Path, Path, int]:
    """
    Merge all FASTQs for the given sample into R1/R2 (or single).
    Returns (out_r1, out_r2, n_files_merged). Missing outputs are returned as None-equivalent Paths (use exists()).
    """
    files: List[Path] = []
    for acc in accessions:
        acc_files = find_fastqs_for_accession(indir, acc)
        if not acc_files:
            log_lines.append(f"[WARN] No FASTQs found for accession {acc} (sample {sample})")
        files.extend(acc_files)

    # Deduplicate (in case the same file matched multiple accessions)
    unique_files = sorted(set(files), key=lambda p: naturalsort_key(str(p)))
    if not unique_files:
        log_lines.append(f"[WARN] Sample {sample}: no files found to merge.")
        return (outdir / f"{sample}_R1.fastq.gz", outdir / f"{sample}_R2.fastq.gz", 0)

    # Split by R1/R2
    r1_files, r2_files, unknown_files = [], [], []
    for f in unique_files:
        side = detect_read_pair_side(f.name)
        if side == 'R1':
            r1_files.append(f)
        elif side == 'R2':
            r2_files.append(f)
        else:
            unknown_files.append(f)

    # If there are unknowns and no obvious pairs, treat as single-end (append to R1)
    if unknown_files and not r2_files and not r1_files:
        r1_files = unknown_files
        unknown_files = []

    # Sort for reproducibility
    r1_files.sort(key=lambda p: naturalsort_key(p.name))
    r2_files.sort(key=lambda p: naturalsort_key(p.name))
    unknown_files.sort(key=lambda p: naturalsort_key(p.name))

    # Prepare outputs
    out_r1 = outdir / f"{sample}_R1.fastq.gz"
    out_r2 = outdir / f"{sample}_R2.fastq.gz"
    out_se = outdir / f"{sample}.fastq.gz"

    n_merged = 0

    # Merge logic
    if r1_files or r2_files:
        # Paired-end outputs
        if r1_files:
            if out_r1.exists() and not force:
                log_lines.append(f"[SKIP] {out_r1} exists. Use --force to overwrite.")
            else:
                log_lines.append(f"[INFO] Merging R1 for {sample} -> {out_r1}")
                log_lines += [f"       + {p}" for p in r1_files]
                stream_concat_to_gz(r1_files, out_r1)
                n_merged += len(r1_files)
        if r2_files:
            if out_r2.exists() and not force:
                log_lines.append(f"[SKIP] {out_r2} exists. Use --force to overwrite.")
            else:
                log_lines.append(f"[INFO] Merging R2 for {sample} -> {out_r2}")
                log_lines += [f"       + {p}" for p in r2_files]
                stream_concat_to_gz(r2_files, out_r2)
                n_merged += len(r2_files)
        if unknown_files:
            # Put leftovers into R1 with a warning
            log_lines.append(f"[WARN] {len(unknown_files)} files for {sample} had unknown read side; appending to R1:")
            log_lines += [f"       ? {p}" for p in unknown_files]
            if out_r1.exists() and not force and not r1_files:
                log_lines.append(f"[SKIP] {out_r1} exists. Use --force to overwrite to include unknown files.")
            else:
                files_to_write = (r1_files + unknown_files) if r1_files else unknown_files
                # If we already wrote R1 above, append unknowns by re-writing with combined list.
                if out_r1.exists() and (r1_files and not force):
                    log_lines.append(f"[WARN] R1 already exists and --force not set; unknown files not appended.")
                else:
                    # Rebuild complete R1 with both sets (force overwrite)
                    if r1_files:
                        # Overwrite to include both; ensure we don't double-count
                        if out_r1.exists():
                            out_r1.unlink()
                    stream_concat_to_gz(files_to_write, out_r1)
                    n_merged += len(unknown_files) if not r1_files else 0  # avoid double counting
    else:
        # Single-end
        if out_se.exists() and not force:
            log_lines.append(f"[SKIP] {out_se} exists. Use --force to overwrite.")
        else:
            log_lines.append(f"[INFO] Merging single-end for {sample} -> {out_se}")
            log_lines += [f"       + {p}" for p in unknown_files]
            stream_concat_to_gz(unknown_files, out_se)
            n_merged += len(unknown_files)

    return (out_r1, out_r2, n_merged)

def main():
    ap = argparse.ArgumentParser(description="Merge SRA accession FASTQs by sample name.")
    ap.add_argument("-m", "--mapping", required=True,
                    help="Two-column file: <SRA_accession> <sample_name> (tab/space/csv).")
    ap.add_argument("-i", "--indir", required=True,
                    help="Input directory containing FASTQ/FASTQ.GZ files (filenames include accession).")
    ap.add_argument("-o", "--outdir", default=None,
                    help="Output directory (default: <indir>/merged_raw_fastq).")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing outputs.")
    args = ap.parse_args()

    indir = Path(args.indir).resolve()
    if not indir.is_dir():
        sys.exit(f"[ERROR] Input directory not found: {indir}")

    mapping_path = Path(args.mapping).resolve()
    if not mapping_path.is_file():
        sys.exit(f"[ERROR] Mapping file not found: {mapping_path}")

    outdir = Path(args.outdir).resolve() if args.outdir else (indir / "merged_raw_fastq")
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        mapping = parse_mapping_file(mapping_path)
    except Exception as e:
        sys.exit(f"[ERROR] {e}")

    by_sample = group_by_sample(mapping)

    log_lines: List[str] = []
    total_files = 0
    total_samples = 0

    for sample, accs in sorted(by_sample.items(), key=lambda kv: naturalsort_key(kv[0])):
        total_samples += 1
        r1_out, r2_out, n = merge_one_sample(sample, accs, indir, outdir, args.force, log_lines)
        total_files += n

    # Write a simple log
    log_path = outdir / "merge_fastqs.log"
    with open(log_path, "w", encoding="utf-8") as lf:
        lf.write("\n".join(log_lines) + ("\n" if log_lines else ""))
        lf.write(f"[SUMMARY] Samples processed: {total_samples}\n")
        lf.write(f"[SUMMARY] Input files merged: {total_files}\n")
        lf.write(f"[SUMMARY] Output directory: {outdir}\n")

    print(f"[DONE] Samples processed: {total_samples}")
    print(f"[DONE] Input files merged: {total_files}")
    print(f"[DONE] Outputs in: {outdir}")
    print(f"[DONE] Log: {log_path}")

if __name__ == "__main__":
    main()
