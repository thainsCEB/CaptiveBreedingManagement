#!/usr/bin/env python3
"""
Merge SRA accession FASTQs by sample name, with optional cleaning of FASTQ headers,
and a single per-sample progress bar written to the log file.

Features
--------
- Output names: <sample>_1.fastq.gz / <sample>_2.fastq.gz (paired) or <sample>.fastq.gz (single)
- Cleaning flag (-c/--clean-fastq) to normalize FASTQ headers:
    1) Remove 'length=.*$' from sequence headers.
    2) Remove '@SRR.*sra.[0-9]*' from sequence headers (case-insensitive).
    3) Replace quality header line with a bare '+'.
- Force flag (-f/--force) to overwrite existing outputs.
- **Per-sample** progress bar in the log: one bar per sample that advances as each source FASTQ
  is merged (across both _1 and _2 streams, or single-end).

Log progress bar format
-----------------------
[##########..........] 10/20 files (50.0%) :: SampleA
A line is appended after each source file completes merging for that sample.
"""

import argparse
import csv
import gzip
import io
import os
import re
import sys
from pathlib import Path
from typing import Callable, Dict, List, Tuple

# -----------------------------
# Patterns & Utilities
# -----------------------------

# Accessions
ACC_RE = re.compile(r'\b(?:SRR|ERR|DRR)\d+\b', re.IGNORECASE)

# Cleaning regexes
SEQ_HEADER_RE_LEN = re.compile(r'\s+length=.*$', re.IGNORECASE)
SEQ_HEADER_RE_SRR = re.compile(r'@SRR.*sra\.[0-9]*', re.IGNORECASE)

BAR_WIDTH = 20  # characters for the progress bar

def naturalsort_key(s: str):
    """Natural sort key (e.g., file1, file2, file10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

def detect_read_pair_side(name: str) -> str:
    """
    Heuristically detect 'R1' or 'R2' from a filename. Return '' if unknown.
    """
    base = name.lower()
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
    """Create parent dirs and open a gzipped output stream."""
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, 'wb', compresslevel=6)

def _progress_line(done: int, total: int, label: str) -> str:
    if total <= 0:
        total = 1
    frac = min(max(done / total, 0.0), 1.0)
    filled = int(BAR_WIDTH * frac)
    bar = '#' * filled + '.' * (BAR_WIDTH - filled)
    pct = f"{frac * 100:.1f}%"
    return f"[{bar}] {done}/{total} files ({pct}) :: {label}"

def stream_concat_to_gz(inputs: List[Path], outpath: Path, clean: bool,
                        file_done_cb: Callable[[], None] | None):
    """
    Concatenate FASTQs into a gzipped file.
    If clean=True, normalize FASTQ headers.
    Calls file_done_cb() once after finishing each source file.
    """
    if clean:
        with open_out_gz(outpath) as gzout:
            for f in inputs:
                opener = gzip.open if f.suffix == '.gz' else open
                with opener(f, 'rt', encoding='utf-8', errors='replace') as fin:
                    while True:
                        header = fin.readline()
                        if not header:
                            break
                        seq = fin.readline()
                        plus = fin.readline()
                        qual = fin.readline()

                        # guard against malformed last record
                        if not seq or not plus or not qual:
                            break

                        # Clean '@' header line
                        header_clean = SEQ_HEADER_RE_LEN.sub('', header.strip())
                        header_clean = SEQ_HEADER_RE_SRR.sub('', header_clean).strip()
                        if not header_clean.startswith('@'):
                            header_clean = '@' + header_clean

                        # Clean '+' line to be only '+'
                        plus_clean = '+\n'

                        # Write cleaned record
                        gzout.write((header_clean + '\n').encode('utf-8'))
                        gzout.write(seq.encode('utf-8'))
                        gzout.write(plus_clean.encode('utf-8'))
                        gzout.write(qual.encode('utf-8'))
                if file_done_cb:
                    file_done_cb()
    else:
        with open_out_gz(outpath) as gzout:
            for f in inputs:
                if f.suffix == '.gz':
                    with gzip.open(f, 'rb') as fin:
                        for chunk in iter(lambda: fin.read(1024 * 1024), b''):
                            gzout.write(chunk)
                else:
                    with open(f, 'rb') as fin:
                        for chunk in iter(lambda: fin.read(1024 * 1024), b''):
                            gzout.write(chunk)
                if file_done_cb:
                    file_done_cb()

def parse_mapping_file(path: Path) -> Dict[str, str]:
    """
    Read a two-column file (delimiter auto-detected among tab/space/comma).
    Returns dict: accession -> sample
    """
    text = path.read_text(encoding='utf-8', errors='replace')
    delimiters = ['\t', ',', None]  # None means split on whitespace
    mapping: Dict[str, str] = {}
    for delim in delimiters:
        try:
            reader = csv.reader(io.StringIO(text), delimiter=delim) if delim else None
            if reader is None:
                rows = [re.split(r'\s+', line.strip(), maxsplit=1)
                        for line in text.splitlines()
                        if line.strip() and not line.lstrip().startswith('#')]
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
    """Find FASTQ/FASTQ.GZ files under indir whose filenames contain the accession."""
    found: List[Path] = []
    acc_l = accession.lower()
    for root, _, files in os.walk(indir):
        for fn in files:
            lower = fn.lower()
            if (lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')) and acc_l in lower):
                found.append(Path(root) / fn)
    return sorted(found, key=lambda p: naturalsort_key(str(p)))

def group_by_sample(mapping: Dict[str, str]) -> Dict[str, List[str]]:
    """Group accessions by sample name, natural-sorted within each sample."""
    by_sample: Dict[str, List[str]] = {}
    for acc, smp in mapping.items():
        by_sample.setdefault(smp, []).append(acc)
    for smp in by_sample:
        by_sample[smp].sort(key=naturalsort_key)
    return by_sample

def merge_one_sample(sample: str, accessions: List[str], indir: Path, outdir: Path,
                     force: bool, clean: bool, log_lines: List[str]) -> Tuple[Path, Path, int]:
    """Merge FASTQs for a sample into _1/_2 or single-end; return (out1, out2, n_files_merged)."""
    files: List[Path] = []
    for acc in accessions:
        acc_files = find_fastqs_for_accession(indir, acc)
        if not acc_files:
            log_lines.append(f"[WARN] No FASTQs found for accession {acc} (sample {sample})")
        files.extend(acc_files)

    unique_files = sorted(set(files), key=lambda p: naturalsort_key(str(p)))
    if not unique_files:
        log_lines.append(f"[WARN] Sample {sample}: no files found to merge.")
        return (outdir / f"{sample}_1.fastq.gz", outdir / f"{sample}_2.fastq.gz", 0)

    r1_files, r2_files, unknown_files = [], [], []
    for f in unique_files:
        side = detect_read_pair_side(f.name)
        if side == 'R1':
            r1_files.append(f)
        elif side == 'R2':
            r2_files.append(f)
        else:
            unknown_files.append(f)

    # Treat as single-end if no obvious pairs but we have unknowns
    if unknown_files and not r2_files and not r1_files:
        r1_files = unknown_files
        unknown_files = []

    r1_files.sort(key=lambda p: naturalsort_key(p.name))
    r2_files.sort(key=lambda p: naturalsort_key(p.name))
    unknown_files.sort(key=lambda p: naturalsort_key(p.name))

    out_r1 = outdir / f"{sample}_1.fastq.gz"
    out_r2 = outdir / f"{sample}_2.fastq.gz"
    out_se = outdir / f"{sample}.fastq.gz"

    # Determine which inputs actually contribute to outputs for this sample
    if r1_files or r2_files:
        files_for_sample = r1_files + r2_files + unknown_files  # unknowns appended to _1 if present
    else:
        files_for_sample = unknown_files

    sample_total = len(files_for_sample)
    sample_done = 0
    if sample_total > 0:
        # Initial 0% line
        log_lines.append(_progress_line(0, sample_total, sample))

    # Shared callback to advance the single per-sample bar
    def file_done_cb():
        nonlocal sample_done
        sample_done += 1
        log_lines.append(_progress_line(sample_done, sample_total, sample))

    n_merged = 0

    if r1_files or r2_files:
        if r1_files:
            if out_r1.exists() and not force:
                log_lines.append(f"[SKIP] {out_r1} exists. Use --force to overwrite.")
            else:
                log_lines.append(f"[INFO] Merging _1 for {sample} -> {out_r1}")
                log_lines += [f"       + {p}" for p in r1_files]
                stream_concat_to_gz(r1_files, out_r1, clean, file_done_cb)
                n_merged += len(r1_files)
        if r2_files:
            if out_r2.exists() and not force:
                log_lines.append(f"[SKIP] {out_r2} exists. Use --force to overwrite.")
            else:
                log_lines.append(f"[INFO] Merging _2 for {sample} -> {out_r2}")
                log_lines += [f"       + {p}" for p in r2_files]
                stream_concat_to_gz(r2_files, out_r2, clean, file_done_cb)
                n_merged += len(r2_files)
        if unknown_files:
            log_lines.append(f"[WARN] {len(unknown_files)} files for {sample} had unknown read side; appending to _1:")
            log_lines += [f"       ? {p}" for p in unknown_files]
            if out_r1.exists() and not force and not r1_files:
                log_lines.append(f"[SKIP] {out_r1} exists. Use --force to overwrite to include unknown files.")
            else:
                files_to_write = (r1_files + unknown_files) if r1_files else unknown_files
                if out_r1.exists() and (r1_files and not force):
                    log_lines.append(f"[WARN] _1 already exists and --force not set; unknown files not appended.")
                else:
                    if r1_files and out_r1.exists():
                        out_r1.unlink()
                    stream_concat_to_gz(files_to_write, out_r1, clean, file_done_cb)
                    n_merged += len(unknown_files) if not r1_files else 0
    else:
        # Single-end
        if out_se.exists() and not force:
            log_lines.append(f"[SKIP] {out_se} exists. Use --force to overwrite.")
        else:
            log_lines.append(f"[INFO] Merging single-end for {sample} -> {out_se}")
            log_lines += [f"       + {p}" for p in unknown_files]
            stream_concat_to_gz(unknown_files, out_se, clean, file_done_cb)
            n_merged += len(unknown_files)

    return (out_r1, out_r2, n_merged)

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Merge SRA accession FASTQs by sample name.")
    ap.add_argument("-m", "--mapping", required=True,
                    help="Two-column file: <SRA_accession> <sample_name> (tab/space/csv).")
    ap.add_argument("-i", "--indir", required=True,
                    help="Input directory containing FASTQ/FASTQ.GZ files.")
    ap.add_argument("-o", "--outdir", default=None,
                    help="Output directory (default: <indir>/merged_raw_fastq).")
    ap.add_argument("-f", "--force", action="store_true",
                    help="Overwrite existing outputs.")
    ap.add_argument("-c", "--clean-fastq", action="store_true",
                    help="Clean FASTQ headers during merge.")
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
        _, _, n = merge_one_sample(sample, accs, indir, outdir, args.force, args.clean_fastq, log_lines)
        total_files += n

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
