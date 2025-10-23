#!/usr/bin/env python3
# Name: assess_mdr.py
# Title: Assess mitochondrial read proportion (MDR) for samples
# Description:
#   Computes MDR (mitochondrial read proportion) per sample using either:
#     (A) BAM/CRAM inputs: count reads mapped to mitochondrial contigs via samtools idxstats
#     (B) FASTQ inputs: subsample reads and align to a mitochondrial reference to estimate MDR
#
#   Modes:
#     1) BAM mode  : provide one or more BAM/CRAM files (-b). Requires samtools in PATH.
#     2) FASTQ mode: provide sample sheet (-s) with 3 columns: sampleID<TAB>read1<TAB>read2
#                    and a mitochondrial reference FASTA (-m/--mito-ref). Requires minimap2 (default)
#                    or bwa-mem2 (-A/--aligner), plus samtools. Optionally subsample (-n/--max-reads).
#
#   Output:
#     TSV to stdout (and optional -o/--out) with columns:
#       SampleID  TotalReads  MitoMappedReads  MDR  Mode  Notes
#
# Author: Taylor Hains
# Date: 2025-10-23
#
# Quick examples:
#   # BAM mode (auto-detect mito contigs: MT, chrM, mitochondrion, etc.)
#   python assess_mdr.py -b sample1.bam sample2.bam -o mdr.tsv
#
#   # FASTQ mode with minimap2 on a subset of reads (per sample)
#   python assess_mdr.py -s samples.tsv -m mtDNA.fasta -n 2000000 -o mdr.tsv
#
#   # FASTQ mode with bwa-mem2:
#   python assess_mdr.py -s samples.tsv -m mtDNA.fasta -A bwa-mem2 -o mdr.tsv
#
import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

MITO_GUESS_NAMES = [
    "MT", "Mt", "mt", "chrM", "chrMT", "mitochondrion", "mitochondria",
    "M", "Mito", "Mitochondrion", "mitogenome", "mitoscaf"
]

def which(x: str) -> Optional[str]:
    return shutil.which(x)

def sh(cmd: List[str], capture: bool = True) -> Tuple[int, str]:
    try:
        if capture:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
            return 0, out
        else:
            subprocess.check_call(cmd)
            return 0, ""
    except subprocess.CalledProcessError as e:
        return e.returncode, e.output if hasattr(e, "output") else ""

def guess_is_mito(contig: str) -> bool:
    c = contig.lower()
    if c in {n.lower() for n in MITO_GUESS_NAMES}:
        return True
    for kw in ["mito", "chrmt", "chrm", "mitochond", "mt-"]:
        if kw in c:
            return True
    return False

# ---------------- BAM MODE ----------------
def mdr_from_bam(bam: Path) -> Tuple[int, int, str]:
    if not which("samtools"):
        raise SystemExit("ERROR: samtools is required for BAM/CRAM mode.")
    # Ensure index exists (best-effort)
    for ext in [".bai", ".csi", ".crai"]:
        if (bam.parent / (bam.name + ext)).exists() or (bam.with_suffix(bam.suffix + ext)).exists():
            break
    else:
        sh(["samtools", "index", str(bam)], capture=False)

    rc, out = sh(["samtools", "idxstats", str(bam)], capture=True)
    if rc != 0:
        return (0, 0, "idxstats_failed")

    total = 0
    mito = 0
    for line in out.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        contig, mapped, unmapped = parts[0], parts[2], parts[3]
        try:
            mapped_i = int(mapped)
            unmapped_i = int(unmapped)
        except ValueError:
            continue
        total += mapped_i + unmapped_i
        if guess_is_mito(contig):
            mito += mapped_i
    note = "ok" if total > 0 else "empty_bam_or_failed"
    return total, mito, note

# -------------- FASTQ MODE --------------
def count_fastq_reads(fq: Path, max_reads: int) -> int:
    import gzip as gz
    opener = gz.open if str(fq).endswith(".gz") else open
    n = 0
    with opener(fq, "rt", errors="ignore") as f:
        for _ in f:
            n += 1
            if n >= max_reads * 4:
                break
    return n // 4

def build_aligner_cmd(aligner: str, mito_ref: Path, r1: Path, r2: Optional[Path]) -> List[str]:
    if aligner == "minimap2":
        if not which("minimap2"):
            raise SystemExit("ERROR: minimap2 not found in PATH (required for FASTQ mode with minimap2).")
        if r2:
            return ["minimap2", "-a", "-x", "sr", str(mito_ref), str(r1), str(r2)]
        else:
            return ["minimap2", "-a", "-x", "sr", str(mito_ref), str(r1)]
    elif aligner == "bwa-mem2":
        if not which("bwa-mem2"):
            raise SystemExit("ERROR: bwa-mem2 not found in PATH (required for FASTQ mode with bwa-mem2).")
        if r2:
            return ["bwa-mem2", "mem", str(mito_ref), str(r1), str(r2)]
        else:
            return ["bwa-mem2", "mem", str(mito_ref), str(r1)]
    else:
        raise SystemExit("ERROR: -A/--aligner must be 'minimap2' or 'bwa-mem2'.")

def mdr_from_fastq(sample_id: str, r1: Path, r2: Optional[Path], mito_ref: Path, aligner: str, max_reads: int) -> Tuple[int, int, str]:
    if not which("samtools"):
        raise SystemExit("ERROR: samtools is required (to count mapped reads from aligner output).")
    total = 0
    if r2:
        total = min(count_fastq_reads(r1, max_reads), count_fastq_reads(r2, max_reads))
    else:
        total = count_fastq_reads(r1, max_reads)
    if total == 0:
        return (0, 0, "no_reads_counted")

    aln = build_aligner_cmd(aligner, mito_ref, r1, r2)
    try:
        aln_p = subprocess.Popen(aln, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
        view_p = subprocess.Popen(["samtools", "view", "-c", "-F", "4", "-"], stdin=aln_p.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        aln_p.stdout.close()
        out_counts, err_counts = view_p.communicate()
        aln_stderr = aln_p.communicate()[1].decode("utf-8", errors="ignore")
    except Exception as e:
        return (total, 0, f"align_or_count_failed:{e}")
    try:
        mapped = int(out_counts.strip())
    except Exception:
        mapped = 0
    note = "ok"
    if "failed" in (err_counts or "").lower() or "error" in (aln_stderr or "").lower():
        note = "align_warn"
    return (total, mapped, note)

# -------------- CLI --------------
def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assess mitochondrial read proportion (MDR) from BAM/CRAM or FASTQ."
    )
    # Modes
    p.add_argument("-b", "--bam", nargs="+", help="BAM/CRAM files for BAM mode (samtools idxstats required).")
    p.add_argument("-s", "--sample-sheet", help="Sample sheet (TSV): sampleID  read1  [read2] for FASTQ mode.")
    p.add_argument("-m", "--mito-ref", help="Mitochondrial reference FASTA (required for FASTQ mode).")
    p.add_argument("-A", "--aligner", choices=["minimap2", "bwa-mem2"], default="minimap2", help="Aligner for FASTQ mode.")
    p.add_argument("-n", "--max-reads", type=int, default=2000000, help="Max reads per sample to consider in FASTQ mode.")
    p.add_argument("-o", "--out", help="Write TSV to this file (otherwise stdout).")
    p.add_argument("-t", "--mdr-threshold", type=float, default=0.001, help="Threshold (fraction) to flag low MDR (< value). e.g., 0.001 = 0.1%%")
    return p.parse_args()

def main():
    args = parse_args()
    rows = []
    header = ["SampleID", "TotalReads", "MitoMappedReads", "MDR", "Mode", "Notes"]

    if args.bam:
        for b in args.bam:
            bam = Path(b)
            sid = bam.stem
            tot, mito, note = mdr_from_bam(bam)
            mdr = (mito / tot) if tot > 0 else 0.0
            flag = "" if mdr >= args.mdr_threshold else "LOW_MDR"
            rows.append([sid, str(tot), str(mito), f"{mdr:.6f}", "BAM", f"{note} {flag}".strip()])

    elif args.sample_sheet:
        if not args.mito_ref:
            raise SystemExit("ERROR: -m/--mito-ref is required for FASTQ mode.")
        mito_ref = Path(args.mito_ref)
        with open(args.sample_sheet) as fh:
            r = csv.reader(fh, delimiter="\t")
            for parts in r:
                if not parts or parts[0].startswith("#"):
                    continue
                if len(parts) < 2:
                    continue
                sid = parts[0].strip()
                r1 = Path(parts[1].strip())
                r2 = Path(parts[2].strip()) if len(parts) >= 3 and parts[2].strip() else None
                tot, mito, note = mdr_from_fastq(sid, r1, r2, mito_ref, args.aligner, args.max_reads)
                mdr = (mito / tot) if tot > 0 else 0.0
                flag = "" if mdr >= args.mdr_threshold else "LOW_MDR"
                rows.append([sid, str(tot), str(mito), f"{mdr:.6f}", "FASTQ", f"{note} {flag}".strip()])
    else:
        raise SystemExit("ERROR: Provide either -b/--bam ... or -s/--sample-sheet ... (with -m/--mito-ref).")

    out_lines = ["\t".join(header)]
    out_lines += ["\t".join(map(str, row)) for row in rows]
    out_text = "\n".join(out_lines) + "\n"

    if args.out:
        Path(args.out).write_text(out_text)
    else:
        sys.stdout.write(out_text)

if __name__ == "__main__":
    main()
