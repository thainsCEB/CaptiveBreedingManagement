#!/usr/bin/env python3
"""
Extract soft-masked (lowercase) intervals from a FASTA and output a merged BED.

- Soft-masked bases are lowercase letters (a-z), including IUPAC ambiguity codes.
- Outputs 0-based, half-open BED intervals.
- Requires: bedtools in PATH.

Usage:
  python softmask_to_bed.py -i ref.softmasked.fa[.gz] -o softmasked.merged.bed
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import IO, Optional

LOWER_RE = re.compile(r'[a-z]')  # any lowercase letter = softmasked

def open_maybe_gzip(path: str) -> IO[bytes]:
    if path.endswith(('.gz', '.bgz')):
        return gzip.open(path, 'rb')
    f = open(path, 'rb')
    head = f.peek(2) if hasattr(f, 'peek') else f.read(2)
    if not hasattr(f, 'peek'):
        f.seek(0)
    if head.startswith(b'\x1f\x8b'):
        f.close()
        return gzip.open(path, 'rb')
    return f

def write_run(out_fh, chrom: str, start0: int, end0: int):
    if end0 > start0:
        out_fh.write(f"{chrom}\t{start0}\t{end0}\n")

def extract_softmask_intervals(fasta_path: str, unsorted_bed_path: str):
    """
    Scan FASTA; write unsorted BED of lowercase runs to unsorted_bed_path.
    """
    with open(unsorted_bed_path, 'w') as out_bed:
        with open_maybe_gzip(fasta_path) as fh:
            chrom: Optional[str] = None
            pos = 0  # 0-based index along current chrom
            in_run = False
            run_start = 0

            for raw in fh:
                line = raw.decode('utf-8', errors='ignore').rstrip('\n\r')
                if not line:
                    continue

                if line.startswith('>'):
                    if in_run and chrom is not None:
                        write_run(out_bed, chrom, run_start, pos)
                    chrom = line[1:].split()[0]
                    pos = 0
                    in_run = False
                    continue

                for ch in line:
                    if ch == '\n' or ch == '\r':
                        continue
                    if LOWER_RE.match(ch):
                        if not in_run:
                            in_run = True
                            run_start = pos
                    else:
                        if in_run:
                            write_run(out_bed, chrom, run_start, pos)
                            in_run = False
                    pos += 1

            if in_run and chrom is not None:
                write_run(out_bed, chrom, run_start, pos)

def run_cmd(cmd: list, stdout_path: str):
    try:
        with open(stdout_path, 'w') as out:
            subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[ERROR] Command failed: {' '.join(cmd)}\n{e.stderr}\n")
        sys.exit(1)

def main():
    p = argparse.ArgumentParser(description="Extract soft-masked intervals from FASTA to merged BED using bedtools.")
    p.add_argument("-i", "--fasta", required=True, help="Soft-masked reference FASTA (.fa[.gz]).")
    p.add_argument("-o", "--out-bed", required=True, help="Output BED file (merged).")
    args = p.parse_args()

    if not shutil.which("bedtools"):
        sys.stderr.write("ERROR: bedtools not found in PATH.\n")
        sys.exit(1)

    with tempfile.TemporaryDirectory(prefix="softmask2bed_") as tmpdir:
        unsorted_bed = os.path.join(tmpdir, "unsorted.bed")
        sorted_bed = os.path.join(tmpdir, "sorted.bed")

        extract_softmask_intervals(args.fasta, unsorted_bed)

        # bedtools sort -> bedtools merge (default behavior, no -d)
        run_cmd(["bedtools", "sort", "-i", unsorted_bed], stdout_path=sorted_bed)
        run_cmd(["bedtools", "merge", "-i", sorted_bed], stdout_path=args.out_bed)

    sys.stderr.write(f"[OK] Wrote merged BED: {args.out_bed}\n")

if __name__ == "__main__":
    main()
