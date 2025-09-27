#!/usr/bin/env python3
"""
Batch ANGSD mitochondrial FASTA pipeline with per-sample temp directories.

Sample sheet format (tab-delimited, no header):
    sample_name<TAB>/absolute/or/relative/path/to/sample.bam

For each sample:
  1) Run ANGSD to make MT FASTA
  2) gunzip the FASTA
  3) strip 'N' characters from sequence ends and unwrap
  4) extract MT gene sequences with bedtools getfasta
  5) rename >MT to >{sample} in the whole mitogenome FASTA
  6) move final files to outdir/mitogenome/
  7) remove per-sample temp directory

Requirements: angsd, bedtools, python3 with Biopython.

Author: Taylor Hains (modified)
"""

import argparse
import subprocess
import os
import sys
import shutil
import pandas as pd
from Bio import SeqIO

def run_cmd(cmd, cwd=None):
    """Run a shell command and exit on failure."""
    pretty = " ".join(cmd)
    print(f"[run] {pretty}")
    result = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"[error] Command failed ({pretty})", file=sys.stderr)
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    return result.stdout

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def clean_and_unwrap_fasta(fasta_in, fasta_out):
    """Trim N from sequence ends and unwrap sequences."""
    with open(fasta_out, "w") as fout:
        for record in SeqIO.parse(fasta_in, "fasta"):
            seq = str(record.seq).upper()
            seq = seq.lstrip("N").rstrip("N")
            fout.write(f">{record.id}\n{seq}\n")

def process_sample(sample, bam, reference, bed, mt_region, threads, outdir):
    # Directories
    mitogenome_dir = ensure_dir(os.path.join(outdir, "mitogenome"))
    tmp_dir = ensure_dir(os.path.join(outdir, f"tmp_{sample}"))

    # ANGSD output prefix lives in tmp_dir
    prefix = os.path.join(tmp_dir, sample)
    fa_gz = f"{prefix}.fa.gz"
    fa = f"{prefix}.fa"
    cleaned_fa = f"{prefix}.cleaned.fa"
    mt_genes_fa_tmp = os.path.join(tmp_dir, f"{sample}.mt_genes.fa")

    print(f"\n[info] Processing {sample}")
    print(f"[info] BAM: {bam}")

    # 1) ANGSD
    angsd_cmd = [
        "angsd",
        "-i", bam,
        "-P", str(threads),
        "-doFasta", "2",
        "-doCounts", "1",
        "-minMapQ", "20",
        "-minQ", "20",
        "-remove_bads", "1",
        "-uniqueOnly", "1",
        "-only_proper_pairs", "0",
        "-ref", reference,
        "-setMinDepthInd", "2",
        "-r", mt_region,
        "-out", prefix
    ]
    run_cmd(angsd_cmd)

    # 2) gunzip
    if os.path.exists(fa_gz):
        run_cmd(["gunzip", "-f", fa_gz])
    else:
        print(f"[warn] Expected ANGSD output not found: {fa_gz}", file=sys.stderr)

    # 3) unwrap + strip Ns at ends
    if os.path.exists(fa):
        clean_and_unwrap_fasta(fa, cleaned_fa)
    else:
        print(f"[error] Missing FASTA for {sample}: {fa}", file=sys.stderr)
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return

    # 4) bedtools getfasta for MT genes
    run_cmd(["bedtools", "getfasta", "-fi", cleaned_fa, "-bed", bed, "-nameOnly", "-fo", mt_genes_fa_tmp])

    # 5) rename >MT to >{sample} in the whole mitogenome FASTA
    with open(cleaned_fa) as fin, open(cleaned_fa + ".renamed", "w") as fout:
        for line in fin:
            if line.startswith(">MT"):
                fout.write(line.replace(">MT", f">{sample}"))
            else:
                fout.write(line)
    final_fa = os.path.join(mitogenome_dir, f"{sample}.fa")
    shutil.move(cleaned_fa + ".renamed", final_fa)

    # 6) move MT genes fasta
    final_mt_genes_fa = os.path.join(mitogenome_dir, f"{sample}.mt_genes.fa")
    shutil.move(mt_genes_fa_tmp, final_mt_genes_fa)

    # 7) clean tmp
    shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"[done] {sample}")
    print(f"       mitogenome: {final_fa}")
    print(f"       MT genes:   {final_mt_genes_fa}")

def main():
    p = argparse.ArgumentParser(description="ANGSD MT pipeline (batch) with per-sample cleanup")
    p.add_argument("-s", "--sample-sheet", required=True,
                   help="Two-column TSV: sample_name<TAB>bam_path")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p.add_argument("-b", "--bed", required=True, help="BED file of MT genes")
    p.add_argument("-m", "--mt-region", required=True,
                   help="MT region name/range (e.g., MT or chrM or 'MT:1-16569')")
    p.add_argument("-t", "--threads", type=int, default=4, help="Threads for ANGSD [4]")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory (will contain 'mitogenome/' and temp dirs)")
    args = p.parse_args()

    # Validate inputs
    for path in (args.reference, args.bed):
        if not os.path.exists(path):
            print(f"[error] Not found: {path}", file=sys.stderr)
            sys.exit(1)

    ensure_dir(args.outdir)

    # Read sample sheet
    try:
        df = pd.read_csv(args.sample_sheet, sep="\t", header=None, names=["sample", "bam"])
    except Exception as e:
        print(f"[error] Could not read sample sheet: {e}", file=sys.stderr)
        sys.exit(1)

    # Iterate samples
    for _, row in df.iterrows():
        sample, bam = str(row["sample"]).strip(), str(row["bam"]).strip()
        if not sample or not bam:
            print(f"[warn] Skipping malformed row: {row.to_dict()}", file=sys.stderr)
            continue
        if not os.path.exists(bam):
            print(f"[warn] BAM not found for {sample}: {bam} (skipping)", file=sys.stderr)
            continue
        process_sample(sample, bam, args.reference, args.bed, args.mt_region, args.threads, args.outdir)

if __name__ == "__main__":
    main()
