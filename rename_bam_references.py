#!/usr/bin/env python3

"""
rename_bam_references.py

Rename @SQ SN (sequence names) in BAM headers using an accession→chromosome map,
processing one or many BAMs (via --bam-list). By default, missing accessions are allowed.

Examples
--------
# Single BAM
rename_bam_references.py -b in.bam -m acc2chrom.tsv -o outdir -x

# Many BAMs from a list file (one path per line, comments allowed)
rename_bam_references.py -L bamlist.txt -m acc2chrom.tsv -o outdir --threads 4

# Disallow missing accessions (fail fast)
rename_bam_references.py -b in.bam -m acc2chrom.tsv --no-allow-missing

Input map format
----------------
Two columns (tab/space-delimited):
    ACCESSION   CHROM
e.g.
    NC_000001.11    chr1
    NC_012920.1     MT

Notes
-----
• BAM only (not CRAM/SAM).
• We rewrite the header with updated SQ SN values and stream all reads to a new BAM.
• The order/lengths of SQ are preserved; only names are changed for keys present in the map.
"""

import sys
import os
import argparse
from pathlib import Path
import multiprocessing as mp
from typing import Dict, List, Tuple, Optional

try:
    import pysam
except ImportError as e:
    sys.stderr.write("[ERROR] pysam is required. Try: pip install pysam\n")
    raise

def load_map(path: str) -> Dict[str, str]:
    m = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Bad map line (need 2 columns): {line}")
            acc, chrom = parts[0], parts[1]
            m[acc] = chrom
    if not m:
        raise ValueError("Mapping file appears empty.")
    return m

def read_bamlist(path: str) -> List[str]:
    bam_paths = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            bam_paths.append(line.split()[0])
    if not bam_paths:
        raise ValueError("No BAMs found in bamlist.")
    return bam_paths

def ensure_outdir(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

def derive_outpath(bam_path: str, outdir: Optional[Path], suffix: str = ".renamed.bam") -> Path:
    p = Path(bam_path)
    out = p.with_suffix("")  # drop .bam
    out = Path(out.name)     # base name
    out = (outdir or p.parent) / f"{out}{suffix}"
    return out

def rename_one(bam_path: str, acc_map: Dict[str, str], outdir: Optional[Path],
               allow_missing: bool, index_output: bool, overwrite: bool,
               dry_run: bool, verbose: bool) -> Tuple[str, Optional[str]]:
    """
    Returns (bam_path, error_message or None)
    """
    try:
        bam_p = Path(bam_path)
        if not bam_p.exists():
            return (bam_path, "input not found")

        out_p = derive_outpath(bam_path, outdir)
        if out_p.exists() and not overwrite:
            if verbose:
                sys.stderr.write(f"[SKIP] {bam_p} -> {out_p} (exists)\n")
            return (str(out_p), None)

        if dry_run:
            sys.stderr.write(f"[DRYRUN] Would process {bam_p} -> {out_p}\n")
            return (str(out_p), None)

        with pysam.AlignmentFile(str(bam_p), "rb") as bam_in:
            header = bam_in.header.to_dict()

            # Update SQ names
            sq_list = header.get("SQ", [])
            name_changes = 0
            for sq in sq_list:
                old = sq.get("SN")
                if old in acc_map:
                    new = acc_map[old]
                    if verbose:
                        sys.stderr.write(f"[HEADER] {bam_p.name}: {old} -> {new}\n")
                    sq["SN"] = new
                    name_changes += 1
                else:
                    if not allow_missing:
                        raise KeyError(f"Reference '{old}' not in map.")
            if verbose:
                sys.stderr.write(f"[INFO] {bam_p.name}: changed {name_changes} SQ names\n")

            ensure_outdir(outdir or bam_p.parent)
            with pysam.AlignmentFile(str(out_p), "wb", header=header) as bam_out:
                for rec in bam_in:
                    bam_out.write(rec)

        if index_output:
            pysam.index(str(out_p))

        return (str(out_p), None)
    except Exception as e:
        return (bam_path, str(e))

def worker(args):
    return rename_one(*args)

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Rename BAM reference names using accession→chrom map.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    g_in = ap.add_argument_group("Inputs")
    g_in.add_argument("-b","--bam", help="Input BAM (one).", default=None)
    g_in.add_argument("-L","--bam-list", help="File with list of BAM paths (one per line).", default=None)
    g_in.add_argument("-m","--map", required=True, help="Accession→chrom mapping TSV/whitespace file.")

    g_out = ap.add_argument_group("Outputs")
    g_out.add_argument("-o","--outdir", default=None, help="Output directory (defaults to input file's dir).")
    g_out.add_argument("-s","--suffix", default=".renamed.bam", help="Suffix for output BAM name.")
    g_out.add_argument("-x","--index", action="store_true", help="Index output BAMs (.bai).")
    g_out.add_argument("-f","--overwrite", action="store_true", help="Overwrite outputs if present.")
    g_out.add_argument("-n","--dry-run", action="store_true", help="Plan only; do not write outputs.")
    g_out.add_argument("-v","--verbose", action="store_true", help="Verbose logging.")

    g_misc = ap.add_argument_group("Behavior")
    allow = g_misc.add_mutually_exclusive_group()
    allow.add_argument("--allow-missing", dest="allow_missing", action="store_true", help="Allow references not present in map (default).")
    allow.add_argument("--no-allow-missing", dest="allow_missing", action="store_false", help="Disallow missing references (fail).")
    ap.set_defaults(allow_missing=True)
    g_misc.add_argument("-t","--threads", type=int, default=1, help="Parallel BAMs to process.")

    args = ap.parse_args(argv)

    # Apply custom suffix if provided
    def _derive(bam_path: str) -> Path:
        p = Path(bam_path)
        out = p.with_suffix("")
        out = Path(out.name)
        base = (Path(args.outdir) if args.outdir else p.parent)
        return base / f"{out}{args.suffix}"

    # Collect BAMs
    bam_inputs: List[str] = []
    if args.bam:
        bam_inputs.append(args.bam)
    if args.bam_list:
        bam_inputs.extend(read_bamlist(args.bam_list))
    if not bam_inputs:
        ap.error("Provide --bam or --bam-list")

    # Validate / create outdir if specified
    outdir_path = Path(args.outdir) if args.outdir else None
    if outdir_path:
        ensure_outdir(outdir_path)

    acc_map = load_map(args.map)

    tasks = []
    for b in bam_inputs:
        tasks.append((
            b, acc_map, outdir_path, args.allow_missing,
            args.index, args.overwrite, args.dry_run, args.verbose
        ))

    if args.threads > 1 and len(tasks) > 1:
        with mp.Pool(processes=args.threads) as pool:
            results = pool.map(worker, tasks)
    else:
        results = list(map(worker, tasks))

    n_err = 0
    for out_or_in, err in results:
        if err is None:
            sys.stderr.write(f"[OK] {out_or_in}\n")
        else:
            sys.stderr.write(f"[FAIL] {out_or_in}: {err}\n")
            n_err += 1

    if n_err:
        sys.exit(1)

if __name__ == "__main__":
    main()
