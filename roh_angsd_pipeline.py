#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ROH / Heterozygosity pipeline (ANGSD + realSFS)

Key features
------------
- Single sample sheet (-S/--sample-sheet): first column = Sample, second column = BAM path.
- Output/working dir (-O/--outdir); SAFs auto-write to OUTDIR/saf/.
- Autosome mode (keeps MT, excludes X/Y/Z/W-like contigs).
- Ancestral FASTA (-a/--anc) → UNFOLDED SFS if ancestral != reference; otherwise folded.
- Robust parser for `realSFS fst stats2` (uses LAST (start,end) tuple; skips headers).
- Defaults: run ANGSD -doSaf and compute global SFS (disable with --no-run-dosaf, --no-make-sfs).
- Resume-aware; `--force` to recompute.
- FAI tiling fallback for Regions.txt when self-FST windowing fails.
- Cleans intermediate files by default; disable with --keep-intermediate.

Outputs
-------
- Per-segment table (no "Number") with **ROH_Length (Mb)**
- Per-sample **F_ROH**
- Combined: `AllSamples_ROH_Length.txt`, `AllSamples_F_ROH.txt`
- Summary: `AllSamples_ROH_Summary.txt` (Sample, n_ROH, Sum_ROH (Mb), F_ROH)

Usage examples
--------------
# Default (doSaf + make-sfs enabled)
python3 roh_angsd_pipeline.py -O ./roh_out -S samples.tsv -R reference.fa -c auto -T 16

# Unfolded SFS with an ancestral FASTA
python3 roh_angsd_pipeline.py -O ./roh_out -S samples.tsv -R reference.fa -a ancestral.fa -c auto -T 16

# Keep intermediates
python3 roh_angsd_pipeline.py -O ./roh_out -S samples.tsv -R reference.fa --keep-intermediate
"""

import argparse
import re
import sys
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional
import pandas as pd
import numpy as np
import glob


# ------------------------ utils ------------------------

def run(cmd: List[str], check: bool = True, capture_output: bool = False, text: bool = True) -> subprocess.CompletedProcess:
    """
    Run a shell command. On failure, raise RuntimeError with captured stdout/stderr
    so the user sees the real error message (e.g., from realSFS).
    """
    print("+", " ".join(cmd), flush=True)
    try:
        return subprocess.run(cmd, check=check, capture_output=capture_output, text=text)
    except subprocess.CalledProcessError as e:
        msg = []
        msg.append(f"[ERROR] Command failed with exit {e.returncode}: {' '.join(e.cmd)}")
        if e.stdout:
            msg.append("[STDOUT]")
            msg.append(e.stdout.strip())
        if e.stderr:
            msg.append("[STDERR]")
            msg.append(e.stderr.strip())
        raise RuntimeError("\n".join(msg)) from e

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def file_line_count(p: Path) -> int:
    if not p.exists():
        return 0
    n = 0
    with p.open("r") as fh:
        for _ in fh:
            n += 1
    return n

def unlink_if_exists(p: Path):
    try:
        if p.exists():
            p.unlink()
            print(f"[CLEAN] removed {p}")
    except Exception as e:
        print(f"[WARN] could not remove {p}: {e}")


# ------------------------ parsing helpers ------------------------

def parse_sliding_line_for_region(s: str) -> Tuple[str, str]:
    """
    Parse a realSFS 'fst stats2' line with possibly three tuples:
      (a,b)(c,d)(start,end)   <TAB> chr <TAB> midPos <TAB> Nsites ...
    We take the **last** (start,end) tuple, chromosome from the 2nd whitespace token.
    """
    toks = s.strip().split()
    if len(toks) < 2:
        raise ValueError(f"Malformed stats2 line (too few cols): {s!r}")
    chrom = toks[1]
    matches = list(re.finditer(r"\((\d+),\s*(\d+)\)", s))
    if not matches:
        raise ValueError(f"No (start,end) tuples in stats2 line: {s!r}")
    start, end = matches[-1].group(1), matches[-1].group(2)
    return chrom, f"{start}-{end}"

def write_regions_from_slidingwindow(sliding_file: Path, out_regions: Path, force: bool = False):
    """
    Build Regions.txt from 'fst stats2' output. Resume-safe:
      - If Regions.txt exists and not force, leave it alone.
    Skips header lines starting with 'region\\tchr' and 'nSites:'.
    Seeds each chromosome with 'chrom:1-100000', then appends one line per window.
    """
    if out_regions.exists() and not force:
        print(f"[SKIP] {out_regions} exists; use --force to rebuild.")
        return

    print(f"[INFO] Reading sliding-window file: {sliding_file}")
    with sliding_file.open() as fh, out_regions.open("w") as out:
        seen_chr = set()
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("region\tchr") or s.startswith("nSites:"):
                continue
            chrom, region = parse_sliding_line_for_region(line)
            if chrom not in seen_chr:
                out.write(f"{chrom}:1-100000\n")
                seen_chr.add(chrom)
            out.write(f"{chrom}:{region}\n")
    print(f"[OK] Wrote regions to {out_regions}")


def compute_hetero_from_sfs_line(s: str) -> float:
    """
    Heterozygosity = middle-bin / sum(SFS). For per-individual SFS from realSFS
    (folded or unfolded in this 2-entry setup), index 1 is the het bin in our usage.
    """
    vals = [float(x) for x in s.strip().split()]
    if len(vals) < 2:
        raise ValueError(f"SFS line has <2 entries: {s}")
    denom = sum(vals)
    return (vals[1] / denom) if denom > 0 else np.nan


# ------------------------ ROH helpers ------------------------

def rle_flags(flags: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if flags.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    diffs = np.diff(flags)
    run_starts = np.r_[0, np.where(diffs != 0)[0] + 1]
    run_ends = np.r_[run_starts[1:], flags.size]
    values = flags[run_starts]
    lengths = run_ends - run_starts
    return values, lengths

def call_roh_and_lengths(df: pd.DataFrame,
                         sample: str,
                         het_cutoff: float,
                         window_mb: float,
                         step_mb: float) -> Tuple[pd.DataFrame, float]:
    """
    Return per-segment ROH lengths (Mb) for `sample`, and F_ROH.
    Output columns: ["Sample", "Chromosome", "ROH_Length (Mb)"]
    """
    out_rows = []
    for chrom, sub in df.groupby("Chromosome", sort=False):
        flags = (sub["Hetero"].values < het_cutoff).astype(int)
        values, lengths = rle_flags(flags)
        for v, L in zip(values, lengths):
            if v == 1 and L > 1:  # drop singletons
                roh_len_mb = L * window_mb - (L - 1) * step_mb
                out_rows.append([sample, chrom, roh_len_mb])

    roh_df = pd.DataFrame(out_rows, columns=["Sample", "Chromosome", "ROH_Length (Mb)"])

    # denominator (Mb) for F_ROH:
    total_span_mb = 0.0
    for chrom, sub in df.groupby("Chromosome", sort=False):
        n = len(sub)
        if n > 0:
            total_span_mb += n * window_mb - max(0, n - 1) * step_mb

    f_roh = (roh_df["ROH_Length (Mb)"].sum() / total_span_mb) if total_span_mb > 0 else np.nan
    return roh_df, f_roh


# ------------------------ per-region SFS (RESUMABLE) ------------------------

def build_per_region_sfs(sample: str,
                         saf_idx: Path,
                         regions_file: Path,
                         threads: int,
                         out_dir: Path,
                         folded: bool,
                         force: bool = False) -> Path:
    """
    Resume logic:
      - If <sample>.sfs exists and not force:
         * Count lines in Regions.txt (R) and in sample.sfs (S).
         * If S == R, return (done).
         * If S < R, skip first S regions and continue appending.
         * If S > R, raise error unless force.
      - If force: delete .sfs and start from scratch.
    """
    out_sfs = out_dir / f"{sample}.sfs"
    R = file_line_count(regions_file)
    if out_sfs.exists() and not force:
        S = file_line_count(out_sfs)
        if S == R:
            print(f"[SKIP] {out_sfs} complete ({S}/{R} regions).")
            return out_sfs
        if S > R:
            raise RuntimeError(f"{out_sfs} has more lines ({S}) than Regions.txt ({R}); use --force.")
        start_idx = S
        print(f"[RESUME] {sample}: Regions {S}/{R} already processed; resuming at {S+1}.")
        mode = "a"
    else:
        if out_sfs.exists():
            out_sfs.unlink()
        start_idx = 0
        mode = "a"

    with regions_file.open() as fh, out_sfs.open(mode) as out:
        for i, line in enumerate(fh):
            if i < start_idx:
                continue
            region = line.strip()
            if not region:
                continue
            cmd = ["realSFS", str(saf_idx), "-r", region, "-P", str(threads)]
            if folded:
                cmd.extend(["-fold", "1"])
            cp = run(cmd, capture_output=True)
            out.write(cp.stdout)
            if ((i + 1) % 200) == 0:
                print(f"[INFO] {sample}: processed {i+1}/{R} regions")
    print(f"[OK] Wrote per-region SFS: {out_sfs}")
    return out_sfs


def paste_regions_and_sfs(regions_file: Path, sfs_file: Path, out_file: Path, force: bool = False):
    """Resume-safe: if out_file exists and not force, and line counts match, skip."""
    if out_file.exists() and not force:
        R = file_line_count(regions_file)
        S = file_line_count(sfs_file)
        O = file_line_count(out_file)
        if O == min(R, S):
            print(f"[SKIP] {out_file} already pasted ({O} lines).")
            return
    with regions_file.open() as fr, sfs_file.open() as fs, out_file.open("w") as out:
        for rline, sline in zip(fr, fs):
            out.write(rline.strip() + "\t" + sline)

def load_region_table(regionSFS_file: Path, sample: str) -> pd.DataFrame:
    rows = []
    with regionSFS_file.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            region_token = parts[0]
            sfs_line = parts[1] if len(parts) > 1 else ""
            m = re.match(r"([^:]+):(\d+)-(\d+)", region_token)
            if not m:
                continue
            chrom = m.group(1)
            region = f"{m.group(2)}-{m.group(3)}"
            hetero = compute_hetero_from_sfs_line(sfs_line)
            rows.append([sample, chrom, region, hetero])
    return pd.DataFrame(rows, columns=["Sample", "Chromosome", "Region", "Hetero"])


# ------------------------ contig filtering (auto) ------------------------

SEXLIKE_REGEX = re.compile(
    r"(?:^|[^A-Za-z0-9])(?:chr)?(?:X|Y|Z|W)(?:\b|[_\-].*|$)",
    re.IGNORECASE,
)
MITO_REGEX = re.compile(
    r"(?:^|.*mitochond.*|^(?:chr)?M(?:T)?$|^MT$)",
    re.IGNORECASE,
)

def build_autosome_rf_from_fai(ref_fa: Path, out_rf: Path, force: bool = False) -> Tuple[List[str], List[str]]:
    fai = Path(str(ref_fa) + ".fai")
    if (not fai.exists()) or force:
        print(f"[INFO] (Re)indexing reference with samtools faidx ...")
        run(["samtools", "faidx", str(ref_fa)])

    kept, excluded = [], []
    with open(fai, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            contig = line.split("\t", 1)[0]
            if MITO_REGEX.search(contig):      # always keep MT
                kept.append(contig); continue
            if SEXLIKE_REGEX.search(contig):   # drop sex-like
                excluded.append(contig)
            else:
                kept.append(contig)

    if not kept:
        sys.exit("[ERROR] No contigs kept for autosomes; check FASTA headers/filters.")

    if out_rf.exists() and not force:
        print(f"[SKIP] {out_rf} exists; use --force to rebuild.")
        return kept, excluded

    with out_rf.open("w") as out:
        for c in kept:
            out.write(c + "\n")

    print(f"[OK] Wrote autosome(+MT) contig list to {out_rf}")
    if excluded:
        print(f"[INFO] Excluded {len(excluded)} sex-like contigs (first 10): {excluded[:10]}")
    return kept, excluded


# ------------------------ FAI tiling fallback ------------------------

def iter_autosome_contigs_from_fai(ref_fa: Path, chrom_region: str, keep_mito: bool = True, force_index: bool = False) -> List[Tuple[str, int]]:
    fai = Path(str(ref_fa) + ".fai")
    if (not fai.exists()) or force_index:
        print(f"[INFO] Indexing reference with samtools faidx ...")
        run(["samtools", "faidx", str(ref_fa)])

    contigs: List[Tuple[str, int]] = []
    with fai.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t")
            name, length = parts[0], int(parts[1])
            low = chrom_region.lower()
            if low == "all":
                contigs.append((name, length))
            elif low == "auto":
                if keep_mito and MITO_REGEX.search(name):
                    contigs.append((name, length))
                elif SEXLIKE_REGEX.search(name):
                    continue
                else:
                    contigs.append((name, length))
            else:
                # explicit contig name (not interval)
                if ":" in chrom_region:
                    # explicit interval; caller shouldn't use this tiler
                    continue
                if name == chrom_region or name == f"chr{chrom_region}" or chrom_region == name.replace("chr", ""):
                    contigs.append((name, length))
    return contigs

def build_regions_from_fai(ref_fa: Path, chrom_region: str, win_bp: int, step_bp: int, out_regions: Path):
    """
    Fallback Regions.txt: tile allowed contigs into [start,end] with given win/step.
    Seeds each contig with 'contig:1-100000' to mimic prior convention.
    """
    contigs = iter_autosome_contigs_from_fai(ref_fa, chrom_region, keep_mito=True, force_index=False)
    if not contigs:
        raise RuntimeError(f"No contigs found for chrom-region={chrom_region}.")
    with out_regions.open("w") as out:
        for name, length in contigs:
            # seed line
            out.write(f"{name}:1-100000\n")
            # tile
            start = 1
            while start <= length:
                end = min(length, start + win_bp - 1)
                out.write(f"{name}:{start}-{end}\n")
                start += step_bp
    print(f"[OK] Wrote tiled regions to {out_regions}")


# ------------------------ ANGSD doSaf (RESUMABLE) ------------------------

def maybe_run_angsd_doSaf(sample: str,
                          bam: Path,
                          ref_fa: Path,
                          anc_fa: Optional[Path],
                          chrom_region: str,
                          sites_file: Optional[Path],
                          threads: int,
                          outprefix: Path,
                          workdir: Path,
                          force: bool = False):
    """
    Resume-safe:
      - If outprefix.saf.idx exists and not force, skip.
    chrom_region:
      - 'all'  -> no -r/-rf
      - 'auto' -> build autosomes.rf and use -rf
      - other  -> use -r "<chrom_region>"
    Ancestral:
      - If anc_fa and anc_fa != ref_fa → we still pass -anc anc_fa (unfolded downstream).
      - Else pass -anc ref_fa.
    """
    saf_idx = Path(str(outprefix) + ".saf.idx")
    if saf_idx.exists() and not force:
        print(f"[SKIP] {saf_idx} exists; use --force to regenerate.")
        return

    anc_arg = anc_fa if (anc_fa and anc_fa != ref_fa) else ref_fa

    cmd = [
        "angsd",
        "-i", str(bam),
        "-anc", str(anc_arg),
        "-ref", str(ref_fa),
        "-doSaf", "1",
        "-GL", "1",
        "-doCounts", "1",
        "-uniqueOnly", "1",
        "-remove_bads", "1",
        "-C", "50",
        "-baq", "2",
        "-minMapQ", "30",
        "-minQ", "30",
        "-P", str(threads),
        "-out", str(outprefix)
    ]

    if chrom_region:
        low = chrom_region.lower()
        if low == "all":
            pass
        elif low == "auto":
            rf_path = workdir / "autosomes.rf"
            build_autosome_rf_from_fai(ref_fa, rf_path, force=force)
            cmd.extend(["-rf", str(rf_path)])
        else:
            cmd.extend(["-r", chrom_region])

    if sites_file is not None:
        cmd.extend(["-sites", str(sites_file)])

    run(cmd)


# ------------------------ sample sheet ------------------------

def read_sample_sheet(sheet_path: Path) -> Tuple[List[str], List[Path]]:
    """
    Read a sample sheet where the first column is sample name and the second column is BAM path.
    Supports tab/comma/space delimiters. Skips empty lines and lines starting with '#'.
    Allows a header row if it contains 'sample' and 'bam'.
    """
    samples: List[str] = []
    bams: List[Path] = []
    with sheet_path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # split by tab/comma/space
            parts = re.split(r"[,\t ]+", s)
            if len(parts) < 2:
                # tolerate a header line like "Sample\tBAM" on first row
                if len(samples) == 0 and ("sample" in s.lower() and "bam" in s.lower()):
                    continue
                raise ValueError(f"Bad line in sample sheet (need 2 columns): {line!r}")
            sample, bam = parts[0], parts[1]
            samples.append(sample)
            bams.append(Path(bam))
    if len(samples) == 0:
        raise ValueError("No samples found in sample sheet.")
    return samples, bams


# ------------------------ cleaning ------------------------

def clean_intermediates(wd: Path, prefix: str, samples: List[str], do_clean: bool):
    if not do_clean:
        print("[INFO] Keeping intermediate files (user requested).")
        return
    print("[INFO] Cleaning intermediate files ...")
    # Global intermediates
    unlink_if_exists(wd / "test.ml")
    unlink_if_exists(wd / "slidingwindow")
    # any files generated by realSFS fst index/stats2 using prefix
    for f in glob.glob(str(wd / f"{prefix}.fst.*")):
        unlink_if_exists(Path(f))
    # per-sample intermediates
    for s in samples:
        unlink_if_exists(wd / f"{s}.sfs")
        unlink_if_exists(wd / f"{s}.regionSFS")
        unlink_if_exists(wd / f"{s}.global.sfs")
    print("[INFO] Intermediate cleanup complete.")


# ------------------------ main ------------------------

def main():
    ap = argparse.ArgumentParser(description="ROH/Heterozygosity pipeline (ANGSD/realSFS) — single sample sheet, autosome filtering, ancestral/unfolded support, resume, FAI fallback, and cleanup.")
    ap.add_argument("-O", "--outdir", required=True, help="Output (working) directory for this run. SAFs are written to OUTDIR/saf/.")
    ap.add_argument("-S", "--sample-sheet", required=True, help="Sample sheet: first column = sample name, second column = BAM path.")
    ap.add_argument("-R", "--ref", help="Reference FASTA (required for doSaf or --chrom-region auto).")
    ap.add_argument("-a", "--anc", help="Ancestral FASTA (optional). If provided and different from --ref, SFS is UNFOLDED.")
    ap.add_argument("-x", "--sites", help="Optional sites file (panel).")
    ap.add_argument("-c", "--chrom-region", default="all",
                    help="Region for ANGSD -r. Use 'all' (default), 'auto' (exclude X/Y/Z/W; include MT), or any explicit region string.")
    ap.add_argument("-T", "--threads", type=int, default=8, help="Threads for ANGSD/realSFS.")
    # defaults inverted: enabled by default; disable with --no-*
    ap.add_argument("--no-run-dosaf", dest="run_dosaf", action="store_false", help="Disable ANGSD -doSaf step (default: enabled).")
    ap.add_argument("--no-make-sfs", dest="make_sfs", action="store_false", help="Disable global SFS/heterozygosity (default: enabled).")
    ap.set_defaults(run_dosaf=True, make_sfs=True)
    ap.add_argument("-k", "--win-bp", type=int, default=100_000, help="Window size (bp) for stats2/tiling.")
    ap.add_argument("-p", "--step-bp", type=int, default=50_000, help="Step (bp) for stats2/tiling.")
    ap.add_argument("-H", "--het-cutoff", type=float, default=0.001435663, help="Per-window heterozygosity cutoff to call ROH.")
    ap.add_argument("-f", "--sliding-self-prefix", default="test", help="Prefix for self-Fst index (e.g., 'test').")
    ap.add_argument("-g", "--regions-file", default="Regions.txt", help="Filename for Regions.txt (created in OUTDIR).")
    ap.add_argument("--force", action="store_true", help="Recompute steps even if outputs exist (disables resume for those steps).")
    ap.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate files (default: cleaned).")
    args = ap.parse_args()

    # Output / working directory
    wd = Path(args.outdir).resolve(); ensure_dir(wd)
    # SAF dir is always a subdir of outdir
    saf_dir = wd / "saf"; ensure_dir(saf_dir)

    # Parse sample sheet (first col=sample, second col=bam)
    samples, bams = read_sample_sheet(Path(args.sample_sheet).resolve())

    ref = Path(args.ref).resolve() if args.ref else None
    anc = Path(args.anc).resolve() if args.anc else None
    folded = not (anc and ref and (anc != ref))  # UNFOLDED only if anc provided and differs from ref

    # doSaf (default ON)
    if args.run_dosaf:
        if not args.ref:
            sys.exit("ANGSD -doSaf needs --ref (disable with --no-run-dosaf if you already have .saf.idx).")
        if len(bams) != len(samples):
            sys.exit("Sample sheet error: number of BAMs does not match number of samples.")
        for sample, bam in zip(samples, bams):
            outprefix = saf_dir / f"{sample}"
            maybe_run_angsd_doSaf(
                sample=sample,
                bam=Path(bam),
                ref_fa=ref,
                anc_fa=anc,
                chrom_region=args.chrom_region,
                sites_file=Path(args.sites).resolve() if args.sites else None,
                threads=args.threads,
                outprefix=outprefix,
                workdir=wd,
                force=args.force
            )
    else:
        if args.chrom_region.lower() == "auto" and not args.ref:
            print("[WARN] --chrom-region auto requested but --no-run-dosaf set and --ref missing; skipping rf construction.")

    # Global SFS / heterozygosity (default ON)
    if args.make_sfs:
        hetero_global_rows = []
        for sample in samples:
            saf_idx = saf_dir / f"{sample}.saf.idx"
            if not saf_idx.exists():
                sys.exit(f"Missing {saf_idx}; run with default doSaf or remove --make-sfs.")
            out_sfs_global = wd / f"{sample}.global.sfs"
            if out_sfs_global.exists() and not args.force:
                print(f"[SKIP] {out_sfs_global} exists.")
                s = out_sfs_global.read_text().strip()
                hetero = compute_hetero_from_sfs_line(s) if s else np.nan
            else:
                cmd = ["realSFS", str(saf_idx), "-P", str(args.threads)]
                if folded:
                    cmd.extend(["-fold", "1"])
                cp = run(cmd, capture_output=True)
                out_sfs_global.write_text(cp.stdout)
                hetero = compute_hetero_from_sfs_line(cp.stdout.strip())
            hetero_global_rows.append([sample, hetero])
        pd.DataFrame(hetero_global_rows, columns=["Sample", "Global_Heterozygosity"])\
          .to_csv(wd / "Heterozygosity.txt", sep="\t", index=False)

    # Self-index + sliding-window (resume-safe with fallback to FAI tiling)
    seed_sample = samples[0]
    seed_saf = saf_dir / f"{seed_sample}.saf.idx"
    if not seed_saf.exists():
        sys.exit(f"Missing {seed_saf}; cannot build sliding windows.")

    test_ml = wd / "test.ml"
    test_fst_idx = wd / f"{args.sliding_self_prefix}.fst.idx"
    sliding_file = wd / "slidingwindow"
    regions_path = wd / args.regions_file

    built_regions = False
    if (test_ml.exists() and test_fst_idx.exists() and sliding_file.exists()) and not args.force:
        print(f"[SKIP] Using existing sliding-window artifacts.")
        # Ensure Regions.txt exists or build from slidingwindow
        if regions_path.exists() and not args.force:
            built_regions = True
        else:
            try:
                write_regions_from_slidingwindow(sliding_file, regions_path, force=args.force)
                built_regions = True
            except Exception as e:
                print("\n[WARN] Failed to parse existing slidingwindow; falling back to FAI tiling.")
                print(str(e))
                if not ref:
                    sys.exit("Ref FASTA required for FAI tiling fallback; provide -R/--ref.")
                build_regions_from_fai(ref, args.chrom_region, args.win_bp, args.step_bp, regions_path)
                built_regions = True
    else:
        # Try the realSFS self-FST path, then fall back to FAI tiling
        try:
            cmd_ml = ["realSFS", str(seed_saf), str(seed_saf), "-P", str(args.threads)]
            if folded:
                cmd_ml.extend(["-fold", "1"])
            cp_ml = run(cmd_ml, capture_output=True)
            test_ml.write_text(cp_ml.stdout)

            run([
                "realSFS", "fst", "index",
                str(seed_saf), str(seed_saf),
                "-sfs", str(test_ml),
                "-whichFst", "1",
                "-fstout", str(wd / args.sliding_self_prefix)
            ])

            cmd_sw = ["realSFS", "fst", "stats2", str(test_fst_idx), "-win", str(args.win_bp), "-step", str(args.step_bp)]
            if folded:
                cmd_sw.extend(["-fold", "1"])
            cp_sw = run(cmd_sw, capture_output=True)
            sliding_file.write_text(cp_sw.stdout)

            write_regions_from_slidingwindow(sliding_file, regions_path, force=args.force)
            built_regions = True

        except Exception as e:
            print("\n[WARN] realSFS self-FST windowing failed; falling back to FAI tiling.")
            print(str(e))  # includes stderr from realSFS via run()
            if not ref:
                sys.exit("Ref FASTA required for FAI tiling fallback; provide -R/--ref.")
            build_regions_from_fai(ref, args.chrom_region, args.win_bp, args.step_bp, regions_path)
            built_regions = True

    if not built_regions and not regions_path.exists():
        raise RuntimeError("Regions.txt not found; cannot continue.")

    # Per-sample SFS → heterozygosity → ROH (resume-safe)
    window_mb = args.win_bp / 1_000_000.0
    step_mb   = args.step_bp / 1_000_000.0
    all_roh_lengths, all_froh = [], []
    per_sample_summary = []  # rows: [Sample, n_ROH, Sum_ROH (Mb), F_ROH]

    def _normalize_roh_df(df_in: pd.DataFrame) -> pd.DataFrame:
        # Back-compat: if older files had 'Number' and/or 'ROH_Length' without units, fix them.
        df2 = df_in.copy()
        if "ROH_Length (Mb)" not in df2.columns and "ROH_Length" in df2.columns:
            df2 = df2.rename(columns={"ROH_Length": "ROH_Length (Mb)"})
        if "Number" in df2.columns:
            df2 = df2.drop(columns=["Number"])
        return df2

    for sample in samples:
        print(f"\n[===] SAMPLE {sample} [===]")
        saf_idx = saf_dir / f"{sample}.saf.idx"
        if not saf_idx.exists():
            sys.exit(f"Missing {saf_idx}")

        # Resume per-region SFS
        sample_sfs = build_per_region_sfs(
            sample=sample,
            saf_idx=saf_idx,
            regions_file=regions_path,
            threads=args.threads,
            out_dir=wd,
            folded=folded,
            force=args.force
        )

        # Paste (resume-safe by line-count)
        regionSFS_file = wd / f"{sample}.regionSFS"
        paste_regions_and_sfs(regions_path, sample_sfs, regionSFS_file, force=args.force)

        # Derive heterozygosity table
        df = load_region_table(regionSFS_file, sample)

        # Save compact table like "<sample>.RoH" (Chromosome,Region,Hetero)
        out_roh_like = wd / f"{sample}.RoH"
        if out_roh_like.exists() and not args.force:
            print(f"[SKIP] {out_roh_like} exists.")
        else:
            df[["Chromosome", "Region", "Hetero"]].to_csv(out_roh_like, sep="\t", header=False, index=False)

        # ROH & F_ROH (skip if already present unless --force)
        out_roh_lengths = wd / f"{sample}_ROH_Lengths.txt"
        out_froh        = wd / f"{sample}_F_ROH.txt"

        if out_roh_lengths.exists() and out_froh.exists() and not args.force:
            print(f"[SKIP] {out_roh_lengths} and {out_froh} exist.")
            # still collect for combined & summary
            try:
                tmp_roh = _normalize_roh_df(pd.read_csv(out_roh_lengths, sep="\t"))
                all_roh_lengths.append(tmp_roh)
                tmp_f  = pd.read_csv(out_froh, sep="\t")
                all_froh.append(list(tmp_f.iloc[0]))
                # per-sample summary
                n_roh  = len(tmp_roh)
                sum_mb = float(tmp_roh["ROH_Length (Mb)"].sum()) if n_roh > 0 else 0.0
                f_roh_val = float(tmp_f["F_ROH"].iloc[0]) if "F_ROH" in tmp_f.columns and len(tmp_f) > 0 else float("nan")
                per_sample_summary.append([sample, n_roh, sum_mb, f_roh_val])
            except Exception:
                print("[WARN] Could not read previous outputs; recomputing this sample.")
                roh_df, f_roh = call_roh_and_lengths(
                    df=df,
                    sample=sample,
                    het_cutoff=args.het_cutoff,
                    window_mb=window_mb,
                    step_mb=step_mb
                )
                roh_df.to_csv(out_roh_lengths, sep="\t", index=False)
                pd.DataFrame([[sample, f_roh]], columns=["Sample", "F_ROH"]).to_csv(out_froh, sep="\t", index=False)
                if not roh_df.empty:
                    all_roh_lengths.append(roh_df)
                all_froh.append([sample, f_roh])
                n_roh  = len(roh_df)
                sum_mb = float(roh_df["ROH_Length (Mb)"].sum()) if n_roh > 0 else 0.0
                per_sample_summary.append([sample, n_roh, sum_mb, float(f_roh)])
        else:
            roh_df, f_roh = call_roh_and_lengths(
                df=df,
                sample=sample,
                het_cutoff=args.het_cutoff,
                window_mb=window_mb,
                step_mb=step_mb
            )
            roh_df.to_csv(out_roh_lengths, sep="\t", index=False)
            pd.DataFrame([[sample, f_roh]], columns=["Sample", "F_ROH"]).to_csv(out_froh, sep="\t", index=False)
            if not roh_df.empty:
                all_roh_lengths.append(roh_df)
            all_froh.append([sample, f_roh])
            n_roh  = len(roh_df)
            sum_mb = float(roh_df["ROH_Length (Mb)"].sum()) if n_roh > 0 else 0.0
            per_sample_summary.append([sample, n_roh, sum_mb, float(f_roh)])

    # Combined outputs (safe to overwrite)
    if all_roh_lengths:
        pd.concat(all_roh_lengths, ignore_index=True).to_csv(wd / "AllSamples_ROH_Length.txt", sep="\t", index=False)
    if all_froh:
        pd.DataFrame(all_froh, columns=["Sample", "F_ROH"]).to_csv(wd / "AllSamples_F_ROH.txt", sep="\t", index=False)

    # Per-sample summary: count, sum of ROH (Mb), and F_ROH
    if per_sample_summary:
        pd.DataFrame(
            per_sample_summary,
            columns=["Sample", "n_ROH", "Sum_ROH (Mb)", "F_ROH"]
        ).to_csv(wd / "AllSamples_ROH_Summary.txt", sep="\t", index=False)

    # Cleanup (default: remove intermediates)
    clean_intermediates(
        wd=wd,
        prefix=args.sliding_self_prefix,
        samples=samples,
        do_clean=(not args.keep_intermediate)
    )

    print("\n[DONE] ROH pipeline completed (resume-aware, with FAI fallback).")
    print(f"Outputs in: {wd}")


if __name__ == "__main__":
    main()
