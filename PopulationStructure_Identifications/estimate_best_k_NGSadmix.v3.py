#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Title: Run NGSadmix across K, estimate best K (Evanno), and extract best Q
# Description:
#   - Runs NGSadmix for K = [minK..maxK] with N replicates each (always with -printInfo 1).
#   - Parses log-likelihoods, computes Evanno DeltaK (faithful to R bestK_evanno), and writes outputs.
#   - Supports --evanno-only (no NGSadmix) and --ngs-only (no Evanno).
#   - Optional ggplot ΔK PNG via Rscript.
# Author: Taylor Hains
# Date: 2025-10-19
# -----------------------------------------------------------------------------

import os
import re
import shutil
import argparse
import shlex
import subprocess
import pandas as pd
import numpy as np
from typing import List, Tuple

def run_cmd(cmd: List[str], log_path: str = None) -> subprocess.CompletedProcess:
    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    print("[cmd]", cmd_str)
    if log_path:
        try:
            with open(log_path, "a") as _lf:
                _lf.write(f"[cmd] {cmd_str}\n")
        except Exception:
            pass
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
    if log_path:
        try:
            with open(log_path, "a") as _lf:
                if proc.stdout:
                    _lf.write(proc.stdout if proc.stdout.endswith("\n") else proc.stdout + "\n")
                if proc.stderr:
                    _lf.write(proc.stderr if proc.stderr.endswith("\n") else proc.stderr + "\n")
        except Exception:
            pass
    return proc

def run_ngsadmix(beagle_file: str, min_k: int, max_k: int, threads: int, replicates: int, outdir: str) -> pd.DataFrame:
    basename = os.path.basename(beagle_file).replace(".beagle.gz", "")
    os.makedirs(outdir, exist_ok=True)
    rows = []

    for k in range(min_k, max_k + 1):
        for rep in range(1, replicates + 1):
            outprefix = f"{outdir}/{basename}.K{k}.rep{rep}"
            log_file  = f"{outprefix}.log"
            print(f"[NGSadmix] K={k} rep={rep}")
            cmd = [
                "NGSadmix", "-likes", beagle_file,
                "-K", str(k),
                "-P", str(threads),
                "-o", outprefix,
                "-printInfo", "1",
                "-seed", "8",
            ]
            proc = run_cmd(cmd, log_file)
            if proc.returncode != 0:
                print(f"[warn] NGSadmix failed for K={k} rep={rep} (exit {proc.returncode})")
                if proc.stdout.strip():
                    print(proc.stdout)
                if proc.stderr.strip():
                    print(proc.stderr)

            # Parse likelihood from log
            logl = np.nan
            try:
                if os.path.exists(log_file):
                    with open(log_file) as fh:
                        txt = fh.read()
                    m = re.search(r"bestLike\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt)
                    if not m:
                        m = re.search(r"like\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt, re.IGNORECASE)
                    if m:
                        logl = float(m.group(1))
            except Exception as e:
                print(f"[warn] Could not parse likelihood from {log_file}: {e}")

            rows.append((k, rep, logl, beagle_file))

    df = pd.DataFrame(rows, columns=["K", "rep", "logL", "beagle_file"])
    df.to_csv(f"{outdir}/likelihoods.tsv", sep="\t", index=False)
    return df

def _is_sequential(int_list: List[int]) -> bool:
    if len(int_list) <= 1:
        return True
    diffs = [b - a for a, b in zip(int_list[:-1], int_list[1:])]
    return all(d == diffs[0] for d in diffs)

def _structure_summary(df: pd.DataFrame, outdir: str) -> Tuple[List[int], pd.DataFrame]:
    g = df.groupby("K")["logL"]
    summary = pd.DataFrame({
        "K": g.mean().index.astype(int),
        "variable": "L(K)",
        "value": g.mean().values,
        "sd": g.std().fillna(0).values,
    }).sort_values("K")
    summary.to_csv(f"{outdir}/structure_summary.tsv", sep="\t", index=False)
    if not summary.empty:
        best_k = int(summary.loc[summary["value"].idxmax(), "K"])
        return [best_k], summary
    return [], summary

def evanno_method(df: pd.DataFrame, outdir: str, multispecies: bool=False) -> Tuple[List[int], pd.DataFrame]:
    """
    Python port of R bestK_evanno():
      - K sequential, equal replicates per K, not all single-replicate
      - Computes L(K), L'(K), L''(K), and delta K = mean(|L''(K)|)/sd(L(K))
    Falls back to structure summary if invalid.
    Writes:
      - evanno_summary.tsv (tidy)
      - evanno_deltaK.tsv  (K, DeltaK)
    """
    if df.empty or "logL" not in df or df["logL"].isna().all():
        print("[evanno] No likelihoods available.")
        pd.DataFrame(columns=["K","DeltaK"]).to_csv(f"{outdir}/evanno_deltaK.tsv", sep="\t", index=False)
        return [], pd.DataFrame(columns=["K","DeltaK"])

    counts = df.groupby("K").size().sort_index()
    Ks = counts.index.astype(int).tolist()
    equal_reps = counts.nunique() == 1
    all_single = (counts.iloc[0] == 1) and counts.eq(1).all()
    sequential = _is_sequential(Ks)

    if (not sequential) or (not equal_reps) or all_single:
        if all_single:
            print("[evanno] Not enough information to compute Evanno statistics (one run per K).")
        else:
            print("[evanno] WARNING: Non-sequential K or unequal replicates. Reverting to structure method.")
        return _structure_summary(df, outdir)

    # Split by K preserving replicate order
    split = {int(k): grp.sort_values("rep")["logL"].to_list() for k, grp in df.groupby("K")}
    k_all = list(range(min(split.keys()), max(split.keys()) + 1))

    rows = []
    for k in k_all:
        LK = np.array(split[k], dtype=float)
        if (k > k_all[0]) and (k < k_all[-1]):
            dK  = np.array(split[k], dtype=float) - np.array(split[k-1], dtype=float)
            ddK = np.array(split[k+1], dtype=float) - 2.0*np.array(split[k], dtype=float) + np.array(split[k-1], dtype=float)
            if np.isclose(np.std(LK, ddof=1), 0.0):
                raise RuntimeError("No deviation between runs of the same K. Evanno statistics cannot be computed.")
            vals = [np.mean(LK), np.mean(dK), np.mean(np.abs(ddK)), np.mean(np.abs(ddK)) / np.std(LK, ddof=1)]
            sds  = [np.std(LK, ddof=1), np.std(dK, ddof=1), np.std(ddK, ddof=1), np.nan]
        else:
            vals = [np.mean(LK), np.nan, np.nan, np.nan]
            sds  = [np.std(LK, ddof=1), np.nan, np.nan, np.nan]

        variables = ["L(K)", "L'(K)", "L''(K)", "delta K"]
        for var, val, sd in zip(variables, vals, sds):
            rows.append((int(k), var, (float(val) if (val is not None and not np.isnan(val)) else np.nan),
                         (float(sd) if (sd is not None and (isinstance(sd, (int, float)) or np.isnan(sd))) else np.nan)))

    summary = pd.DataFrame(rows, columns=["K","variable","value","sd"]).sort_values(["K","variable"])
    summary.to_csv(f"{outdir}/evanno_summary.tsv", sep="\t", index=False)

    delta = summary[summary["variable"] == "delta K"][["K","value"]].rename(columns={"value":"DeltaK"}).dropna()
    delta.to_csv(f"{outdir}/evanno_deltaK.tsv", sep="\t", index=False)

    if delta.empty:
        print("[evanno] Insufficient data for DeltaK (need interior Ks with nonzero SD).")
        return [], delta

    if multispecies:
        best_ks = delta.sort_values("DeltaK", ascending=False).head(5)["K"].astype(int).tolist()
        print(f"[evanno] Top ΔK (multi): {best_ks}")
    else:
        best_ks = [int(delta.loc[delta["DeltaK"].idxmax(), "K"])]
        print(f"[evanno] Best K = {best_ks[0]}")
    return best_ks, delta

def extract_best_q_matrices(outdir: str, df_likelihoods: pd.DataFrame, best_ks: List[int], output_prefix: str="bestQ") -> List[str]:
    """
    For each K in best_ks, find the replicate with the highest logL and copy:
      - {basename}.K{K}.rep{rep}.qopt    -> {output_prefix}.K{K}.qopt
      - {basename}.K{K}.rep{rep}.fopt.gz -> {output_prefix}.K{K}.fopt.gz
      - {basename}.K{K}.rep{rep}.filter  -> {output_prefix}.K{K}.filter
    Returns a list of destination file paths successfully written.
    """
    saved = []
    if df_likelihoods.empty or not best_ks:
        return saved

    basename = os.path.basename(df_likelihoods["beagle_file"].iloc[0]).replace(".beagle.gz", "")
    for k in sorted(set(best_ks)):
        subset = df_likelihoods[(df_likelihoods["K"] == k) & df_likelihoods["logL"].notna()]
        if subset.empty:
            print(f"[save] No logL for K={k}; skip copying.")
            continue
        rep = int(subset.loc[subset["logL"].idxmax(), "rep"])
        base = os.path.join(outdir, f"{basename}.K{k}.rep{rep}")

        # qopt
        src = f"{base}.qopt"
        dst = os.path.join(outdir, f"{output_prefix}.K{k}.qopt")
        try:
            shutil.copyfile(src, dst)
            print(f"[save] Best Q for K={k} (rep={rep}) → {dst}")
            saved.append(dst)
        except FileNotFoundError:
            print(f"[save] Missing file: {src}")

        # fopt.gz
        src = f"{base}.fopt.gz"
        dst = os.path.join(outdir, f"{output_prefix}.K{k}.fopt.gz")
        try:
            shutil.copyfile(src, dst)
            print(f"[save] Best fopt for K={k} (rep={rep}) → {dst}")
            saved.append(dst)
        except FileNotFoundError:
            print(f"[save] Missing file: {src}")

        # filter
        src = f"{base}.filter"
        dst = os.path.join(outdir, f"{output_prefix}.K{k}.filter")
        try:
            shutil.copyfile(src, dst)
            print(f"[save] Best filter for K={k} (rep={rep}) → {dst}")
            saved.append(dst)
        except FileNotFoundError:
            print(f"[save] Missing file: {src}")

    return saved

def write_and_run_evanno_plot(rscript_bin: str, delta_tsv: str, out_png: str, width: float=6.0, height: float=4.0) -> None:
    """
    Create an R script that reads evanno_deltaK.tsv and writes a ΔK PNG using ggplot2.
    """
    import tempfile
    r_code = """
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript evanno_plot.R <evanno_deltaK.tsv> <out.png> [width] [height]")
}
in_tsv <- args[1]
out_png <- args[2]
w <- ifelse(length(args) >= 3, as.numeric(args[3]), 6.0)
h <- ifelse(length(args) >= 4, as.numeric(args[4]), 4.0)

df <- read_tsv(in_tsv, show_col_types = FALSE, progress = FALSE)
if (!all(c("K","DeltaK") %in% colnames(df))) {
  stop("Input file must contain columns: K, DeltaK")
}
df <- df %>% arrange(K)

p2 <- ggplot(df, aes(x = K, y = DeltaK, group = 1)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(df$K))) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  labs(title = "Evanno \u0394K", x = "K", y = "\u0394K") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )

ggsave(filename = out_png, plot = p2, width = w, height = h, dpi = 300, units = "in")
"""
    with tempfile.TemporaryDirectory() as td:
        r_path = os.path.join(td, "evanno_plot.R")
        with open(r_path, "w") as f:
            f.write(r_code)
        cmd = [rscript_bin, r_path, delta_tsv, out_png, str(width), str(height)]
        print(f"[plot] Running: {' '.join(cmd)}")
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            print(proc.stdout)
            print(proc.stderr)
            raise RuntimeError(f"R plotting failed with exit code {proc.returncode}")
        else:
            if proc.stdout.strip():
                print(proc.stdout)
            if proc.stderr.strip():
                print(proc.stderr)
            print(f"[plot] Wrote {out_png}")




def parse_evanno_list(list_path: str, outdir: str) -> pd.DataFrame:
    """Read a whitespace-delimited list: K rep logfile_path (3 columns).
    Parse each logfile for bestLike/like and return a DataFrame with columns K, rep, logL, beagle_file.
    The beagle_file column is synthesized from the logfile basename so downstream file extraction works.
    Also writes likelihoods.from_list.tsv in outdir.
    """
    rows = []
    import os, re
    with open(list_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                print(f"[evanno-list] Skip line (need exactly 3 columns: K rep logfile): {line.strip()}")
                continue
            try:
                k = int(parts[0]); rep = int(parts[1])
            except ValueError:
                print(f"[evanno-list] Non-integer K/rep; skip: {line.strip()}")
                continue
            log_path = parts[2]
            base_log = os.path.basename(log_path)
            if base_log.endswith(".log"):
                base_log = base_log[:-4]
            m2 = re.search(r"^(.*)\.K%s\.rep%s$" % (k, rep), base_log)
            if m2:
                basename = m2.group(1)
            else:
                basename = base_log.split(".K")[0]
            beagle_file = basename + ".beagle.gz"
            logl = None
            try:
                with open(log_path, "r") as lf:
                    txt = lf.read()
                m = re.search(r"bestLike\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt)
                if not m:
                    m = re.search(r"like\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt, re.IGNORECASE)
                if m:
                    logl = float(m.group(1))
            except Exception as e:
                print(f"[evanno-list] Could not read {log_path}: {e}")
            rows.append((k, rep, logl, beagle_file))
    df = pd.DataFrame(rows, columns=["K","rep","logL","beagle_file"])
    if not df.empty:
        outp = os.path.join(outdir, "likelihoods.from_list.tsv")
        df.to_csv(outp, sep="\t", index=False)
        print(f"[evanno-list] Wrote {outp}")
    return df

def load_or_parse_likelihoods(outdir: str, beagle_file: str, min_k: int, max_k: int, replicates: int) -> pd.DataFrame:
    """
    Load likelihoods.tsv if present; otherwise parse existing .log files in outdir.
    Returns a DataFrame with columns: K, rep, logL, beagle_file
    """
    like_path = os.path.join(outdir, "likelihoods.tsv")
    if os.path.exists(like_path):
        try:
            df = pd.read_csv(like_path, sep="\t")
            if {"K","rep","logL","beagle_file"}.issubset(df.columns):
                return df
        except Exception:
            pass
    # parse logs
    import glob
    basename = os.path.basename(beagle_file).replace(".beagle.gz", "")
    rows = []
    for k in range(min_k, max_k + 1):
        for rep in range(1, replicates + 1):
            log_file = os.path.join(outdir, f"{basename}.K{k}.rep{rep}.log")
            if not os.path.exists(log_file):
                patt = os.path.join(outdir, f"*K{k}.rep{rep}.log")
                matches = glob.glob(patt)
                if matches:
                    log_file = matches[0]
                else:
                    continue
            logl = np.nan
            try:
                with open(log_file) as fh:
                    txt = fh.read()
                m = re.search(r"bestLike\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt)
                if not m:
                    m = re.search(r"like\s*[:=]\s*(-?\d+(?:\.\d+)?)", txt, re.IGNORECASE)
                if m:
                    logl = float(m.group(1))
            except Exception:
                pass
            rows.append((k, rep, logl, beagle_file))
    df = pd.DataFrame(rows, columns=["K", "rep", "logL", "beagle_file"])
    if not df.empty:
        df.to_csv(like_path, sep="\t", index=False)
    return df

def main():
    parser = argparse.ArgumentParser(
        description="NGSadmix K scan with Evanno best-K and best .qopt extraction (no ADMIXTURE)."
    )
    parser.add_argument("-l","--likes", required=True, help="Input .beagle.gz likelihood file")
    parser.add_argument("-k","--min-k", type=int, required=True, help="Minimum K")
    parser.add_argument("-K","--max-k", type=int, required=True, help="Maximum K")
    parser.add_argument("-r","--reps", type=int, default=5, help="Replicates per K (default: 5)")
    parser.add_argument("-P","--threads", type=int, default=4, help="Threads for NGSadmix (default: 4)")
    parser.add_argument("-o","--outdir", required=True, help="Output directory")
    parser.add_argument("-m","--multispecies", action="store_true", help="Report top 5 ΔK values instead of single best K")

    # Plotting
    parser.add_argument("-p","--plot", action="store_true", help="Generate Evanno ΔK plot via Rscript")
    parser.add_argument("-R","--rscript", default="Rscript", help="Path to Rscript binary (default: Rscript)")
    parser.add_argument("-w","--plot-width", type=float, default=6.0, help="Plot width in inches (default: 6)")
    parser.add_argument("-H","--plot-height", type=float, default=4.0, help="Plot height in inches (default: 4)")
    parser.add_argument("-O","--plot-out", default=None, help="Optional path for the PNG plot (default: <outdir>/evanno_deltaK.png)")

    # Exclusive modes
    parser.add_argument("-E","--evanno-only", action="store_true", help="Only compute Evanno/ΔK from existing .log files in outdir; do not run NGSadmix")
    parser.add_argument("-N","--ngs-only", action="store_true", help="Only run NGSadmix; skip Evanno/ΔK and plotting")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Run NGSadmix unless evanno-only
    if not args.evanno_only:
        df = run_ngsadmix(args.likes, args.min_k, args.max_k, args.threads, args.reps, args.outdir)
    else:
        df = load_or_parse_likelihoods(args.outdir, args.likes, args.min_k, args.max_k, args.reps)

    # Mode checks & Evanno
    if args.evanno_only and args.ngs_only:
        raise SystemExit("Cannot use --evanno-only and --ngs-only together. Choose one.")

    if args.ngs_only:
        print("[mode] NGSadmix-only: will not compute Evanno/ΔK.")
        evanno_best_ks = []
    else:
        evanno_best_ks, _ = evanno_method(df, args.outdir, multispecies=args.multispecies)

    # Write bestK.txt
    bestk_path = os.path.join(args.outdir, "bestK.txt")
    with open(bestk_path, "w") as f:
        if evanno_best_ks:
            f.write(f"Evanno:{','.join(map(str, evanno_best_ks))}\n")
        else:
            f.write("Evanno:\n")
    print(f"[write] {bestk_path}")

    # Copy best artifacts (.qopt/.fopt.gz/.filter)
    if evanno_best_ks:
        extract_best_q_matrices(args.outdir, df, evanno_best_ks)

    # Optional plotting (skip if NGS-only)
    if args.plot and (not args.ngs_only):
        delta_tsv = os.path.join(args.outdir, "evanno_deltaK.tsv")
        out_png = args.plot_out if args.plot_out else os.path.join(args.outdir, "evanno_deltaK.png")
        try:
            write_and_run_evanno_plot(args.rscript, delta_tsv, out_png, width=args.plot_width, height=args.plot_height)
        except Exception as e:
            print(f"[plot] Skipping plot due to error: {e}")

    print(f"[done] Outputs in: {args.outdir}/")

if __name__ == "__main__":
    main()
