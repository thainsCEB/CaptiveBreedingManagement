#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Title: Run NGSadmix across K, estimate best K (Evanno), and extract best Q
# Description:
#   - Runs NGSadmix for K = [minK..maxK] with N replicates each (always with -printInfo 1).
#   - Parses log-likelihoods, computes Evanno DeltaK, and writes evanno_deltaK.tsv.
#   - Writes bestK.txt (Evanno result) and copies best .qopt for chosen K(s).
# Author: Taylor Hains
# Date: 2025-10-17
# -----------------------------------------------------------------------------

import os
import re
import shutil
import argparse
import subprocess
import pandas as pd
import numpy as np

def run_cmd(cmd, log_path=None):
    """Run a command; if log_path is provided, stderr goes to that file."""
    if log_path:
        with open(log_path, "w") as lf:
            return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=lf,
                                  check=False, text=True)
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          check=False, text=True)

def run_ngsadmix(beagle_file, min_k, max_k, threads, replicates, outdir):
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
                "-printInfo", "1",
                "-o", outprefix
            ]
            run_cmd(cmd, log_path=log_file)

            # Parse the final likelihood from the log (last line containing "like=")
            logl = None
            try:
                with open(log_file) as fh:
                    lines = fh.readlines()
                like_lines = [ln for ln in lines if "like=" in ln]
                if like_lines:
                    m = re.search(r'like=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', like_lines[-1])
                    if m:
                        logl = float(m.group(1))
            except Exception as e:
                print(f"[warn] Could not parse likelihood from {log_file}: {e}")

            rows.append((k, rep, logl, beagle_file))

    df = pd.DataFrame(rows, columns=["K", "rep", "logL", "beagle_file"])
    df.to_csv(f"{outdir}/likelihoods.tsv", sep="\t", index=False)
    return df

def evanno_method(df, outdir, multispecies=False):
    means = df.groupby("K")["logL"].mean()
    sds   = df.groupby("K")["logL"].std()

    if means.index.size == 0:
        print("[evanno] No likelihoods available.")
        pd.DataFrame(columns=["K","DeltaK"]).to_csv(f"{outdir}/evanno_deltaK.tsv", sep="\t", index=False)
        return [], pd.DataFrame(columns=["K","DeltaK"])

    delta_k = []
    k_min, k_max = int(min(means.index)), int(max(means.index))
    for k in range(k_min + 1, k_max):
        l_km1 = means.get(k - 1, np.nan)
        l_k   = means.get(k,     np.nan)
        l_kp1 = means.get(k + 1, np.nan)
        sd_k  = sds.get(k,       np.nan)
        if np.all(pd.notna([l_km1, l_k, l_kp1, sd_k])) and sd_k > 0:
            delta_k.append((k, abs(l_kp1 - 2*l_k + l_km1) / sd_k))

    delta_df = pd.DataFrame(delta_k, columns=["K","DeltaK"])
    delta_df.to_csv(f"{outdir}/evanno_deltaK.tsv", sep="\t", index=False)

    if delta_df.empty:
        print("[evanno] Insufficient data for DeltaK (need at least K-1, K, K+1 with nonzero SD).")
        return [], delta_df

    if multispecies:
        best_ks = delta_df.sort_values("DeltaK", ascending=False).head(5)["K"].astype(int).tolist()
        print(f"[evanno] Top ΔK (multi): {best_ks}")
    else:
        best_ks = [int(delta_df.loc[delta_df["DeltaK"].idxmax(), "K"])]
        print(f"[evanno] Best K = {best_ks[0]}")
    return best_ks, delta_df

def extract_best_artifacts(outdir, df_likelihoods, best_ks, output_prefix="bestQ"):
    """
    For each K in best_ks, find the replicate with the highest logL and copy:
      - {basename}.K{K}.rep{rep}.qopt   -> {output_prefix}.K{K}.qopt
      - {basename}.K{K}.rep{rep}.fopt.gz-> {output_prefix}.K{K}.fopt.gz
      - {basename}.K{K}.rep{rep}.filter -> {output_prefix}.K{K}.filter
    Returns a list of destination file paths successfully written.
    """
    saved = []
    if df_likelihoods.empty or not best_ks:
        return saved

    basename = os.path.basename(df_likelihoods["beagle_file"].iloc[0]).replace(".beagle.gz", "")
    for k in sorted(set(best_ks)):
        subset = df_likelihoods[(df_likelihoods["K"] == k) & df_likelihoods["logL"].notna()]
        if subset.empty:
            print(f"[save] No logL for K={k}; skip copying artifacts")
            continue

        rep = int(subset.loc[subset["logL"].idxmax(), "rep"])
        base = os.path.join(outdir, f"{basename}.K{k}.rep{rep}")
        srcs = {
            ".qopt":     f"{base}.qopt",
            ".fopt.gz":  f"{base}.fopt.gz",
            ".filter":   f"{base}.filter",
        }
        dsts = {
            ".qopt":     os.path.join(outdir, f"{output_prefix}.K{k}.qopt"),
            ".fopt.gz":  os.path.join(outdir, f"{output_prefix}.K{k}.fopt.gz"),
            ".filter":   os.path.join(outdir, f"{output_prefix}.K{k}.filter"),
        }

        for ext, src in srcs.items():
            dst = dsts[ext]
            try:
                if os.path.exists(src):
                    shutil.copyfile(src, dst)
                    print(f"[save] Best {ext} for K={k} (rep={rep}) → {dst}")
                    saved.append(dst)
                else:
                    print(f"[save] Missing {src} (skipping)")
            except Exception as e:
                print(f"[warn] Failed to copy {src} → {dst}: {e}")

    return saved


def main():
    ap = argparse.ArgumentParser(
        description="NGSadmix K scan with Evanno best-K and best .qopt extraction (no ADMIXTURE)."
    )
    # Short-hand flags
    ap.add_argument("-b","--beagle", required=True, help="Input .beagle.gz for NGSadmix")
    ap.add_argument("-i","--minK", type=int, default=1, help="Minimum K (default: 1)")
    ap.add_argument("-k","--maxK", type=int, required=True, help="Maximum K")
    ap.add_argument("-t","--threads", type=int, default=4, help="Threads (default: 4)")
    ap.add_argument("-r","--reps", type=int, default=3, help="Replicates per K for NGSadmix (default: 3)")
    ap.add_argument("-o","--outdir", default="NGSadmix_runs", help="Output directory (default: NGSadmix_runs)")
    ap.add_argument("-m","--multispecies", action="store_true",
                    help="Evanno multi-species mode: report top ΔK values (top 5)")

    args = ap.parse_args()

    if args.minK < 1 or args.maxK < args.minK:
        raise SystemExit(f"[error] Invalid K range: minK={args.minK}, maxK={args.maxK}")

    os.makedirs(args.outdir, exist_ok=True)

    # NGSadmix + Evanno
    df = run_ngsadmix(args.beagle, args.minK, args.maxK, args.threads, args.reps, args.outdir)
    evanno_best_ks, _ = evanno_method(df, args.outdir, multispecies=args.multispecies)

    # Save summary bestK.txt (Evanno only)
    bestk_path = os.path.join(args.outdir, "bestK.txt")
    with open(bestk_path, "w") as f:
        if evanno_best_ks:
            f.write(f"Evanno:{','.join(map(str, evanno_best_ks))}\n")
        else:
            f.write("Evanno:\n")
    print(f"[write] {bestk_path}")

    # Copy best .qopt for any chosen K from Evanno
    extract_best_q_matrices(args.outdir, df, evanno_best_ks)

    print(f"[done] Outputs in: {args.outdir}/")

if __name__ == "__main__":
    main()
