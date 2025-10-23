#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estimate_best_k_NGSadmix.evanno.py
Run NGSadmix across K values, collect log-likelihoods, and compute Best K using
the Evanno method exactly mirroring the logic in the provided R implementation.

Author: Taylor Hains
Date: 2025-10-22

Notes:
- Evanno requires SEQUENTIAL K values (no gaps) and the SAME number of replicates per K.
- For any K where the within-K standard deviation of log-likelihoods is zero,
  Evanno statistics cannot be computed (will raise RuntimeError), matching the R behavior.
- Outputs:
    <outdir>/likelihoods.tsv                     (K, rep, logL)
    <outdir>/evanno.summary.tsv                  (K, variable, value, sd)  # L(K), L'(K), L''(K), delta K
    <outdir>/evanno.deltaK.tsv                   (K, DeltaK)               # convenience copy of delta K
    <outdir>/loglikelihood_plot.png              (mean L(K) vs K)
    <outdir>/evanno_plot.png                     (Delta K vs K)
"""
import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_ngsadmix(beagle_file: str, max_k: int, threads: int, replicates: int, outdir: str) -> pd.DataFrame:
    basename = os.path.basename(beagle_file).replace('.beagle.gz', '')
    os.makedirs(outdir, exist_ok=True)
    log_likelihoods = []

    for k in range(1, max_k + 1):
        for rep in range(1, replicates + 1):
            outprefix = f"{outdir}/{basename}.K{k}.rep{rep}"
            log_file = f"{outprefix}.log"
            print(f"[NGSadmix] Running K={k} Rep={rep}")
            with open(log_file, 'w') as lf:
                subprocess.run(
                    ["NGSadmix", "-likes", beagle_file, "-K", str(k), "-P", str(threads), "-o", outprefix],
                    stderr=lf, stdout=subprocess.DEVNULL, check=False
                )

            # Extract final log-likelihood from the log
            try:
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                like_lines = [l for l in lines if 'like=' in l]
                if like_lines:
                    # example fragment: "... like= -12345.6789 ..."
                    last = like_lines[-1]
                    val_str = last.split('like=')[1].split()[0]
                    logl = float(val_str)
                    log_likelihoods.append((k, rep, logl))
                else:
                    print(f"[warn] No 'like=' line found in {log_file}; skipping entry.", file=sys.stderr)
            except Exception as e:
                print(f"[warn] Failed to parse log-likelihood from {log_file}: {e}", file=sys.stderr)

    if not log_likelihoods:
        raise RuntimeError("No log-likelihoods were parsed. Check that NGSadmix ran and produced logs.")

    df = pd.DataFrame(log_likelihoods, columns=["K", "rep", "logL"]).sort_values(["K", "rep"])
    df.to_csv(f"{outdir}/likelihoods.tsv", sep="\t", index=False)
    return df


def _check_evanno_requirements(df: pd.DataFrame):
    # Require sequential Ks and equal number of replicates for each K
    counts = df.groupby('K')['rep'].nunique().to_dict()
    ks_sorted = sorted(counts.keys())
    sequential = all((ks_sorted[i+1] - ks_sorted[i] == 1) for i in range(len(ks_sorted)-1))
    equal_reps = len(set(counts.values())) == 1
    if not sequential or not equal_reps:
        raise RuntimeError("K values are not sequential or there is an uneven number of runs per K. "
                           "Evanno statistics cannot be computed.")


def evanno_method_strict(df: pd.DataFrame, outdir: str) -> pd.DataFrame:
    """
    Mirror the R bestK_evanno() logic:
      For each K in [minK..maxK], compute:
        L(K)      = mean of logL per K
        L'(K)     = mean of (logL_K - logL_{K-1}) [undefined/NA at bounds]
        L''(K)    = mean of abs(logL_{K+1} - 2*logL_K + logL_{K-1}) [NA at bounds]
        delta K   = mean(abs(L''(K))) / sd(L(K)), requires sd(L(K)) > 0; otherwise error
    """
    _check_evanno_requirements(df)

    # Split by K
    split = {k: g['logL'].reset_index(drop=True) for k, g in df.groupby('K')}
    ks = sorted(split.keys())
    kmin, kmax = ks[0], ks[-1]

    rows = []
    for k in ks:
        LK = split[k]  # vector of logL for this K
        mean_LK = float(LK.mean())
        sd_LK = float(LK.std(ddof=1)) if len(LK) > 1 else 0.0

        if k > kmin and k < kmax:
            dK_vec = split[k] - split[k-1]
            ddK_vec = split[k+1] - 2*split[k] + split[k-1]

            mean_dK = float(dK_vec.mean())
            # R version uses mean(abs(ddK))
            mean_abs_ddK = float(np.abs(ddK_vec).mean())

            if sd_LK == 0.0:
                raise RuntimeError("No deviation between runs of the same K. Evanno statistics cannot be computed.")

            deltaK = mean_abs_ddK / sd_LK
            sd_dK = float(dK_vec.std(ddof=1)) if len(dK_vec) > 1 else np.nan
            sd_ddK = float(ddK_vec.std(ddof=1)) if len(ddK_vec) > 1 else np.nan
        else:
            mean_dK = np.nan
            mean_abs_ddK = np.nan
            deltaK = np.nan
            sd_dK = np.nan
            sd_ddK = np.nan

        rows.extend([
            (k, "L(K)", mean_LK, sd_LK),
            (k, "L'(K)", mean_dK, sd_dK),
            (k, "L''(K)", mean_abs_ddK, sd_ddK),
            (k, "delta K", deltaK, np.nan),
        ])

    summary = pd.DataFrame(rows, columns=["K", "variable", "value", "sd"])
    summary.to_csv(f"{outdir}/evanno.summary.tsv", sep="\t", index=False)

    # Convenience DeltaK table (drop NaNs)
    delta_df = summary[(summary["variable"] == "delta K") & (~summary["value"].isna())][["K", "value"]].copy()
    delta_df = delta_df.rename(columns={"value": "DeltaK"}).sort_values("K")
    delta_df.to_csv(f"{outdir}/evanno.deltaK.tsv", sep="\t", index=False)

    # Best K based on max DeltaK
    if not delta_df.empty:
        best_k = int(delta_df.loc[delta_df["DeltaK"].idxmax(), "K"])
    else:
        best_k = None

    return summary, delta_df, best_k


def plot_results(df: pd.DataFrame, evanno_summary: pd.DataFrame, delta_df: pd.DataFrame, outdir: str):
    # Plot mean L(K) vs K
    means = df.groupby("K")["logL"].mean()
    plt.figure()
    means.plot(marker='o', title="Mean Log-Likelihood vs K")
    plt.xlabel("K")
    plt.ylabel("Mean log-likelihood")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{outdir}/loglikelihood_plot.png", dpi=200)
    plt.close()

    # Plot Delta K vs K
    if not delta_df.empty:
        plt.figure()
        delta_df.set_index("K")["DeltaK"].plot(marker='o', title="Evanno Delta K")
        plt.xlabel("K")
        plt.ylabel("Delta K")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{outdir}/evanno_plot.png", dpi=200)
        plt.close()


def main():
    ap = argparse.ArgumentParser(description="Run NGSadmix across K values and estimate best K using the Evanno method (strict, R-matched).")
    ap.add_argument("--beagle", required=True, help="Input .beagle.gz file")
    ap.add_argument("--maxK", type=int, required=True, help="Maximum K to test (starts at K=1)")
    ap.add_argument("--threads", type=int, default=4, help="Threads for NGSadmix")
    ap.add_argument("--reps", type=int, default=3, help="Number of replicates per K")
    ap.add_argument("--outdir", default="NGSadmix_runs", help="Output directory")
    args = ap.parse_args()

    df = run_ngsadmix(args.beagle, args.maxK, args.threads, args.reps, args.outdir)

    # Compute Evanno (strict checks)
    try:
        evanno_summary, delta_df, best_k = evanno_method_strict(df, args.outdir)
    except RuntimeError as e:
        # Mirror strict R behavior by failing loudly with a clear message
        print(f"[ERROR] {e}", file=sys.stderr)
        print(f"Likelihoods saved to: {args.outdir}/likelihoods.tsv")
        sys.exit(1)

    # Plots
    plot_results(df, evanno_summary, delta_df, args.outdir)

    # Report
    if best_k is not None:
        print(f"\nBest K (Evanno): {best_k}")
        print("Files written:")
        print(f"  - {args.outdir}/likelihoods.tsv")
        print(f"  - {args.outdir}/evanno.summary.tsv")
        print(f"  - {args.outdir}/evanno.deltaK.tsv")
        print(f"  - {args.outdir}/loglikelihood_plot.png")
        print(f"  - {args.outdir}/evanno_plot.png")
    else:
        print("Delta K not computable (likely due to boundary-only values). Check inputs.")


if __name__ == "__main__":
    main()
