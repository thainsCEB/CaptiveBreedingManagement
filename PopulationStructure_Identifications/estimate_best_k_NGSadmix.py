#!/usr/bin/env python3

import os
import subprocess
import argparse
import pandas as pd
import numpy as np
import shutil


def run_ngsadmix(beagle_file, max_k, threads, replicates, outdir):
    basename = os.path.basename(beagle_file).replace('.beagle.gz', '')
    os.makedirs(outdir, exist_ok=True)
    log_likelihoods = []

    for k in range(1, max_k + 1):
        for rep in range(1, replicates + 1):
            outprefix = f"{outdir}/{basename}.K{k}.rep{rep}"
            log_file = f"{outprefix}.log"
            print(f"Running K={k} Rep={rep}")
            subprocess.run([
                "NGSadmix", "-likes", beagle_file,
                "-K", str(k), "-P", str(threads),
                "-printInfo", "1",  # Added here
                "-o", outprefix
            ], stderr=open(log_file, 'w'))

            with open(log_file, 'r') as f:
                lines = f.readlines()
            like_lines = [l for l in lines if 'like=' in l]
            if like_lines:
                logl = float(like_lines[-1].split('like=')[1].split()[0])
                log_likelihoods.append((k, rep, logl, beagle_file))

    df = pd.DataFrame(log_likelihoods, columns=["K", "rep", "logL", "beagle_file"])
    df.to_csv(f"{outdir}/likelihoods.tsv", sep="\t", index=False)
    return df


def evanno_method(df, outdir, is_multispecies):
    means = df.groupby("K")["logL"].mean()
    sds = df.groupby("K")["logL"].std()
    delta_k = []

    for k in range(2, means.index.max()):
        l_km1 = means.get(k - 1, np.nan)
        l_k = means.get(k, np.nan)
        l_kp1 = means.get(k + 1, np.nan)
        sd_k = sds.get(k, np.nan)
        if pd.notna(l_km1) and pd.notna(l_k) and pd.notna(l_kp1) and sd_k > 0:
            delta = abs(l_kp1 - 2 * l_k + l_km1) / sd_k
            delta_k.append((k, delta))

    delta_df = pd.DataFrame(delta_k, columns=["K", "DeltaK"])
    delta_df.to_csv(f"{outdir}/evanno_deltaK.tsv", sep="\t", index=False)

    if not delta_df.empty:
        if is_multispecies:
            top5 = delta_df.sort_values(by="DeltaK", ascending=False).head(5)
            print("\nTop 5 K values based on DeltaK (multi-species mode):")
            print(top5)
            best_ks = top5["K"].tolist()
        else:
            best_k = delta_df.loc[delta_df["DeltaK"].idxmax(), "K"]
            print(f"\nBest K based on Evanno method (single-species mode): {best_k}")
            best_ks = [best_k]
    else:
        best_ks = []
        print("Insufficient data to compute Evanno DeltaK.")

    return means, delta_df, best_ks


def extract_best_q_matrices(outdir, df_likelihoods, best_ks, output_prefix="bestQ"):
    saved_files = []
    for k in best_ks:
        subset = df_likelihoods[df_likelihoods["K"] == k]
        best_rep = subset.loc[subset["logL"].idxmax(), "rep"]
        basename = os.path.basename(df_likelihoods["beagle_file"].iloc[0]).replace('.beagle.gz', '')
        source_q = os.path.join(outdir, f"{basename}.K{k}.rep{int(best_rep)}.qopt")
        dest_q = os.path.join(outdir, f"{output_prefix}.K{k}.qopt")

        try:
            shutil.copyfile(source_q, dest_q)
            saved_files.append(dest_q)
            print(f"Saved best Q matrix for K={k} (rep={int(best_rep)}): {dest_q}")
        except FileNotFoundError:
            print(f"Warning: Q matrix not found for K={k}, replicate={int(best_rep)}")

    return saved_files


def main():
    parser = argparse.ArgumentParser(description="Run NGSadmix across K values and estimate best K using Evanno method.")
    parser.add_argument("--beagle", required=True, help="Input .beagle.gz file")
    parser.add_argument("--maxK", type=int, required=True, help="Maximum K to test")
    parser.add_argument("--threads", type=int, default=4, help="Threads for NGSadmix")
    parser.add_argument("--reps", type=int, default=3, help="Number of replicates per K")
    parser.add_argument("--outdir", default="NGSadmix_runs", help="Output directory")
    parser.add_argument("--multispecies", action="store_true", help="Use multi-species mode to report 4th highest DeltaK")

    args = parser.parse_args()
    df = run_ngsadmix(args.beagle, args.maxK, args.threads, args.reps, args.outdir)
    means, delta_df, best_ks = evanno_method(df, args.outdir, args.multispecies)
    q_files = extract_best_q_matrices(args.outdir, df, best_ks)

    if best_ks:
        print(f"\nSuggested K(s): {', '.join(map(str, best_ks))}")
    print(f"Results saved in: {args.outdir}/")


if __name__ == "__main__":
    main()
