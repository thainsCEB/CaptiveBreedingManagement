#!/usr/bin/env python3

import os
import subprocess
import argparse
import re
import pandas as pd

def run_admixture(plink_prefix, max_k, threads, outdir):
    os.makedirs(outdir, exist_ok=True)
    cv_errors = []

    for k in range(1, max_k + 1):
        print(f"Running ADMIXTURE for K={k}")
        out_log = os.path.join(outdir, f"{os.path.basename(plink_prefix)}.K{k}.log")
        with open(out_log, "w") as logf:
            subprocess.run(
                ["admixture", "--cv", f"-j{threads}", f"{plink_prefix}.bed", str(k)],
                stdout=logf,
                stderr=subprocess.STDOUT,
                cwd=outdir,
            )

        # Parse the CV error
        with open(out_log, "r") as logf:
            for line in logf:
                match = re.search(r"CV error.*:\s+([\d\.]+)", line)
                if match:
                    cv_error = float(match.group(1))
                    cv_errors.append((k, cv_error))
                    break

    # Convert to DataFrame
    df = pd.DataFrame(cv_errors, columns=["K", "CV_Error"])
    df.sort_values("CV_Error", inplace=True)
    df.to_csv(os.path.join(outdir, "cv_errors.tsv"), sep="\t", index=False)

    return df

def main():
    parser = argparse.ArgumentParser(description="Estimate best K using ADMIXTURE and CV error.")
    parser.add_argument("-p", "--plink", required=True, help="Prefix of PLINK files (no .bed/.bim/.fam extension)")
    parser.add_argument("-k", "--max-k", type=int, required=True, help="Maximum K to test (starts at 1)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("-o", "--outdir", default="admixture_out", help="Output directory")
    parser.add_argument("--species-mode", choices=["single", "multi"], default="single",
                        help="Species mode: 'single' for best K, 'multi' for top 5 Ks")

    args = parser.parse_args()

    df = run_admixture(args.plink, args.max_k, args.threads, args.outdir)

    print("\nCross-validation errors:")
    print(df.to_string(index=False))

    if args.species_mode == "single":
        best_k = df.iloc[0]["K"]
        print(f"\nBest K (single species mode): K={int(best_k)}")
    else:
        print(f"\nTop 5 K values (multi-species mode):")
        print(df.head(5).to_string(index=False))

if __name__ == "__main__":
    main()
