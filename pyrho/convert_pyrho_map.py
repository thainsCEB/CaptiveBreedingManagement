#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import glob
import re

def extract_scaffold_id(filename):
    """Extract scaffold ID as the last '.'-separated chunk before .rmap"""
    return filename.rsplit(".rmap", 1)[0].split(".")[-1]

def is_sex_chromosome(name):
    return bool(re.search(r'\b[XYWZ]\b', name, re.IGNORECASE)) or any(x in name.upper() for x in ["CHRZ", "CHRX", "CHRY", "CHRW", "Z", "W", "X", "Y"])

def convert_rmap_to_df(input_file, scaffold):
    df = pd.read_csv(input_file, sep=r"\s+", header=None, names=["start", "end", "r"])
    df["window_size"] = df["end"] - df["start"]
    df["cM"] = 100 * df["window_size"] * df["r"]
    df["cM_per_Mb"] = df["r"] * 1e8
    df["cumulative_cM"] = df["cM"].cumsum()
    df.insert(0, "scaffold", scaffold)
    return df[["scaffold", "start", "end", "cM_per_Mb", "cumulative_cM", "window_size"]]

def main():
    parser = argparse.ArgumentParser(description="Convert .rmap files to recombination maps with BEDGRAPH format and stats.")
    parser.add_argument("--in-dir", required=True, help="Input directory with .rmap files")
    parser.add_argument("--out-dir", required=True, help="Directory for converted scaffold outputs")
    parser.add_argument("--prefix", required=True, help="Prefix to include in output filenames")
    parser.add_argument("--concat-output", required=True, help="Path to final concatenated output file")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    rmap_files = glob.glob(os.path.join(args.in_dir, "*.rmap"))
    if not rmap_files:
        print("[ERROR] No .rmap files found.")
        return

    all_dfs = []
    sex_chroms = set()
    auto_chroms = set()

    print("\n[INFO] Mean cM/Mb per scaffold:\n")

    for path in rmap_files:
        filename = os.path.basename(path)
        try:
            scaffold = extract_scaffold_id(filename)
            outfile = os.path.join(args.out_dir, f"{scaffold}_{args.prefix}_converted.txt")

            df = convert_rmap_to_df(path, scaffold)
            df.to_csv(outfile, sep="\t", index=False, header=False)
            all_dfs.append(df)

            mean_rate = df["cM_per_Mb"].mean()
            print(f"  {scaffold:<15}  {mean_rate:.4f} cM/Mb")

            if is_sex_chromosome(scaffold):
                sex_chroms.add(scaffold)
            else:
                auto_chroms.add(scaffold)

        except Exception as e:
            print(f"[WARN] Failed to process {filename}: {e}")

    if all_dfs:
        merged = pd.concat(all_dfs, ignore_index=True)
        merged.to_csv(args.concat_output, sep="\t", index=False)

        genome_mean = merged["cM_per_Mb"].mean()

        # Autosomal and sex chromosome separation
        auto_df = merged[merged["scaffold"].isin(auto_chroms)]
        sex_df = merged[merged["scaffold"].isin(sex_chroms)]

        print("\n[INFO] Genome-wide mean recombination rate: {:.4f} cM/Mb".format(genome_mean))
        if not auto_df.empty:
            print("[INFO] Autosomal mean recombination rate:  {:.4f} cM/Mb".format(auto_df["cM_per_Mb"].mean()))
        if not sex_df.empty:
            print("[INFO] Sex chromosome mean recombination rate: {:.4f} cM/Mb".format(sex_df["cM_per_Mb"].mean()))

        print(f"[INFO] Concatenated output written to: {args.concat_output}\n")
    else:
        print("[ERROR] No valid .rmap files processed.")

if __name__ == "__main__":
    main()
