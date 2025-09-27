#!/usr/bin/env python3

import argparse
import subprocess
import os
import gzip
import pandas as pd

def run(cmd, desc):
    print(f"[RUNNING] {desc}")
    subprocess.run(cmd, shell=True, check=True)

def decompress_gz(gz_path, out_path):
    print(f"[INFO] Decompressing {gz_path} to {out_path}")
    with gzip.open(gz_path, 'rt') as fin, open(out_path, 'w') as fout:
        fout.writelines(fin)

def parse_ref_beagle(tmp_beagle, filter_file, output_file="tmp.ref"):
    print("[INFO] Filtering and building tmp.ref")
    tmp = pd.read_csv(tmp_beagle, sep="\t", header=None, dtype=str)
    filter_df = pd.read_csv(filter_file, sep="\t", header=0, dtype=str)
    tmp = tmp[tmp.iloc[:, 0].isin(filter_df.iloc[:, 0])].copy()

    tmp["chr"] = tmp.iloc[:, 0].apply(lambda x: x.split("_")[0])
    tmp["pos"] = tmp.iloc[:, 0].apply(lambda x: x.split("_")[1])
    tmp["name"] = tmp.iloc[:, 0]
    tmp["A0_freq"] = tmp.iloc[:, 2]
    tmp["A1"] = tmp.iloc[:, 1]
    tmp["id"] = tmp.iloc[:, 0]

    out_df = tmp[["id", "chr", "pos", "name", "A0_freq", "A1"]]
    out_df.to_csv(output_file, sep="\t", index=False, header=False)

def generate_sites_and_bed(refpanel_file, prefix):
    print(f"[INFO] Generating .sites and .bed files from {refpanel_file}")
    df = pd.read_csv(refpanel_file, sep="\s+|	", engine="python")
    sites_file = f"refPanel_{prefix}.refPanel.sites"
    bed_file = f"refPanel_{prefix}.refPanel.bed"

    df[["chr", "pos"]].to_csv(sites_file, sep="\t", index=False, header=False)

    bed_df = pd.DataFrame({
        "chr": df["chr"],
        "start": df["pos"].astype(int) - 1,
        "end": df["pos"].astype(int)
    })
    bed_df.to_csv(bed_file, sep="\t", index=False, header=False)

def build_from_admixture(bim_file, p_file, q_file, prefix, K):
    print("[INFO] Building reference panel and ancestry stats from ADMIXTURE input")
    ref_out = f"refPanel_{prefix}.P.txt"
    nind_out = f"nInd_{prefix}.Q.txt"

    with open(ref_out, "w") as f:
        f.write("id chr pos name A0_freq A1 " + " ".join([f"K{i+1}" for i in range(K)]) + "\n")

    tmp_txt = f"{prefix}.tmp.txt"
    with open(tmp_txt, "w") as tmp_out, open(bim_file, "r") as bim:
        for line in bim:
            fields = line.strip().split()
            if len(fields) >= 6:
                id_ = f"{fields[0]}_{fields[3]}"
                tmp_out.write(f"{id_} {fields[0]} {fields[3]} {fields[1]} {fields[5]} {fields[4]}\n")

    run(f"paste -d' ' {tmp_txt} {p_file} >> {ref_out}", "Merging .bim and .P file")
    os.remove(tmp_txt)

    with open(nind_out, "w") as f:
        qopt = pd.read_csv(q_file, sep=" ", header=None)
        if qopt.shape[1] < K:
            raise ValueError(f"[ERROR] Found only {qopt.shape[1]} columns in qfile, but --K={K} was requested.")
        f.write(" ".join([f"K{i+1}" for i in range(K)]) + "\n")
        f.write(" ".join(qopt.iloc[:, :K].sum().astype(str).tolist()) + "\n")

    generate_sites_and_bed(ref_out, prefix)
    print("[DONE] ADMIXTURE ref panel and cluster summary built.")

def main():
    parser = argparse.ArgumentParser(description="Build reference panel from NGSadmix or ADMIXTURE results")
    parser.add_argument("--mode", choices=["ngsadmix", "admixture"], default="ngsadmix", help="Choose input mode")
    parser.add_argument("--K", type=int, required=True, help="Number of ancestry clusters")

    parser.add_argument("--beagle", help="Input .beagle.gz file (NGSadmix mode)")
    parser.add_argument("--filter", help=".beagle.filter file (NGSadmix mode)")
    parser.add_argument("--fopt", help=".fopt.gz file (NGSadmix mode)")
    parser.add_argument("--qopt", help=".qopt file (NGSadmix mode)")

    parser.add_argument("--bim", help=".bim file (ADMIXTURE mode)")
    parser.add_argument("--pfile", help=".P file (ADMIXTURE mode)")
    parser.add_argument("--qfile", help=".Q file (ADMIXTURE mode)")

    parser.add_argument("--prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    if args.mode == "ngsadmix":
        if not all([args.beagle, args.filter, args.fopt, args.qopt]):
            parser.error("NGSadmix mode requires --beagle, --filter, --fopt, --qopt")

        tmp_beagle = "tmp.beagle"
        tmp_fopt = "tmp.fopt"
        decompress_gz(args.beagle, tmp_beagle)
        decompress_gz(args.fopt, tmp_fopt)

        out_refpanel = f"refPanel_{args.prefix}.refPanel.txt"
        out_nind = f"nInd_{args.prefix}.refPanel.txt"

        with open(out_refpanel, "w") as f:
            f.write("id chr pos name A0_freq A1 " + " ".join([f"K{i+1}" for i in range(args.K)]) + "\n")

        parse_ref_beagle(tmp_beagle, args.filter, "tmp.ref")
        run(f"paste tmp.ref {tmp_fopt} >> {out_refpanel}", "Combining ref and fopt")
        run("rm tmp.ref tmp.beagle tmp.fopt", "Cleaning up")

        with open(out_nind, "w") as f:
            qopt = pd.read_csv(args.qopt, sep=" ", header=None)
            if qopt.shape[1] < args.K:
                raise ValueError(f"[ERROR] Found only {qopt.shape[1]} columns in qopt, but --K={args.K} was requested.")
            f.write(" ".join([f"K{i+1}" for i in range(args.K)]) + "\n")
            f.write(" ".join(qopt.iloc[:, :args.K].sum().astype(str).tolist()) + "\n")

        generate_sites_and_bed(out_refpanel, args.prefix)
        print("[DONE] NGSadmix reference panel and cluster summary built.")

    elif args.mode == "admixture":
        if not all([args.bim, args.pfile, args.qfile]):
            parser.error("ADMIXTURE mode requires --bim, --pfile, and --qfile")
        build_from_admixture(args.bim, args.pfile, args.qfile, args.prefix, args.K)

if __name__ == "__main__":
    main()
