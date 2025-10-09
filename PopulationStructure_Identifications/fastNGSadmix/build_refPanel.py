#!/usr/bin/env python3

import argparse
import subprocess
import os
import gzip
import pandas as pd


def _strip_known_suffixes(path: str) -> str:
    """
    Return the path with known suffixes removed to get the common prefix.
    Known suffixes include: .beagle.gz, .beagle, .fopt.gz, .qopt, .beagle.filter, .mafs.gz, .mafs
    """
    for suf in [".beagle.gz", ".beagle", ".fopt.gz", ".qopt", ".beagle.filter", ".mafs.gz", ".mafs"]:
        if path.endswith(suf):
            return path[:-len(suf)]
    # If none matched, remove trailing extensions generically
    for suf in [".gz", ".txt", ".tsv"]:
        if path.endswith(suf):
            return path[:-len(suf)]
    return path

def _first_existing(candidates):
    for c in candidates:
        if os.path.exists(c):
            return c
    return None

def infer_ngsadmix_paths(args):
    """
    Infer paths for beagle, fopt, qopt, filter, mafs from either --in-prefix or any one of the individual files.
    Respects already-provided args; only fills missing ones.
    Returns a dict with keys: beagle, fopt, qopt, filter, mafs
    """
    prefix = None
    # Priority: in-prefix, else derive from any provided file
    if getattr(args, "in_prefix", None):
        prefix = args.in_prefix
    else:
        for key in ["beagle","fopt","qopt","filter","mafs"]:
            val = getattr(args, key, None)
            if val:
                prefix = _strip_known_suffixes(val)
                break

    if not prefix:
        raise SystemExit("[ERROR] Provide --in-prefix OR one of --beagle/--fopt/--qopt/--filter/--mafs to infer the rest.")

    paths = {
        "beagle": getattr(args, "beagle", None),
        "fopt": getattr(args, "fopt", None),
        "qopt": getattr(args, "qopt", None),
        "filter": getattr(args, "filter", None),
        "mafs": getattr(args, "mafs", None),
    }

    # Propose candidates from prefix
    beagle_cands = [prefix + ".beagle.gz", prefix + ".beagle"]
    fopt_cands   = [prefix + ".fopt.gz", prefix + ".fopt"]
    qopt_cands   = [prefix + ".qopt"]
    filter_cands = [prefix + ".beagle.filter", prefix + ".filter"]
    mafs_cands   = [prefix + ".mafs.gz", prefix + ".mafs"]

    if paths["beagle"] is None:
        paths["beagle"] = _first_existing(beagle_cands) or beagle_cands[0]
    if paths["fopt"] is None:
        paths["fopt"] = _first_existing(fopt_cands) or fopt_cands[0]
    if paths["qopt"] is None:
        paths["qopt"] = _first_existing(qopt_cands) or qopt_cands[0]
    if paths["filter"] is None:
        paths["filter"] = _first_existing(filter_cands) or filter_cands[0]
    if paths["mafs"] is None:
        paths["mafs"] = _first_existing(mafs_cands) or mafs_cands[0]

    return paths, prefix
def run(cmd, desc):
    print(f"[RUNNING] {desc}")
    subprocess.run(cmd, shell=True, check=True)

def decompress_gz(gz_path, out_path):
    print(f"[INFO] Decompressing {gz_path} to {out_path}")
    with gzip.open(gz_path, 'rt') as fin, open(out_path, 'w') as fout:
        fout.writelines(fin)

def parse_ref_beagle(tmp_beagle, filter_file, output_file="tmp.ref") -> bool:
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


def generate_sites_and_bed(refpanel_file, prefix, mafs_file=None):
    """
    Write:
      - refPanel_{prefix}.refPanel.sites => STRICTLY from .mafs(.gz) columns 1-4:
            chr  pos  major  minor
      - refPanel_{prefix}.refPanel.bed   => from ref panel chr/pos: chr  start  end (0-based)
    If mafs_file is missing/unreadable, raise an error (to ensure alleles come from MAFS only).
    """
    import os, gzip
    import pandas as pd

    print(f"[INFO] Generating .sites and .bed files from {refpanel_file}")
    if not mafs_file or not os.path.exists(mafs_file):
        raise FileNotFoundError("[ERROR] --mafs not provided or not found; the .sites must be built from MAFS.")

    # Read only columns 1-4 (1-based) = 0-3 (0-based): chromo, position, major, minor
    if str(mafs_file).endswith(".gz"):
        with gzip.open(mafs_file, "rt") as fin:
            mafs_df = pd.read_csv(fin, sep=r"\s+|\t", engine="python", header=0, usecols=[0,1,2,3])
    else:
        mafs_df = pd.read_csv(mafs_file, sep=r"\s+|\t", engine="python", header=0, usecols=[0,1,2,3])
    mafs_df.columns = ["chr","pos","major","minor"]

    sites_file = f"refPanel_{prefix}.refPanel.sites"
    bed_file = f"refPanel_{prefix}.refPanel.bed"

    # Write sites strictly from MAFS
    mafs_df.to_csv(sites_file, sep="\t", index=False, header=False)
    print(f"[INFO] Wrote sites from MAFS (chr pos major minor) to {sites_file}")

    # BED from the ref panel table
    df = pd.read_csv(refpanel_file, sep=r"\s+|\t", engine="python", dtype=str)
    cols_lower = {c.lower(): c for c in df.columns}
    if not all(k in cols_lower for k in ("chr","pos")):
        raise ValueError("[ERROR] ref panel table must contain 'chr' and 'pos' columns to make BED.")
    chr_col = cols_lower["chr"]; pos_col = cols_lower["pos"]
    bed_df = pd.DataFrame({
        "chr": df[chr_col],
        "start": df[pos_col].astype(int) - 1,
        "end": df[pos_col].astype(int)
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

    parser.add_argument("--beagle", help="Input .beagle(.gz) file (NGSadmix mode)")
    parser.add_argument("--in-prefix", help="Prefix to infer --beagle, --fopt, --qopt, --filter, and --mafs (NGSadmix mode).")
    parser.add_argument("--filter", help=".beagle.filter file (NGSadmix mode)")
    parser.add_argument("--fopt", help=".fopt.gz file (NGSadmix mode)")
    parser.add_argument("--qopt", help=".qopt file (NGSadmix mode)")

    parser.add_argument("--bim", help=".bim file (ADMIXTURE mode)")
    parser.add_argument("--pfile", help=".P file (ADMIXTURE mode)")
    parser.add_argument("--qfile", help=".Q file (ADMIXTURE mode)")

    parser.add_argument("--prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    if args.mode == "ngsadmix":
        # Infer all required file paths from --in-prefix or any provided file
        paths, inferred_prefix = infer_ngsadmix_paths(args)
        args.beagle = paths["beagle"]
        args.fopt   = paths["fopt"]
        args.qopt   = paths["qopt"]
        args.filter = paths["filter"]
        args.mafs   = paths["mafs"]

        # Now validate that the inferred/provided files exist where needed
        if not all([args.beagle, args.fopt, args.qopt, args.filter]):
            raise SystemExit("[ERROR] Could not infer all required files. Provide --in-prefix or explicit paths.")

        tmp_beagle = "tmp.beagle"
        tmp_fopt = "tmp.fopt"
        decompress_gz(args.beagle, tmp_beagle)
        decompress_gz(args.fopt, tmp_fopt)

        out_refpanel = f"refPanel_{args.prefix}.refPanel.txt"
        # Infer mafs_path from fopt/beagle prefixes if not provided
        def _strip_known_suffixes(path: str) -> str:
            for suf in [".beagle.gz",".beagle",".fopt.gz",".fopt",".qopt",".beagle.filter",".filter",".mafs.gz",".mafs"]:
                if path.endswith(suf):
                    return path[:-len(suf)]
            return path
        def _first_existing(cands):
            import os
            for c in cands:
                if os.path.exists(c):
                    return c
            return None
        _pref_fopt = _strip_known_suffixes(args.fopt) if getattr(args, "fopt", None) else None
        _pref_beag = _strip_known_suffixes(args.beagle) if getattr(args, "beagle", None) else None
        mafs_path = getattr(args, "mafs", None)
        if not mafs_path:
            if _pref_fopt:
                mafs_path = _first_existing([_pref_fopt + ".mafs.gz", _pref_fopt + ".mafs"])
            if not mafs_path and _pref_beag:
                mafs_path = _first_existing([_pref_beag + ".mafs.gz", _pref_beag + ".mafs"])
        if not mafs_path:
            raise SystemExit("[ERROR] Could not infer .mafs(.gz). Provide --mafs or place it next to your .fopt/.beagle with the same prefix.")
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

        generate_sites_and_bed(out_refpanel, args.prefix, mafs_path)
        print("[DONE] NGSadmix reference panel and cluster summary built.")

    elif args.mode == "admixture":
        if not all([args.bim, args.pfile, args.qfile]):
            parser.error("ADMIXTURE mode requires --bim, --pfile, and --qfile")
        build_from_admixture(args.bim, args.pfile, args.qfile, args.prefix, args.K)

if __name__ == "__main__":
    main()
