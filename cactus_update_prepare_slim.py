#!/usr/bin/env python3
"""
cactus_update_prepare_exact.py

Create a minimal Cactus update input and run exactly:
  cactus-update-prepare add node <in_hal> <out_prefix>/input.txt \
    --genome <new_name> \
    --outDir <out_prefix>/steps \
    --jobStore <out_prefix>/jobstore \
    --cactus-prepare-options ' --maxCores N --defaultMemory MEM --defaultDisk DISK '

Modes:
  sister   : add NEW as sister to an existing leaf (--target-leaf)
  ancestor : add NEW under a named internal node (--parent-node)
  replace  : replace assembly for an existing leaf (no tree edit)

Outputs:
  <out_prefix>/input.txt        # tree (line 1) + "<new_name>\t<staged_fa>" (line 2)
  <out_prefix>/inputs/          # FASTA symlink
  <out_prefix>/steps/           # prepare's steps dir (created)
  <out_prefix>/jobstore/        # jobstore dir (created)
  <out_prefix>.sh               # emitted pipeline (executable)
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime

LEAF_TOKEN = r"[A-Za-z0-9_.-]+"


# -------- Tree helpers --------

def _strip_lengths(t: str) -> str:
    return re.sub(r":[0-9eE.+-]+", "", t)

def _ensure_semicolon(t: str) -> str:
    t = t.strip()
    return t if t.endswith(";") else t + ";"

def _count_leaves(t: str) -> int:
    return len(re.findall(rf"({LEAF_TOKEN})(?=[,);])", _strip_lengths(t)))

def add_sister(tree: str, target_leaf: str, new_leaf: str) -> str:
    pat = re.compile(rf"(?<![A-Za-z0-9_.-]){re.escape(target_leaf)}(?=[,):])")
    if not pat.search(tree):
        sys.exit(f"[ERROR] target leaf '{target_leaf}' not found in HAL tree.")
    return pat.sub(f"({target_leaf},{new_leaf})", tree, count=1)

def graft_under(tree: str, internal_name: str, new_leaf: str) -> str:
    m = re.search(rf"(?<![A-Za-z0-9_.-]){re.escape(internal_name)}\(", tree)
    if not m:
        sys.exit(f"[ERROR] internal node '{internal_name}(' not found in HAL tree.")
    start, end = m.start(), m.end()
    return tree[:start] + f"({internal_name},{new_leaf})" + tree[end-1:]


# -------- HAL helpers --------

def hal_tree(hal: str, halstats_bin: str) -> str:
    if not shutil.which(halstats_bin):
        sys.exit(f"[ERROR] halStats not found: {halstats_bin}")
    try:
        r = subprocess.run([halstats_bin, "--tree", hal], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        sys.exit(f"[ERROR] halStats --tree failed ({e.returncode})")
    t = r.stdout.strip()
    if not t:
        sys.exit("[ERROR] halStats returned empty tree.")
    return _ensure_semicolon(t)

def hal_genomes(hal: str, halstats_bin: str):
    if not shutil.which(halstats_bin):
        sys.exit(f"[ERROR] halStats not found: {halstats_bin}")
    try:
        r = subprocess.run([halstats_bin, "--genomes", hal], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        sys.exit(f"[ERROR] halStats --genomes failed ({e.returncode})")
    return set(x.strip() for x in r.stdout.split() if x.strip())


# -------- IO helpers --------

def safe_symlink(src: str, dst: str):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)

def write_input_txt(path: str, tree: str, name: str, fasta_path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(_ensure_semicolon(tree) + "\n")
        fh.write(f"{name}\t{fasta_path}\n")


# -------- Main --------

def main():
    ap = argparse.ArgumentParser(description="Exact cactus-update-prepare wrapper for add node.")
    ap.add_argument("--in-hal", required=True, help="Existing HAL to update (e.g., evolverMammals.hal).")
    ap.add_argument("--out-prefix", required=True, help="Output prefix directory base (e.g., newMammals).")
    ap.add_argument("--new-name", required=True, help="Name of the new (or existing, for replace) genome.")
    ap.add_argument("--new-fasta", required=True, help="Path to new FASTA/2bit.")

    ap.add_argument("--halstats-bin", default="halStats", help="Path to halStats [default: halStats].")
    ap.add_argument("--prepare-bin", default="cactus-update-prepare", help="Path to cactus-update-prepare [default: cactus-update-prepare].")
    ap.add_argument("--max-cores", type=int, default=12, help="--maxCores (passed via --cactus-prepare-options).")
    ap.add_argument("--default-memory", default="200G", help="--defaultMemory (via --cactus-prepare-options).")
    ap.add_argument("--default-disk", default="400G", help="--defaultDisk (via --cactus-prepare-options).")

    sub = ap.add_subparsers(dest="mode", required=True)

    sp = sub.add_parser("sister", help="Add NEW as sister to an existing leaf.")
    sp.add_argument("--target-leaf", required=True)

    an = sub.add_parser("ancestor", help="Add NEW under a named internal (ancestor) node.")
    an.add_argument("--parent-node", required=True)

    rp = sub.add_parser("replace", help="Replace assembly for an existing leaf (no tree edit).")
    # no extra args

    args = ap.parse_args()

    outdir = args.out_prefix.rstrip("/")

    steps_dir = os.path.join(outdir, "steps")
    jobstore_dir = os.path.join(outdir, "jobstore")
    inputs_dir = os.path.join(outdir, "inputs")
    input_txt = os.path.join(outdir, "input.txt")
    run_sh = outdir + ".sh" if not args.out_prefix.endswith(".sh") else args.out_prefix

    for d in (outdir, steps_dir, jobstore_dir, inputs_dir):
        os.makedirs(d, exist_ok=True)

    # stage FASTA in inputs/ with stable name
    # keep original extension if recognizable
    base_ext = os.path.splitext(args.new_fasta)[1]
    if args.new_fasta.endswith(".fa.gz") or args.new_fasta.endswith(".fasta.gz"):
        staged = os.path.join(inputs_dir, f"{args.new_name}.fa.gz")
    elif base_ext:
        staged = os.path.join(inputs_dir, f"{args.new_name}{base_ext}")
    else:
        staged = os.path.join(inputs_dir, f"{args.new_name}.fa")
    safe_symlink(args.new_fasta, staged)

    # read current HAL tree and genomes
    tree = hal_tree(args.in_hal, args.halstats_bin)
    genomes = hal_genomes(args.in_hal, args.halstats_bin)

    # edit tree based on mode
    if args.mode == "sister":
        if args.target_leaf not in genomes:
            sys.exit(f"[ERROR] target leaf '{args.target_leaf}' not present in HAL.")
        if args.new_name in genomes:
            sys.exit(f"[ERROR] new-name '{args.new_name}' already exists in HAL; choose another or use 'replace'.")
        new_tree = add_sister(tree, args.target_leaf, args.new_name)

    elif args.mode == "ancestor":
        if args.new_name in genomes:
            sys.exit(f"[ERROR] new-name '{args.new_name}' already exists in HAL; choose another or use 'replace'.")
        new_tree = graft_under(tree, args.parent_node, args.new_name)

    else:  # replace
        if args.new_name not in genomes:
            sys.exit(f"[ERROR] replace mode expects existing leaf name as --new-name; '{args.new_name}' not found in HAL.")
        new_tree = tree  # no edit

    if _count_leaves(new_tree) < 2:
        sys.exit("[ERROR] Edited tree has <2 leaves; aborting.")

    # write input.txt with tree first line + new/replaced mapping
    write_input_txt(input_txt, new_tree, args.new_name, staged)

    # build exact prepare command
    if not shutil.which(args.prepare_bin):
        sys.exit(f"[ERROR] cactus-update-prepare not found: {args.prepare_bin}")

# decide which name goes on --genome (placement anchor)
	if args.mode == "sister":
    	genome_flag_value = args.target_leaf
	elif args.mode == "ancestor":
    	genome_flag_value = args.parent_node
	else:  # replace
    	genome_flag_value = args.new_name  # replacing an existing leaf

	cpo = f" --maxCores {args.max_cores} --defaultMemory {args.default_memory} --defaultDisk {args.default_disk} "
	cmd = [
    	args.prepare_bin,
    	"add", "node",
    	args.in_hal,
    	input_txt,
    	"--genome", genome_flag_value,       # <-- fixed here
    	"--outDir", steps_dir,
    	"--jobStore", jobstore_dir,
    	"--cactus-prepare-options", cpo,
	]

    print("[INFO] Running:", " ".join(cmd))
    try:
        r = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        sys.exit(f"[ERROR] cactus-update-prepare failed ({e.returncode})")

    header = [
        f"# date : {datetime.now()}",
        f"## cactus-update-prepare : {shutil.which(args.prepare_bin) or args.prepare_bin}",
        f"## HAL : {args.in_hal}",
        f"## outdir : {outdir}",
        f"## steps : {steps_dir}",
        f"## jobstore : {jobstore_dir}",
        "## exact command:",
        "## " + " ".join(cmd),
        "",
    ]
    with open(run_sh, "w") as fh:
        fh.write("\n".join(header))
        fh.write(r.stdout or "")
    os.chmod(run_sh, 0o755)

    print(f"[OK] Wrote {input_txt}")
    print(f"[OK] Staged FASTA  → {staged}")
    print(f"[OK] Pipeline      → {run_sh}")
    print(f"[HINT] Run: bash {run_sh}")

if __name__ == "__main__":
    main()
