#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download & merge SRA FASTQs per sample, choosing defline by headerType.
(Updated to use **fasterq-dump**)

INPUT TABLE (3 columns, TSV/CSV/whitespace OK):
  accession   sampleName   headerType

where headerType should be: 'instrument', 'numeric', or 'unknown'.
- instrument     -> --seq-defline '@$sn[_$rn]/$ri'    (preserves instrument-style)
- numeric/unknown-> --seq-defline '@$ac:$si/$ri'      (forces @ACCESSION:SPOTID/ri)

What it does
------------
1) Downloads each accession into a tmp dir via fasterq-dump:
   -e {fqd_threads} --split-files --skip-technical --seq-defline ... --qual-defline '+'
   (Note: fasterq-dump writes UNCOMPRESSED .fastq files; this script compresses them to .fastq.gz)
2) For each sample:
   - If a single accession → rename to: raw_fastq/<sampleName>_1.fastq.gz and ..._2.fastq.gz
   - If multiple accessions → merge each read file in sorted accession order
     into those same outputs (concatenate gz streams).
3) Optional cleanup of tmp dir.

Usage
-----
python download_and_merge_sra_fqd.py \
  -t runs.3col.tsv \
  -o OUTDIR \
  --jobs 4 \
  --fqd-threads 4 \
  --sratoolkit /path/to/sratoolkit \
  --tmpdir /scratch/with/space \
  --clean-up

Notes
-----
- Output directory layout:
    OUTDIR/
      raw_fastq/
      tmp_sra/   (deleted with --clean-up)
- Output filenames end with _1.fastq.gz / _2.fastq.gz
- Ensure **sufficient free space** in --tmpdir (fasterq-dump writes large temporary files).
- If outputs exist and --force not set, the sample is skipped.
"""
import argparse, csv, os, re, shutil, subprocess, sys
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Tuple, Dict

SRA_RUN_RE = re.compile(r'^[SED]RR[0-9]+$', re.IGNORECASE)

def parse_args():
    ap = argparse.ArgumentParser(description="Download & merge SRA FASTQs by sample (3-col table) using fasterq-dump.")
    ap.add_argument("-t","--table", required=True, help="3-column table: accession sampleName headerType")
    ap.add_argument("-o","--outdir", required=True, help="Base output dir; will create raw_fastq/ and tmp_sra/")
    ap.add_argument("--tmpdir", default="", help="Temporary dir for fasterq-dump (-t). Defaults to OUTDIR/tmp_sra")
    ap.add_argument("--sratoolkit", default="", help="Path to SRA toolkit root **or** to fasterq-dump binary")
    ap.add_argument("--jobs", type=int, default=4, help="Parallel downloads (default: 4)")
    ap.add_argument("--fqd-threads","-e", type=int, default=4, help="Threads for fasterq-dump -e (default: 4)")
    ap.add_argument("--clean-up","-c", action="store_true", help="Delete tmp_sra after merge")
    ap.add_argument("--force", action="store_true", help="Overwrite existing final FASTQs")
    return ap.parse_args()

def resolve_fasterq_dump(hint: str) -> str:
    # If hint is a directory, assume bin/fasterq-dump inside
    if hint:
        p = Path(hint)
        if p.is_dir():
            cand = p / "bin" / "fasterq-dump"
            if cand.exists() and os.access(cand, os.X_OK):
                return str(cand)
        # If hint looks like an executable path
        if Path(hint).exists() and os.access(hint, os.X_OK):
            return str(hint)
    exe = shutil.which("fasterq-dump")
    if not exe:
        sys.exit("ERROR: fasterq-dump not found (use --sratoolkit or add to PATH).")
    return exe

def read_table(path: str):
    with open(path, newline="") as fh:
        text = fh.read()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return []
    # detect delimiter on the first non-header data row
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    rows = list(csv.reader(lines, delimiter=delim)) if delim else [re.split(r"\s+", ln) for ln in lines]
    triples = []
    for i,row in enumerate(rows):
        if len(row) < 3:
            continue
        acc, sample, htype = row[0].strip(), row[1].strip(), row[2].strip().lower()
        # skip header row
        if i == 0 and (row[0].lower() == "accession" or not SRA_RUN_RE.match(acc)):
            continue
        if htype not in {"instrument","numeric","unknown"}:
            htype = "unknown"
        triples.append((acc, sample, htype))
    return triples

def defline_for(htype: str) -> str:
    if htype == "instrument":
        # Keep instrument-style names if present
        return "@$sn[_$rn]/$ri"
    # numeric / unknown
    return "@$ac:$si/$ri"

def ensure_dirs(base: Path):
    final_dir = base / "raw_fastq"
    tmp_dir = base / "tmp_sra"
    final_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return final_dir, tmp_dir

def run_checked(cmd: list, log_prefix=""):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{log_prefix} failed (rc={p.returncode})\nCMD: {' '.join(cmd)}\nSTDERR:\n{p.stderr}")
    return p

def compress_fastq(path: Path):
    """Compress a .fastq file to .fastq.gz using pigz if available, else gzip."""
    if not path.exists():
        raise FileNotFoundError(str(path))
    pigz = shutil.which("pigz")
    if pigz:
        subprocess.run([pigz, "-f", str(path)], check=True)
    else:
        subprocess.run(["gzip", "-f", str(path)], check=True)
    gz = path.with_suffix(path.suffix + ".gz")
    if not gz.exists():
        raise FileNotFoundError(str(gz))
    return gz

def download_one(fasterq_dump: str, acc: str, htype: str, tmp_dir: Path, fqd_threads: int):
    defline = defline_for(htype)
    cmd = [
        fasterq_dump,
        "-e", str(fqd_threads),
        "--split-files",
        "--skip-technical",
        "-t", str(tmp_dir),
        "--seq-defline", defline,
        "--qual-defline", "+",
        "-O", str(tmp_dir),
        acc
    ]
    print(f"[download] {acc} ({htype}) :: {' '.join(cmd)}", file=sys.stderr)
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        last = p.stderr.strip().splitlines()[-1:] or [p.stderr.strip()]
        print(f"[error] {acc} rc={p.returncode} :: {last[0]}", file=sys.stderr)
        if "storage exhausted" in p.stderr.lower() or "no space left" in p.stderr.lower():
            print("[hint] Temp/output disk is likely full. Try: --tmpdir /path/with/space, reduce --fqd-threads, or point --outdir to a larger disk.", file=sys.stderr)
        return False, acc

    # fasterq-dump writes uncompressed files: {acc}_1.fastq, {acc}_2.fastq (if paired)
    r1_plain = tmp_dir / f"{acc}_1.fastq"
    r2_plain = tmp_dir / f"{acc}_2.fastq"
    ok_any = False
    if r1_plain.exists():
        compress_fastq(r1_plain); ok_any = True
    if r2_plain.exists():
        compress_fastq(r2_plain); ok_any = True
    if not ok_any:
        print(f"[warn] {acc}: no output FASTQs found", file=sys.stderr)
        return False, acc
    return True, acc

def merge_or_rename(sample: str, accs, tmp_dir: Path, final_dir: Path, force: bool):
    out1 = final_dir / f"{sample}_1.fastq.gz"
    out2 = final_dir / f"{sample}_2.fastq.gz"
    if (out1.exists() or out2.exists()) and not force:
        print(f"[skip] {sample}: outputs exist (use --force to overwrite)", file=sys.stderr)
        return

    # collect per-read files in sorted accession order
    accs_sorted = sorted(accs)
    r1_files, r2_files = [], []
    for acc in accs_sorted:
        f1 = tmp_dir / f"{acc}_1.fastq.gz"
        f2 = tmp_dir / f"{acc}_2.fastq.gz"
        if f1.exists(): r1_files.append(str(f1))
        if f2.exists(): r2_files.append(str(f2))

    if not r1_files and not r2_files:
        print(f"[warn] {sample}: no FASTQs found in tmp for its accessions", file=sys.stderr)
        return

    # single-accession rename; else concatenate gz streams
    if len(accs_sorted) == 1:
        if r1_files:
            shutil.move(r1_files[0], out1)
        if r2_files:
            shutil.move(r2_files[0], out2)
        print(f"[rename] {sample}: wrote {out1.name if r1_files else ''} {out2.name if r2_files else ''}", file=sys.stderr)
    else:
        if r1_files:
            with open(out1, "wb") as w:
                for fp in r1_files:
                    with open(fp, "rb") as r:
                        shutil.copyfileobj(r, w)
        if r2_files:
            with open(out2, "wb") as w:
                for fp in r2_files:
                    with open(fp, "rb") as r:
                        shutil.copyfileobj(r, w)
        print(f"[merge] {sample}: merged {len(r1_files)} R1 and {len(r2_files)} R2 parts", file=sys.stderr)

def main():
    args = parse_args()
    base = Path(args.outdir).resolve()
    final_dir, tmp_dir = ensure_dirs(base)
    if args.tmpdir:
        tmp_dir = Path(args.tmpdir).resolve()
        tmp_dir.mkdir(parents=True, exist_ok=True)
    fasterq_dump = resolve_fasterq_dump(args.sratoolkit)

    triples = read_table(args.table)
    if not triples:
        sys.exit("ERROR: no rows parsed from 3-column table.")

    # Map sample -> list of accessions; also record header type per accession
    sample2accs = defaultdict(list)
    acc2type = {}
    for acc, sample, htype in triples:
        sample2accs[sample].append(acc)
        acc2type[acc] = htype

    # Download all accessions in parallel
    accs = sorted(acc2type.keys())
    print(f"[plan] accessions: {len(accs)}  samples: {len(sample2accs)}", file=sys.stderr)
    ok_map = {}

    # Control parallelism with --jobs; each job runs fasterq-dump with -e args.fqd_threads
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = [ex.submit(download_one, fasterq_dump, acc, acc2type[acc], tmp_dir, args.fqd_threads) for acc in accs]
        for fut in as_completed(futs):
            ok, acc = fut.result()
            ok_map[acc] = ok

    # Merge per sample (only using successful downloads)
    for sample, acc_list in sample2accs.items():
        acc_ok = [a for a in acc_list if ok_map.get(a, False)]
        if not acc_ok:
            print(f"[warn] {sample}: no successful downloads", file=sys.stderr)
            continue
        merge_or_rename(sample, acc_ok, tmp_dir, final_dir, args.force)

    if args.clean_up:
        try:
            shutil.rmtree(tmp_dir)
            print(f"[clean] removed {tmp_dir}", file=sys.stderr)
        except Exception as e:
            print(f"[warn] could not remove tmp dir: {e}", file=sys.stderr)

    # Print final dir path for downstream tooling
    print(str(final_dir))

if __name__ == "__main__":
    main()
