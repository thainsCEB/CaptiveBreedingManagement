#!/usr/bin/env python3
"""
Run fastNGSadmix per sample, auto-resume, then merge Q matrices and (optionally) write a summary.

Key features
------------
- Auto-resume: samples with existing Q files in OUTDIR are skipped.
- Use --force / -F to re-run all samples.
- --id-list / -i : file with sample IDs in column 1 (header OK)
- --like-dir / -l: directory containing {sample}.beagle.gz
- --nname   / -n : -Nname file for fastNGSadmix
- --fname   / -f : -fname file for fastNGSadmix
- --outdir  / -o : output directory; creates OUTDIR/logs/ for per-sample stdout/stderr
- Parallel execution (--procs / -p)
- Merges per-sample Q to OUTDIR/combined_Q.tsv: [sample_id, pop1..popK]
- Optional summary (--summary / -y) to OUTDIR/summary.tsv:
    [sample_id, <taxon1>(%).. <taxonK>(%), Putative_taxon, Status]
  Putative_taxon uses --taxa-map / -T labels (2-col TSV: pop / label)
  Status = 'pure' if maxQ >= 0.9999 else 'hybrid'
  For hybrids: taxa with >=10% are listed alphabetically and joined with ' x ' (else top two).
"""

import argparse
import csv
import os
import re
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple

PURE_THRESHOLD = 0.9999
HYBRID_LIST_THRESHOLD = 0.005  # taxa >= 0.5% are included in hybrid label

# ------------------------- CLI ------------------------- #

def parse_args():
    ap = argparse.ArgumentParser(description="Run fastNGSadmix per ID and merge Q matrices.")
    # Inputs
    ap.add_argument("--id-list", "-i", required=True,
                    help="File with sample IDs in column 1 (header allowed).")
    ap.add_argument("--nname", "-n", required=True,
                    help="-Nname file for fastNGSadmix.")
    ap.add_argument("--fname", "-f", required=True,
                    help="-fname file for fastNGSadmix.")
    ap.add_argument("--like-dir", "-l", required=True,
                    help="Directory with {sample}.beagle.gz files.")
    # fastNGSadmix options
    ap.add_argument("--which-pops", "-w", default="all",
                    help="Value for -whichPops (default: all).")
    ap.add_argument("--out-suffix", "-O", default=".fastNGSadmix_chimpsv1",
                    help="Suffix for per-sample output basename (default: .fastNGSadmix_chimpsv1).")
    # Paths / execution
    ap.add_argument("--outdir", "-o", required=True,
                    help="Output directory (creates logs/ inside).")
    ap.add_argument("--procs", "-p", type=int, default=4,
                    help="Parallel jobs (default: 4).")
    ap.add_argument("--force", "-F", action="store_true",
                    help="Re-run all samples even if Q files already exist in OUTDIR.")
    ap.add_argument("--dry-run", "-d", action="store_true",
                    help="Print commands but do not execute.")
    # Outputs
    ap.add_argument("--merged-q", "-m", default=None,
                    help="Path for merged Q TSV (default: OUTDIR/combined_Q.tsv).")
    ap.add_argument("--summary", "-y", action="store_true",
                    help="Write OUTDIR/summary.tsv with %% per taxon + Putative_taxon + Status.")
    ap.add_argument("--taxa-map", "-T", default=None,
                    help="2-col TSV mapping pop(1..K or pop1..popK) -> taxon label (used in summary).")
    return ap.parse_args()

# ------------------------- IO helpers ------------------------- #

def read_ids(path: str) -> List[str]:
    ids = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            # Skip header-like first line
            if i == 1 and any(h in line.lower() for h in ("sample", "id")):
                continue
            ids.append(line.split("\t")[0])
    return ids

def first_existing_q(base: str) -> Optional[str]:
    for ext in (".Q", ".qopt", ".qposterior_means"):
        p = f"{base}{ext}"
        if os.path.isfile(p):
            return p
    return None

def q_exists_for_sample(sample_id: str, out_suffix: str, outdir: str) -> bool:
    base = os.path.join(outdir, f"{sample_id}{out_suffix}")
    return first_existing_q(base) is not None

# ------------------------- fastNGSadmix runner ------------------------- #

def run_fastngs_for_sample(sample_id: str, like_dir: str, nname: str, fname: str,
                           which_pops: str, out_suffix: str, outdir: str,
                           dry_run: bool = False, log_dir: Optional[str] = None) -> Tuple[int, str]:
    likes = os.path.join(like_dir, f"{sample_id}.beagle.gz")
    outbase = os.path.join(outdir, f"{sample_id}{out_suffix}")

    cmd = [
        "fastNGSadmix",
        "-likes", likes,
        "-Nname", nname,
        "-fname", fname,
        "-whichPops", str(which_pops),
        "-out", outbase,
    ]

    if dry_run:
        print("[DRY-RUN]", " ".join(cmd))
        return 0, sample_id

    if not os.path.isfile(likes):
        return 2, f"{sample_id}: missing likes file {likes}"

    os.makedirs(outdir, exist_ok=True)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    if log_dir:
        log_out = os.path.join(log_dir, f"{sample_id}.stdout.log")
        log_err = os.path.join(log_dir, f"{sample_id}.stderr.log")
        with open(log_out, "w") as out, open(log_err, "w") as err:
            ret = subprocess.run(cmd, stdout=out, stderr=err)
    else:
        ret = subprocess.run(cmd)

    if ret.returncode != 0:
        return ret.returncode, f"{sample_id}: fastNGSadmix failed (code {ret.returncode})"
    return 0, sample_id

# ------------------------- Parsing utilities ------------------------- #

_float_re = re.compile(r"""
    ^[+-]?(
        (\d+(\.\d*)?)|(\.\d+)
    )([eE][+-]?\d+)?$
""", re.VERBOSE)

def _parse_float_tokens(tokens: List[str]) -> List[float]:
    """Return floats from tokens; ignore non-numeric tokens like 'K1'."""
    vals = []
    for t in tokens:
        tt = t.strip()
        if _float_re.match(tt):
            try:
                vals.append(float(tt))
            except ValueError:
                pass
    return vals

def _read_q_values(qpath: str) -> Optional[List[float]]:
    """
    Robust reader for fastNGSadmix Q-like files:
    - Skips empty/comment/header lines (e.g., 'K1 K2 ...')
    - Returns first line that yields >=1 float
    """
    with open(qpath) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            tokens = re.split(r"\s+", ln)
            vals = _parse_float_tokens(tokens)
            if vals:
                return vals
    return None

# ------------------------- Merge & Summary ------------------------- #

def merge_q(ids: List[str], out_suffix: str, outdir: str, merged_out_path: Optional[str] = None):
    rows = []  # (sample_id, [q1..qK])
    K = None
    missing = []
    for sid in ids:
        base = os.path.join(outdir, f"{sid}{out_suffix}")
        qpath = first_existing_q(base)
        if qpath is None:
            missing.append(sid)
            continue
        vals = _read_q_values(qpath)
        if vals is None:
            missing.append(sid)
            continue
        if K is None:
            K = len(vals)
        elif K != len(vals):
            raise SystemExit(f"[ERROR] Mixed K: {sid} has {len(vals)}, earlier was {K}")
        rows.append((sid, vals))

    if K is None:
        raise SystemExit("[ERROR] No Q data found to merge.")

    merged_out = merged_out_path or os.path.join(outdir, "combined_Q.tsv")
    with open(merged_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id"] + [f"pop{i+1}" for i in range(K)])
        for sid, vals in rows:
            w.writerow([sid] + [f"{v:.6f}" for v in vals])

    print(f"[OK] Wrote {merged_out} with {len(rows)} samples; K={K}")
    if missing:
        sys.stderr.write(f"[WARN] Missing Q for {len(missing)} samples: "
                         f"{', '.join(missing[:15])}{' ...' if len(missing) > 15 else ''}\n")
    return rows, K, merged_out

def read_taxa_map(path: str, K_expected: Optional[int] = None) -> Optional[List[str]]:
    """
    Read 2-col TSV: [pop, label], where pop is 1..K or 'pop1'..'popK'.
    Returns list of length K (labels), filling missing with 'pop{i}'.
    """
    if not path:
        return None
    mapping = {}
    with open(path) as f:
        rdr = csv.reader(f, delimiter="\t")
        for row in rdr:
            if not row or not row[0].strip():
                continue
            key = row[0].strip()
            label = row[1].strip() if len(row) > 1 else ""
            m = re.match(r"^(?:pop)?(\d+)$", key, flags=re.IGNORECASE)
            if not m:
                raise SystemExit(f"[ERROR] taxa-map bad key '{key}' (use 1..K or pop1..popK)")
            idx = int(m.group(1))  # 1-based
            mapping[idx] = label
    if K_expected:
        return [mapping.get(i+1, f"pop{i+1}") for i in range(K_expected)]
    keys = sorted(mapping.keys())
    return [mapping[k] for k in keys]

def _hybrid_label(vals: List[float], labels: List[str]) -> str:
    """
    For hybrids, list taxa alphabetically (by label) whose proportion >= HYBRID_LIST_THRESHOLD,
    joined with ' x '. If none meet the threshold, use the top two taxa by proportion.
    """
    parts = [(labels[i], vals[i]) for i in range(len(vals))]
    keep = [name for name, v in parts if v >= HYBRID_LIST_THRESHOLD]
    if not keep:
        top2_idx = sorted(range(len(vals)), key=lambda i: vals[i], reverse=True)[:2]
        keep = [labels[i] for i in top2_idx]
    keep = sorted(set(keep), key=lambda s: s.lower())
    return " x ".join(keep)

def write_summary(rows, K, outdir, taxa_labels: Optional[List[str]] = None) -> str:
    """
    rows: [(sample_id, [q1..qK])]
    taxa_labels: list[str] length K; if None, use pop1..K
    """
    labels = taxa_labels if (taxa_labels and len(taxa_labels) == K) else [f"pop{i+1}" for i in range(K)]
    outp = os.path.join(outdir, "summary.tsv")
    header = ["sample_id"] + [f"{lab}(%)" for lab in labels] + ["Putative_taxon", "Status"]
    with open(outp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for sid, vals in rows:
            perc = [v * 100.0 for v in vals]
            max_i = max(range(K), key=lambda i: vals[i])
            max_q = vals[max_i]
            if max_q >= PURE_THRESHOLD:
                putative = labels[max_i]
                status = "pure"
            else:
                putative = _hybrid_label(vals, labels)
                status = "hybrid"
            w.writerow([sid] + [f"{p:.2f}" for p in perc] + [putative, status])
    print(f"[OK] Wrote {outp}")
    return outp

# ------------------------- Main ------------------------- #

def main():
    args = parse_args()

    # Basic validations
    if not isinstance(args.out_suffix, str) or not args.out_suffix:
        raise SystemExit("--out-suffix must be a non-empty string")

    ids = read_ids(args.id_list)
    if not ids:
        raise SystemExit("[ERROR] No IDs found.")

    os.makedirs(args.outdir, exist_ok=True)
    log_dir = os.path.join(args.outdir, "logs")
    if not args.dry_run:
        os.makedirs(log_dir, exist_ok=True)

    # Decide which samples to run (auto-resume)
    if args.force:
        to_run = ids[:]
        already = []
    else:
        already = [sid for sid in ids if q_exists_for_sample(sid, args.out_suffix, args.outdir)]
        to_run = [sid for sid in ids if sid not in already]

    print(f"[info] Total: {len(ids)} | already have Q: {len(already)} | to run now: {len(to_run)}")
    if already:
        print(f"[info] Resuming — skipping existing Q for: {', '.join(already[:20])}"
              f"{' ...' if len(already) > 20 else ''}")

    # Run fastNGSadmix in parallel for missing ones (or all if --force)
    n_ok = n_fail = 0
    if to_run and not args.dry_run:
        with ThreadPoolExecutor(max_workers=args.procs) as ex:
            futures = [ex.submit(
                run_fastngs_for_sample,
                sid, args.like_dir, args.nname, args.fname,
                args.which_pops, args.out_suffix, args.outdir,
                args.dry_run, log_dir
            ) for sid in to_run]
            for fut in as_completed(futures):
                code, msg = fut.result()
                if code == 0:
                    n_ok += 1
                else:
                    n_fail += 1
                    sys.stderr.write(f"[ERR] {msg}\n")
    print(f"[done] ran_now={n_ok}, failed_now={n_fail}, skipped_existing={len(already)}")

    if args.dry_run:
        return

    # Merge Qs for ALL ids, using whatever exists
    rows, K, merged_path = merge_q(ids, args.out_suffix, args.outdir, args.merged_q)

    # Summary
    if args.summary:
        taxa_labels = read_taxa_map(args.taxa_map, K_expected=K) if args.taxa_map else None
        write_summary(rows, K, args.outdir, taxa_labels)

if __name__ == "__main__":
    main()
