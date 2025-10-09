#!/usr/bin/env python3
import os, re, sys, argparse, shutil, statistics, datetime
from pathlib import Path

LOG_PATTERNS = {
    "cv": [
        re.compile(r"CV\s*error\s*[:=]\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
        re.compile(r"cross[-\s]?validation\s*error\s*[:=]\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
    ],
    "llh": [
        re.compile(r"(?:like|log[-\s]?lik|loglik)\s*[:=]\s*(-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
        re.compile(r"best\s+like(?:lihood)?\s*[:=]\s*(-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
        re.compile(r"Final\s+log[-\s]?lik(?:elihood)?\s*[:=]\s*(-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"),
    ],
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Estimate best K from NGSadmix runs; copy .fopt and .filter for the best K value(s) with flexible naming.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--indir", "-i", default=".", help="Input directory to scan for logs, .fopt and .filter files")
    p.add_argument("--outdir", "-o", default="bestK", help="Directory to write outputs into")
    gk = p.add_argument_group("K selection")
    gk.add_argument("--k-list", "-K", type=str, default="", help="Comma-separated Ks to consider (e.g., 2,3,4,5). If empty, infer from files named like *K<k>*")
    gk.add_argument("--k-min", type=int, default=None, help="Minimum K to consider (inclusive) if inferring from files")
    gk.add_argument("--k-max", type=int, default=None, help="Maximum K to consider (inclusive) if inferring from files")
    p.add_argument("--log-glob", default="*.log", help="Glob to find log files used to score K (contains CV or LLH)")
    p.add_argument("--stat", choices=["cv","llh"], default="cv", help="Statistic to optimize: cv=minimize, llh=maximize")
    p.add_argument("--aggregate", choices=["mean","median","min","max","last"], default="mean", help="Aggregation across replicates")
    p.add_argument("--n-best", type=int, default=1, help="How many top K values to keep")
    p.add_argument("--descending", action="store_true", help="Force descending sort (useful if a custom stat where higher is better)")
    gname = p.add_argument_group("Naming / renaming")
    gname.add_argument("--bestk-prefix", default="", help="String to prepend to all emitted best-K filenames (e.g., 'chimps_')")
    gname.add_argument("--label-template", default="{prefix}bestK_K{K}", help="Filename label template. Available fields: {prefix}, {K}")
    gname.add_argument("--dir-name", default=None, help="If set, override the outdir's leaf name used in labels; does not change file system path")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing files in outdir")
    p.add_argument("--dry-run", action="store_true", help="Only report what would be done")
    return p.parse_args()

def discover_Ks(indir, k_list, k_min, k_max):
    if k_list:
        return sorted({int(x.strip()) for x in k_list.split(",") if x.strip()})
    # infer from filenames like *K5*.*
    ks = set()
    for path in Path(indir).rglob("*"):
        if path.is_file():
            m = re.search(r"[^\w]K(\d+)[^\w]", f"_{path.name}_")
            if m:
                ks.add(int(m.group(1)))
    if k_min is not None:
        ks = {k for k in ks if k >= k_min}
    if k_max is not None:
        ks = {k for k in ks if k <= k_max}
    return sorted(ks)

def extract_numbers_from_log(log_path, stat_key):
    pats = LOG_PATTERNS[stat_key]
    vals = []
    try:
        with open(log_path, "r", errors="ignore") as fh:
            for line in fh:
                for pat in pats:
                    m = pat.search(line)
                    if m:
                        try:
                            vals.append(float(m.group(1)))
                        except ValueError:
                            pass
    except Exception:
        pass
    return vals

def score_K(indir, Ks, log_glob, stat, aggregate):
    # collect per-K values
    scores = {}
    for k in Ks:
        vals = []
        # prefer files that explicitly include K<k> in name
        pri = list(Path(indir).rglob(f"*K{k}*"))
        logs = [p for p in pri if p.is_file() and p.match(log_glob)]
        # if none, fall back to any logs in indir
        if not logs:
            logs = [p for p in Path(indir).rglob(log_glob) if p.is_file()]
        for lg in logs:
            v = extract_numbers_from_log(lg, stat)
            vals.extend(v)
        if not vals:
            # no values found; give None to skip later
            scores[k] = None
            continue
        if aggregate == "mean":
            agg = statistics.fmean(vals)
        elif aggregate == "median":
            agg = statistics.median(vals)
        elif aggregate == "min":
            agg = min(vals)
        elif aggregate == "max":
            agg = max(vals)
        elif aggregate == "last":
            agg = vals[-1]
        scores[k] = agg
    return scores

def choose_best(scores, stat, n_best, force_desc=False):
    # Filter out None
    items = [(k,v) for k,v in scores.items() if v is not None]
    if not items:
        return []
    # Lower is better for CV; higher is better for LLH
    reverse = (stat == "llh")
    if force_desc:
        reverse = True
    items.sort(key=lambda kv: kv[1], reverse=reverse)
    return [k for k,_ in items[:n_best]]

def most_recent(paths):
    if not paths:
        return None
    paths = sorted(paths, key=lambda p: p.stat().st_mtime, reverse=True)
    return paths[0]

def find_artifacts_for_K(indir, K):
    # try common layouts
    patterns = [
        f"*K{K}*.fopt", f"*K{K}*.filter",
        f"K{K}/*.fopt", f"K{K}/*.filter",
        f"*K{K}*/*.fopt", f"*K{K}*/*.filter",
    ]
    fopt_candidates, filt_candidates = [], []
    for pat in patterns:
        for p in Path(indir).rglob(pat):
            if p.suffix == ".fopt":
                fopt_candidates.append(p)
            elif p.suffix == ".filter" or p.name.endswith(".filter"):
                filt_candidates.append(p)
    return most_recent(fopt_candidates), most_recent(filt_candidates)

def main():
    args = parse_args()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    Ks = discover_Ks(indir, args.k_list, args.k_min, args.k_max)
    if not Ks:
        print("ERROR: Could not infer any K values. Provide --k-list or ensure filenames contain 'K<k>'.", file=sys.stderr)
        sys.exit(2)

    scores = score_K(indir, Ks, args.log_glob, args.stat, args.aggregate)
    bestKs = choose_best(scores, args.stat, args.n_best, force_desc=args.descending)
    if not bestKs:
        print("ERROR: Failed to compute scores for any K. Do your logs include CV error or log-likelihood?", file=sys.stderr)
        sys.exit(3)

    # write a summary
    summary_path = outdir / "bestK_summary.tsv"
    with open(summary_path, "w") as fh:
        print("#time", "K", "score", "stat", "aggregate", sep="\t", file=fh)
        now = datetime.datetime.now().isoformat(timespec="seconds")
        for k in sorted(scores):
            print(now, k, scores[k] if scores[k] is not None else "NA", args.stat, args.aggregate, sep="\t", file=fh)

    # Copy artifacts
    emitted = []
    for k in bestKs:
        label = args.label_template.format(prefix=args.bestk_prefix, K=k)
        fopt_src, filt_src = find_artifacts_for_K(indir, k)
        if not fopt_src and not filt_src:
            print(f"WARNING: No .fopt or .filter found for K={k}. Skipping copy.", file=sys.stderr)
            continue
        if fopt_src:
            dst = outdir / f"{label}.fopt"
            if dst.exists() and not args.overwrite:
                print(f"NOTE: {dst} exists; use --overwrite to replace.", file=sys.stderr)
            elif not args.dry_run:
                shutil.copy2(fopt_src, dst)
            emitted.append(("fopt", k, str(fopt_src), str(dst)))
        if filt_src:
            dst = outdir / f"{label}.filter"
            if dst.exists() and not args.overwrite:
                print(f"NOTE: {dst} exists; use --overwrite to replace.", file=sys.stderr)
            elif not args.dry_run:
                shutil.copy2(filt_src, dst)
            emitted.append(("filter", k, str(filt_src), str(dst)))

    # Emit manifest
    manifest = outdir / "bestK_manifest.tsv"
    with open(manifest, "w") as fh:
        print("type", "K", "source", "dest", sep="\t", file=fh)
        for row in emitted:
            print(*row, sep="\t", file=fh)

    # Report
    print(f"[bestK] Candidates scored: {', '.join(map(str, Ks))}")
    print(f"[bestK] Best K(s): {', '.join(map(str, bestKs))}  (stat={args.stat}, agg={args.aggregate})")
    print(f"[bestK] Summary: {summary_path}")
    print(f"[bestK] Manifest: {manifest}")
    if emitted:
        print(f"[bestK] Copied {len(emitted)} artifact(s) into {outdir}")
    else:
        print(f"[bestK] No artifacts copied; check warnings above.")

if __name__ == "__main__":
    main()
