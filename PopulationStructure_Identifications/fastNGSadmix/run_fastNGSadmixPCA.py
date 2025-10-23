#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch runner for FastNGSAdmixPCA.R (on PATH).

Modes:
  - global (default): NO -populations flag; outputs -> {outdir}/global/
  - taxa            : pass -populations using taxa (with alias/matching); outputs -> {outdir}/taxa/<taxa>/
  - populations     : pass -populations using population (with alias/matching); outputs -> {outdir}/populations/<population>/

Flat layout in all modes (no per-ID subdirs). Files are prefixed by sample ID.

Sample list (tab-delimited):  id    taxa    population

Alias TSV format (tab-delimited):
  <key>\t<label1,label2,label3>
Example:
  Eastern\tEasternDRCSouthRwanda,EasternTanzania,EasternDRCNorthUganda
"""

import argparse
import csv
import os
import re
import sys
import shlex
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# ---------- global tee logger ----------
class TeeLogger:
    def __init__(self, path: Optional[str] = None):
        self.lock = threading.Lock()
        self.fh = open(path, "w") if path else None  # overwrite by default
    def close(self):
        if self.fh:
            self.fh.close()
            self.fh = None
    def log_text(self, text: str, end: str = "\n"):
        with self.lock:
            sys.stdout.write(text + end)
            sys.stdout.flush()
            if self.fh:
                self.fh.write(text + end)
                self.fh.flush()
    def log_bytes_line(self, b: bytes):
        # decode R output safely and mirror it
        s = b.decode("utf-8", errors="replace")
        with self.lock:
            sys.stdout.write(s)
            sys.stdout.flush()
            if self.fh:
                self.fh.write(s)
                self.fh.flush()

# ---------- utils ----------
def ensure_dir(p: str):
    if p and not os.path.isdir(p):
        os.makedirs(p, exist_ok=True)

def file_exists_nonempty(p: str) -> bool:
    return os.path.isfile(p) and os.path.getsize(p) > 0

def read_sample_list(path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str,str]] = []
    with open(path, 'r', newline='') as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        need = ['id','taxa','population']
        miss = [c for c in need if c not in rdr.fieldnames]
        if miss:
            raise ValueError(f"Sample list missing columns: {miss}. Need: id\\ttaxa\\tpopulation")
        for row in rdr:
            rows.append({'id':(row['id'] or '').strip(),
                         'taxa':(row['taxa'] or '').strip(),
                         'population':(row['population'] or '').strip()})
    seen = set()
    for r in rows:
        if not r['id']:
            raise ValueError("Empty sample ID encountered.")
        if r['id'] in seen:
            raise ValueError(f"Duplicate sample ID: {r['id']}")
        seen.add(r['id'])
    return rows

def filter_samples(samples: List[Dict[str,str]],
                   taxa_filter: Optional[str],
                   pop_filter: Optional[str]) -> List[Dict[str,str]]:
    out = samples
    if taxa_filter:
        keep = {t.strip() for t in taxa_filter.replace(';', ',').split(',') if t.strip()}
        out = [s for s in out if s['taxa'] in keep]
    if pop_filter:
        keep = {p.strip() for p in pop_filter.replace(';', ',').split(',') if p.strip()}
        out = [s for s in out if s['population'] in keep]
    return out

def out_dir_and_prefix(mode: str, outdir: str, sample: Dict[str,str]) -> Tuple[str, str]:
    """Flat layout in all modes; returns (directory, out_prefix)."""
    if mode == "global":
        d = os.path.join(outdir, "global")
    elif mode == "taxa":
        d = os.path.join(outdir, "taxa", sample['taxa'])
    else:  # populations
        d = os.path.join(outdir, "populations", sample['population'])
    return d, os.path.join(d, sample['id'])  # prefix == {dir}/{id}

def verify_inputs(likes: str, qopt: str, ref_prefix: str):
    if not os.path.isfile(likes): raise FileNotFoundError(f"likes not found: {likes}")
    if not os.path.isfile(qopt):  raise FileNotFoundError(f"qopt not found: {qopt}")
    for ext in ('.bed','.bim','.fam'):
        if not os.path.isfile(ref_prefix+ext):
            raise FileNotFoundError(f"Missing reference PLINK file: {ref_prefix+ext}")

# ---------- alias + qopt helpers ----------
def load_alias_file(path: Optional[str]) -> Dict[str, List[str]]:
    """Load alias TSV: key<TAB>csv_labels → dict[key] = [labels…]."""
    if not path:
        return {}
    m: Dict[str, List[str]] = {}
    with open(path, 'r') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n\r').split('\t', 1)
            if len(parts) != 2:
                continue
            key, csv_labels = parts
            labs = [x.strip() for x in csv_labels.split(',') if x.strip()]
            if labs:
                m[key] = labs
    return m

def read_qopt_header(qopt_path: str) -> List[str]:
    """Return qopt column names (first line, tab-delimited), excluding common non-pop headers."""
    with open(qopt_path, 'r') as fh:
        header = fh.readline().rstrip("\n\r")
    cols = [c.strip() for c in header.split('\t') if c.strip()]
    bad = {"ind", "indv", "id", "sample", "samples", "name"}
    return [c for c in cols if c.lower() not in bad]

def resolve_populations_list(qopt_path: str,
                             group_value: str,
                             match_mode: str = "auto",
                             alias_map: Optional[Dict[str, List[str]]] = None) -> List[str]:
    """Resolve a LIST of qopt columns to pass to -populations for taxa/populations modes."""
    if alias_map and group_value in alias_map:
        return alias_map[group_value]

    names = read_qopt_header(qopt_path)

    def exact():    return [c for c in names if c == group_value]
    def prefix():   return [c for c in names if c.startswith(group_value)]
    def contains(): return [c for c in names if group_value in c]
    def regex():
        try:
            rx = re.compile(group_value)
        except re.error as e:
            raise ValueError(f"Invalid regex for '{group_value}': {e}")
        return [c for c in names if rx.search(c)]

    mm = (match_mode or "auto").lower()
    if   mm == "exact":    matched = exact()
    elif mm == "prefix":   matched = prefix()
    elif mm == "contains": matched = contains()
    elif mm == "regex":    matched = regex()
    elif mm == "auto":     matched = exact() or prefix()
    else:
        raise ValueError(f"Unknown --pop-match: {match_mode}")

    if not matched:
        avail = ", ".join(names)
        raise RuntimeError(f"No qopt columns match '{group_value}' (mode={match_mode}). Available: {avail}")
    return matched

# ---------- R invocation ----------
def build_r_cmd(r_bin: str, r_script: str, ref_prefix: str,
                likes: str, qopt: str, out_prefix: str,
                r_flags: Dict[str,str],
                populations_list: Optional[List[str]] = None) -> Tuple[List[str], str]:
    """
    Build the FastNGSAdmixPCA.R command.
    - If populations_list is provided (taxa/populations modes), include -populations (comma-joined).
    - If None (global mode), do NOT include -populations.
    Returns (cmd_list, printable_pops_string).
    """
    cmd = [r_bin, r_script,
           "-likes", likes,
           "-ref", ref_prefix,
           "-qopt", qopt,
           "-out", out_prefix]

    pops_str = ""
    if populations_list:
        pops_str = ",".join(populations_list)
        cmd += ["-populations", pops_str]

    forward = {
        "PCs":"-PCs",
        "multiCores":"-multiCores",
        "saveCovar":"-saveCovar",
        "ngsTools":"-ngsTools",
        "onlyPrior":"-onlyPrior",
        "overlapRef":"-overlapRef",
        "doPlots":"-doPlots",
        "withChr":"-withChr"
    }
    for k,f in forward.items():
        v = r_flags.get(k)
        if v is None or v == "": continue
        cmd += [f, str(v)]
    return cmd, pops_str

def run_one(sample: Dict[str,str], mode: str,
            likes_dir: str, qopt_dir: str, likes_ext: str, qopt_ext: str,
            outdir: str, r_bin: str, r_script: str, ref_prefix: str,
            r_flags: Dict[str,str], resume: bool, dryrun: bool,
            pop_match: str, alias_map: Dict[str, List[str]],
            logger: TeeLogger) -> Dict[str,str]:
    sid = sample['id']

    # Mode-specific dir membership checks (layout only)
    if mode == "taxa" and not sample['taxa']:
        logger.log_text(f"[SKIP] {sid} (taxa) -> missing taxa")
        return {"id":sid,"mode":mode,"status":"skipped_missing_group","rc":0}
    if mode == "populations" and not sample['population']:
        logger.log_text(f"[SKIP] {sid} (populations) -> missing population")
        return {"id":sid,"mode":mode,"status":"skipped_missing_group","rc":0}

    # IO
    likes = os.path.join(likes_dir, f"{sid}{likes_ext}")
    qopt  = os.path.join(qopt_dir,  f"{sid}{qopt_ext}")
    out_dir, out_prefix = out_dir_and_prefix(mode, outdir, sample)
    ensure_dir(out_dir)

    if resume and file_exists_nonempty(f"{out_prefix}_covar.txt"):
        logger.log_text(f"[SKIP] {sid} ({mode}) -> existing {out_prefix}_covar.txt")
        return {"id":sid,"mode":mode,"status":"skipped_existing","rc":0}

    try:
        verify_inputs(likes, qopt, ref_prefix)
    except Exception as e:
        logger.log_text(f"[ERROR] {sid} ({mode}) inputs: {e}")
        return {"id":sid,"mode":mode,"status":"failed","rc":1}

    # Populations list: only for taxa/populations modes; use alias if provided
    populations_list: Optional[List[str]] = None
    if mode == "taxa":
        try:
            populations_list = resolve_populations_list(qopt, sample['taxa'], pop_match, alias_map)
        except Exception as e:
            logger.log_text(f"[ERROR] {sid} (taxa) resolve populations: {e}")
            return {"id":sid,"mode":mode,"status":"failed","rc":1}
    elif mode == "populations":
        try:
            populations_list = resolve_populations_list(qopt, sample['population'], pop_match, alias_map)
        except Exception as e:
            logger.log_text(f"[ERROR] {sid} (populations) resolve populations: {e}")
            return {"id":sid,"mode":mode,"status":"failed","rc":1}

    cmd_list, pops_str = build_r_cmd(r_bin, r_script, ref_prefix, likes, qopt, out_prefix, r_flags, populations_list)
    cmd_str  = " ".join(shlex.quote(c) for c in cmd_list)

    if dryrun:
        if pops_str:
            logger.log_text(f"[DRYRUN][{mode}] {sid} | pops={pops_str}")
        else:
            logger.log_text(f"[DRYRUN][{mode}] {sid}")
        logger.log_text(f"  CMD: {cmd_str}")
        return {"id":sid,"mode":mode,"status":"dryrun","rc":0}

    if pops_str:
        logger.log_text(f"[START][{mode}] {sid} | pops={pops_str}")
    else:
        logger.log_text(f"[START][{mode}] {sid}")

    log_path = f"{out_prefix}.FastNGSAdmixPCA.{mode}.log"
    with open(log_path, "w") as logf:
        logf.write(f"# Launched: {datetime.now().isoformat()}\n# CMD: {cmd_str}\n\n"); logf.flush()
        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(proc.stdout.readline, b""):
            if not line: break
            # tee to terminal + global log, and also into per-sample log
            logger.log_bytes_line(line)
            logf.write(line.decode("utf-8", errors="replace"))
            logf.flush()
        rc = proc.wait()

    status = "ok" if rc==0 else "failed"
    logger.log_text(f"[END][{mode}] {sid} -> {status} (rc={rc})")
    return {"id":sid,"mode":mode,"status":status,"rc":rc}

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Batch launcher for FastNGSAdmixPCA.R with flat outputs; default 'global'. Uses -populations only in taxa/populations modes, with optional alias mapping."
    )
    # IO roots
    ap.add_argument("-l","--likes-dir", required=True, help="Dir with {id}{likes-ext} Beagle files")
    ap.add_argument("-q","--qopt-dir",  required=True, help="Dir with {id}{qopt-ext} qopt files")
    ap.add_argument("-o","--outdir",    required=True, help="Root output directory")
    ap.add_argument("-x","--likes-ext", default=".beagle.gz", help="Beagle suffix (default .beagle.gz)")
    ap.add_argument("-X","--qopt-ext",  default=".qopt", help="qopt suffix (default .qopt)")

    # Inputs
    ap.add_argument("-s","--sample-list", required=True, help="TSV: id\\ttaxa\\tpopulation")

    # Mode (no 'all'); default global
    ap.add_argument("-m","--mode", choices=["global","taxa","populations"], default="global",
                    help="Directory layout mode. -populations is only added in taxa/populations modes.")

    # Optional filters to subset which samples run
    ap.add_argument("-t","--taxa-filter", help="Only run samples with taxa in this CSV/semicolon list")
    ap.add_argument("-p","--population-filter", help="Only run samples with population in this CSV/semicolon list")

    # Matching behavior & alias for taxa/populations modes
    ap.add_argument("-M","--pop-match", choices=["auto","exact","prefix","contains","regex"], default="auto",
                    help="Match rule for mapping group label to qopt columns (default auto: exact then prefix)")
    ap.add_argument("-A","--alias", help="Alias TSV mapping group→CSV labels (e.g., Eastern → EasternDRC...,EasternTanzania,...)")

    # R & ref
    ap.add_argument("-b","--r-binary", default="Rscript", help="Rscript binary (default Rscript)")
    ap.add_argument("-S","--r-script", default="FastNGSAdmixPCA.R", help="R script on PATH (default FastNGSAdmixPCA.R)")
    ap.add_argument("-R","--ref-prefix", required=True, help="PLINK binary prefix for reference panel")

    # Forward R options (passed through unchanged)
    ap.add_argument("-P","--PCs", default="1,2,3", help="PCs to plot (default 1,2,3)")
    ap.add_argument("-c","--multiCores", default="1", help="Cores within R")
    ap.add_argument("--saveCovar", default="0", choices=["0","1"])
    ap.add_argument("--ngsTools",  default="0", choices=["0","1"])
    ap.add_argument("--onlyPrior", default="0", choices=["0","1"])
    ap.add_argument("--overlapRef", default="1", choices=["0","1"])
    ap.add_argument("--doPlots",   default="1", choices=["0","1"])
    ap.add_argument("--withChr",   default="0", choices=["0","1"])

    # Execution controls
    ap.add_argument("-j","--jobs", type=int, default=2, help="Parallel workers (default 2)")
    ap.add_argument("-u","--resume", action="store_true", help="Skip if *_covar.txt exists")
    ap.add_argument("-n","--dry-run", action="store_true", help="Print commands only")

    # Global log file
    ap.add_argument("-L","--log", help="Path to a global log file capturing ALL stdout (in addition to per-sample logs)")

    args = ap.parse_args()

    # Set up global logger
    logger = TeeLogger(args.log)

    try:
        # Load & filter samples
        samples = read_sample_list(args.sample_list)
        samples = filter_samples(samples, args.taxa_filter, args.population_filter)
        if not samples:
            logger.log_text("No samples to run after filtering.")
            sys.exit(0)

        # Drop unlabeled ones for layout-only modes
        if args.mode == "taxa":
            samples = [s for s in samples if s['taxa']]
        elif args.mode == "populations":
            samples = [s for s in samples if s['population']]
        if not samples:
            logger.log_text("No samples left after mode-based restriction.")
            sys.exit(0)

        # Forward flags to R + alias map
        alias_map = load_alias_file(args.alias) if args.alias else {}
        r_flags = dict(
            PCs=args.PCs, multiCores=args.multiCores,
            saveCovar=args.saveCovar, ngsTools=args.ngsTools, onlyPrior=args.onlyPrior,
            overlapRef=args.overlapRef, doPlots=args.doPlots, withChr=args.withChr
        )

        total_runs = len(samples)
        logger.log_text(f"[INFO] Scheduling {total_runs} run(s) (mode={args.mode}, samples={len(samples)}).")

        ensure_dir(args.outdir)
        results = []; completed = 0

        def submit(ex):
            futs = []
            for s in samples:
                futs.append(ex.submit(
                    run_one, s, args.mode,
                    args.likes_dir, args.qopt_dir, args.likes_ext, args.qopt_ext,
                    args.outdir, args.r_binary, args.r_script, args.ref_prefix,
                    r_flags, args.resume, args.dry_run, args.pop_match, alias_map, logger
                ))
            return futs

        if args.jobs <= 1:
            class _Dummy:
                def submit(self, fn, *a, **k):
                    class F:
                        def result(self_inner): return fn(*a, **k)
                    return F()
            ex = _Dummy(); futs = submit(ex)
            for f in futs:
                res = f.result(); results.append(res); completed += 1
                logger.log_text(f"[PROGRESS] completed {completed}/{total_runs}")
        else:
            with ThreadPoolExecutor(max_workers=args.jobs) as ex:
                for f in as_completed(submit(ex)):
                    res = f.result(); results.append(res); completed += 1
                    logger.log_text(f"[PROGRESS] completed {completed}/{total_runs}")

        failed  = sum(1 for r in results if r['status']=="failed")
        skipped = sum(1 for r in results if r['status'] in ("skipped_existing","skipped_missing_group"))
        dry     = sum(1 for r in results if r['status']=="dryrun")
        ran     = total_runs - skipped - dry
        logger.log_text("\n[summary] total_runs={}, ran_now={}, failed_now={}, skipped={}, dryrun={}".format(
            total_runs, ran, failed, skipped, dry
        ))

    finally:
        logger.close()

if __name__ == "__main__":
    main()
