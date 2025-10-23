#!/usr/bin/env python3
"""
ld_prune_beagle_resumable_cli.py

Backwards-compatible CLI with added features:
- Resumable steps with per-step logs
- Defaults: ngsLD --max_kb_dist 1000, --min_maf 0.05
- prune_graph: uses --in/--out and applies weight filter:
    --weight-field column_7
    --weight-filter "column_3 > {min_dist} && column_7 > {min_r2}"
  (so -d/--min-dist and -r/--min-r2 still matter, but only in the filter)
- Threads: -j (chromosome-level parallel jobs), -t (threads per job; ngsLD and prune_graph via -n/OMP)
- Prefix-first layout (and inference): prefer {PREFIX}.{CHR}.suffix; write {PREFIX}.{CHR}.pos by default
- Position derivation: .pos -> from .mafs.gz -> from .beagle.gz (marker chr:pos or chr_pos)
- Sites outputs per chromosome (default): 
    {CHR}.unlinked_ids.txt   (chr_pos)
    {CHR}.unlinked_sites.tsv (chr<TAB>pos)
    {CHR}.sites.tsv          (chr<TAB>pos<TAB>major<TAB>minor) via MAFS; fallback to Beagle alleles
- Merging:
    Full merged GL:        {gl_basename}.beagle.gz
    Pruned merged GL:      {pruned_gl_basename}.beagle.gz (rows for union of unlinked sites)
    Pruned merged MAFS:    {pruned_gl_basename}.mafs.gz   (rows for union of unlinked sites)
    Merged unlinked lists: merged.unlinked_ids.txt / merged.unlinked_sites.tsv / merged.sites.tsv
  Include prior runs' pruned outputs via --external-pruned-list (full paths to .ld.pruned).

Verbose is the default (no -v flag).
"""

import argparse, os, sys, subprocess, gzip, shutil, re
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

def chrom_key(c: str):
    # Natural numeric order; map common non-numeric labels if encountered
    m = re.match(r'^(?:chr)?(\d+)$', c, re.IGNORECASE)
    if m:
        return (0, int(m.group(1)))
    # Map X/Y/MT to large numbers to follow numeric autosomes
    c_lower = c.lower().lstrip('chr')
    mapping = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
    return (1, mapping.get(c_lower, 1_000_000), c_lower)

PREFIX = ""  # set after parsing (explicit -x or inferred)

# ---------- Utilities ----------

def shlex_join(cmd: List[str]) -> str:
    import shlex
    return " ".join(shlex.quote(c) for c in cmd)

def run_cmd(cmd: List[str], log_path: Path, env=None) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as logf:
        logf.write(f"[CMD] {shlex_join(cmd)}\n\n")
        logf.flush()
        p = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env)
        return p.wait()

def perchrd(perchr_dir: Path, chr_name: str, suffix: str) -> Path:
    """Prefer {PREFIX}.{CHR}{suffix}; fallback to {CHR}{suffix} if not found."""
    if PREFIX:
        p2 = perchr_dir / f"{PREFIX}.{chr_name}{suffix}"
        if p2.exists():
            return p2
    p1 = perchr_dir / f"{chr_name}{suffix}"
    return p2 if PREFIX and not p1.exists() else p1  # return p2 (may not exist) to keep consistent style, else p1

def derive_pos_from_mafs(mafs_gz: Path, out_pos: Path) -> None:
    with gzip.open(mafs_gz, "rt") as f, open(out_pos, "w") as out:
        header = f.readline().strip().split()
        def idx(name, default):
            try: return header.index(name)
            except Exception: return default
        chrom_i = idx("chromo", 0)
        pos_i   = idx("position", 1)
        for line in f:
            if not line.strip(): continue
            parts = line.split()
            out.write(f"{parts[chrom_i]}\t{parts[pos_i]}\n")

def derive_pos_from_beagle(beagle_gz: Path, out_pos: Path) -> None:
    with gzip.open(beagle_gz, "rt") as f, open(out_pos, "w") as out:
        _ = f.readline()  # header
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split()
            marker = parts[0]
            m = re.match(r'^(\S+)[:_](\d+)$', marker)
            if m:
                chrom, pos = m.group(1), int(m.group(2))
                out.write(f"{chrom}\t{pos}\n")


def split_beagle_to_perchr(beagle_gz: Path, outdir: Path) -> list:
    """
    Split a merged Beagle.gz (marker column formatted as 'CHR:POS' or 'CHR_POS')
    into per-chromosome Beagle.gz files in `outdir` named '{CHR}.beagle.gz'.

    Returns the list of chromosome names encountered (in first-seen order).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    writers = {}  # chr -> gzip text writer
    chrs = []
    with gzip.open(beagle_gz, "rt") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"{beagle_gz} appears empty or lacks a header")
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            marker = parts[0]
            m = re.match(r'^([^:_\s]+)[:_](\d+)$', marker)
            if not m:
                # If marker column isn't in expected form, try last ':' split as fallback
                if ':' in marker:
                    chrom = marker.rsplit(':', 1)[0]
                elif '_' in marker:
                    chrom = marker.rsplit('_', 1)[0]
                else:
                    raise ValueError(f"Cannot parse CHR from marker '{marker}' (expected CHR:POS or CHR_POS)")
            else:
                chrom = m.group(1)
            w = writers.get(chrom)
            if w is None:
                outpath = outdir / f"{chrom}.beagle.gz"
                w = gzip.open(outpath, "wt")
                w.write(header)
                writers[chrom] = w
                chrs.append(chrom)
            writers[chrom].write(line)
    # Close writers
    for w in writers.values():
        w.close()
    return chrs



def split_mafs_to_perchr(mafs_gz: Path, outdir: Path) -> list:
    """
    Split a merged ANGSD .mafs.gz into per-chromosome .mafs.gz files named '{CHR}.mafs.gz'.
    Expects a header with 'chromo' and 'position' columns.
    Returns the list of chromosome names encountered.
    """
    import gzip
    outdir.mkdir(parents=True, exist_ok=True)
    writers = {}  # chr -> gzip text writer
    chrs = []
    with gzip.open(mafs_gz, "rt") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"{mafs_gz} appears empty or lacks a header")
        hdr = header.strip().split()
        # Preferred header names; fallbacks for older outputs
        chrom_keys = ["chromo", "chromosome"]
        pos_keys = ["position", "pos"]
        chrom_i = next((hdr.index(k) for k in chrom_keys if k in hdr), None)
        pos_i = next((hdr.index(k) for k in pos_keys if k in hdr), None)
        if chrom_i is None or pos_i is None:
            raise ValueError(f"Could not find chromo/position columns in {mafs_gz} header: {hdr}")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            chrom = parts[chrom_i]
            g = writers.get(chrom)
            if g is None:
                outpath = outdir / f"{chrom}.mafs.gz"
                g = gzip.open(outpath, "wt")
                g.write(header)
                writers[chrom] = g
                chrs.append(chrom)
            writers[chrom].write(line)
    for g in writers.values():
        g.close()
    return chrs

def ensure_pos(perchr_dir: Path, outdir: Path, chr_name: str, verbose: bool) -> Path:
    """
    Return {outdir}/{CHR}.pos; if missing, derive from MAFS or Beagle and write it in outdir.
    """
    pos_out = outdir / f"{chr_name}.pos"
    if pos_out.exists():
        return pos_out

    # Source: prefer MAFS (perchr_dir, prefix-aware), else Beagle
    beagle = perchrd(perchr_dir, chr_name, ".beagle.gz")
    if beagle.exists():
        if verbose: print(f"[info] deriving {pos_out} from {beagle.name}")
        derive_pos_from_beagle(beagle, pos_out)
        return pos_out

    mafs = perchrd(perchr_dir, chr_name, ".mafs.gz")
    if mafs.exists():
        if verbose: print(f"[info] deriving {pos_out} from {mafs.name}")
        derive_pos_from_mafs(mafs, pos_out)
        return pos_out

    raise FileNotFoundError(f"{chr_name}: missing .pos in outdir and no .mafs.gz/.beagle.gz to derive positions")

def parse_unlinked_positions(pruned_path: Path) -> List[Tuple[str,int]]:
    out = []
    with open(pruned_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"): continue
            parts = s.split()
            if len(parts) >= 2 and parts[1].isdigit():
                out.append((parts[0], int(parts[1]))); continue
            m = re.match(r'^(\S+)[:_ ](\d+)$', s)
            if m: out.append((m.group(1), int(m.group(2)))); continue
            m2 = re.match(r'^(\S+)_(\d+)$', s)
            if m2: out.append((m2.group(1), int(m2.group(2)))); continue
    return out

def write_unlinked_formats(outdir: Path, chr_name: str, sites: List[Tuple[str,int]]) -> Tuple[Path,Path]:
    ids_path = outdir / f"{chr_name}.unlinked_ids.txt"
    tsv_path = outdir / f"{chr_name}.unlinked_sites.tsv"
    with open(ids_path, "w") as f1, open(tsv_path, "w") as f2:
        for chrom, pos in sites:
            f1.write(f"{chrom}_{pos}\n")
            f2.write(f"{chrom}\t{pos}\n")
    return ids_path, tsv_path

def load_mafs_index(mafs_gz: Path) -> Dict[int, Tuple[str,int,str,str]]:
    idx: Dict[int, Tuple[str,int,str,str]] = {}
    with gzip.open(mafs_gz, "rt") as f:
        header = f.readline().strip().split()
        def col(name, default):
            try: return header.index(name)
            except Exception: return default
        chrom_i = col("chromo", 0)
        pos_i   = col("position", 1)
        major_i = col("major",   2)
        minor_i = col("minor",   3)
        for line in f:
            if not line.strip(): continue
            parts = line.split()
            try: pos = int(parts[pos_i])
            except Exception: continue
            idx[pos] = (parts[chrom_i], pos, parts[major_i], parts[minor_i])
    return idx

def write_sites_4col(outdir: Path, chr_name: str, pruned_sites: List[Tuple[str,int]], mafs_gz: Path) -> Path:
    out_path = outdir / f"{chr_name}.sites.tsv"
    index = load_mafs_index(mafs_gz)
    with open(out_path, "w") as out:
        for chrom, pos in pruned_sites:
            if pos in index:
                c, p, maj, mino = index[pos]
                out.write(f"{c}\t{p}\t{maj}\t{mino}\n")
    return out_path

def write_sites_4col_from_beagle(outdir: Path, chr_name: str, pruned_sites: List[Tuple[str,int]], beagle_gz: Path) -> Path:
    out_path = outdir / f"{chr_name}.sites.tsv"
    want = {(c,p) for c,p in pruned_sites}
    with gzip.open(beagle_gz, "rt") as f, open(out_path, "w") as out:
        _ = f.readline()
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            marker, a1, a2 = parts[0], parts[1], parts[2]
            m = re.match(r'^(\S+)[:_](\d+)$', marker)
            if not m: continue
            chrom, pos = m.group(1), int(m.group(2))
            if (chrom, pos) in want:
                out.write(f"{chrom}\t{pos}\t{a1}\t{a2}\n")
    return out_path



def merge_gz_mafs(mafs_paths: List[Path], out_gz: Path) -> Path:
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_gz, "wt") as w:
        header_written = False
        for p in mafs_paths:
            if not p.exists(): continue
            with gzip.open(p, "rt") as r:
                header = r.readline()
                if not header_written:
                    w.write(header); header_written = True
                shutil.copyfileobj(r, w)
    return out_gz

def merge_plain(paths: List[Path], out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as w:
        for p in paths:
            if not p.exists(): continue
            with open(p, "r") as r:
                shutil.copyfileobj(r, w)
    return out_path

def read_pruned_sites_file(pruned_path: Path) -> List[Tuple[str,int]]:
    return parse_unlinked_positions(pruned_path)

def build_union_pruned_sites(chrs: List[str], outdir: Path, external_list_file: Optional[Path]) -> Dict[str, Set[int]]:
    union: Dict[str, Set[int]] = {}
    # Current run
    for c in chrs:
        pr = outdir / f"{c}.ld.pruned"
        if pr.exists():
            for chrom, pos in parse_unlinked_positions(pr):
                union.setdefault(chrom, set()).add(pos)
    # External
    if external_list_file and external_list_file.exists():
        for line in external_list_file.read_text().splitlines():
            fp = line.strip()
            if not fp: continue
            pth = Path(fp)
            if not pth.exists(): continue
            for chrom, pos in parse_unlinked_positions(pth):
                union.setdefault(chrom, set()).add(pos)
    return union

def write_merged_beagle_full(out_beagle_gz: Path, chrs: List[str], perchr_dir: Path) -> Path:
    out_beagle_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_beagle_gz, "wt") as w:
        header_written = False
        for c in chrs:
            bpath = perchrd(perchr_dir, c, ".beagle.gz")
            if not bpath.exists(): continue
            with gzip.open(bpath, "rt") as r:
                header = r.readline()
                if not header_written:
                    w.write(header); header_written = True
                shutil.copyfileobj(r, w)
    return out_beagle_gz

def write_merged_beagle_pruned(out_beagle_gz: Path, chrs: List[str], perchr_dir: Path, union_sites: Dict[str,Set[int]]) -> Path:
    out_beagle_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_beagle_gz, "wt") as w:
        header_written = False
        for c in chrs:
            s = union_sites.get(c, set())
            if not s: continue
            bpath = perchrd(perchr_dir, c, ".beagle.gz")
            if not bpath.exists(): continue
            with gzip.open(bpath, "rt") as r:
                header = r.readline()
                if not header_written:
                    w.write(header); header_written = True
                for line in r:
                    parts = line.strip().split()
                    if not parts: continue
                    marker = parts[0]
                    m = re.match(r'^(\S+)[:_](\d+)$', marker)
                    if not m: continue
                    chrom, pos = m.group(1), int(m.group(2))
                    if pos in union_sites.get(chrom, set()):
                        w.write(line)
    return out_beagle_gz

def write_merged_pruned_mafs(out_mafs_gz: Path, chrs: List[str], perchr_dir: Path, union_sites: Dict[str,Set[int]]) -> Path:
    return _write_merged_pruned_mafs_impl(out_mafs_gz, chrs, perchr_dir, union_sites)

def write_merged_unlinked_lists(outdir: Path, union_sites: Dict[str,Set[int]]):
    ids_path = outdir / 'merged.unlinked_ids.txt'
    tsv_path = outdir / 'merged.unlinked_sites.tsv'
    # Collect and sort
    rows = []
    for chrom in sorted(union_sites.keys(), key=chrom_key):
        for pos in sorted(union_sites[chrom]):
            rows.append((chrom, pos))
    with open(ids_path, 'w') as f1, open(tsv_path, 'w') as f2:
        for chrom, pos in rows:
            f1.write(f"{chrom}_{pos}\n")
            f2.write(f"{chrom}\t{pos}\n")
    return ids_path, tsv_path

def _write_merged_pruned_mafs_impl(out_mafs_gz: Path, chrs: List[str], perchr_dir: Path, union_sites: Dict[str,Set[int]]) -> Path:
    out_mafs_gz.parent.mkdir(parents=True, exist_ok=True)
    header_written = False
    with gzip.open(out_mafs_gz, "wt") as w:
        for c in chrs:
            s = union_sites.get(c, set())
            if not s: continue
            mpath = perchrd(perchr_dir, c, ".mafs.gz")
            if not mpath.exists(): continue
            with gzip.open(mpath, "rt") as r:
                header = r.readline()
                if not header_written:
                    w.write(header); header_written = True
                for line in r:
                    parts = line.strip().split()
                    if not parts: continue
                    chrom = parts[0]
                    try: pos = int(parts[1])
                    except Exception: continue
                    if pos in union_sites.get(chrom, set()):
                        w.write(line)
    return out_mafs_gz

def ngsld_cmd(beagle_gz: Path, pos: Path, out_ld: Path, n_ind: Optional[int], max_kb_dist: int, min_maf: float, threads: int, n_sites: int) -> List[str]:
    cmd = [
        "ngsLD",
        "--geno", str(beagle_gz),
        "--probs", "1",
        "--pos", str(pos),
        "--out", str(out_ld),
        "--max_kb_dist", str(max_kb_dist),
        "--min_maf", str(min_maf),
        "--n_threads", str(threads if threads and threads > 0 else 1),
        "--n_sites", str(n_sites),
    ]
    if n_ind is not None:
        cmd.extend(["--n_ind", str(n_ind)])
    return cmd

def prune_graph_cmd(ld_tsv: Path, out_pruned: Path, min_dist: int, min_r2: float, threads: int, verbose: bool) -> List[str]:
    cmd = [
        "prune_graph",
        "--in", str(ld_tsv),
        "--out", str(out_pruned),
        "--weight-field", "column_7",
        "--weight-filter", f"column_3 > {min_dist} && column_7 > {min_r2}",
    ]
    if threads and threads > 0:
        cmd.extend(["-n", str(threads)])
    if verbose:
        cmd.append("-v")
    return cmd

def step_markers(base: Path, tag: str) -> Tuple[Path, Path]:
    return base.with_suffix(base.suffix + f".{tag}.running"), base.with_suffix(base.suffix + f".{tag}.done")


def discover_chrs(perchr_dir: Path, chr_list_file: Optional[Path], reverse_sort: bool) -> List[str]:
    if chr_list_file:
        chrs = [c.strip() for c in chr_list_file.read_text().splitlines() if c.strip()]
    else:
        chrs = []
        if PREFIX:
            for p in perchr_dir.glob(f"{PREFIX}.*.beagle.gz"):
                base = p.name[:-len(".beagle.gz")]
                if base.startswith(PREFIX + "."):
                    name = base[len(PREFIX)+1:]
                    chrs.append(name)
        if not chrs:
            for p in perchr_dir.glob("*.beagle.gz"):
                base = p.name[:-len(".beagle.gz")]
                if PREFIX and base.startswith(PREFIX + "."):
                    name = base[len(PREFIX)+1:]
                else:
                    name = base
                chrs.append(name)
    if not chrs:
        raise SystemExit("No chromosomes found. Provide -c/--chr-list or place per-chr Beagle files in --perchr-dir.")
    # De-dup + natural numeric sort
    chrs = list(dict.fromkeys(chrs))
    chrs.sort(key=chrom_key, reverse=reverse_sort)
    return chrs


def infer_prefix(perchr_dir: Path, chrs: List[str]) -> str:
    """Infer a prefix by looking for PREFIX.{CHR}.beagle.gz/.mafs.gz/.pos; return '' if none found."""
    candidates = []
    for c in chrs:
        for suf in (".beagle.gz", ".mafs.gz", ".pos"):
            for p in perchr_dir.glob(f"*.{c}{suf}"):
                base = p.name
                cut = base.rfind(f".{c}{suf}")
                pref = base[:cut] if cut > 0 else ""
                if pref:
                    candidates.append(pref)
    if not candidates:
        for p in perchr_dir.glob("*.beagle.gz"):
            base = p.name[:-len(".beagle.gz")]
            cut = base.rfind(".")
            if cut > 0:
                candidates.append(base[:cut])
    if not candidates:
        return ""
    from collections import Counter
    return Counter(candidates).most_common(1)[0][0]

# ---------- Pipeline ----------

def run_ngsld_for_chr(args, chr_name: str, beagle: Path, pos: Path) -> Path:
    ld_out = args.outdir / f"{chr_name}.ld"
    log = args.outdir / "logs" / f"ngsLD_{chr_name}.log"
    running, done = step_markers(ld_out, "ngsLD")

    if ld_out.exists() and not args.overwrite:
        print(f"[skip] ngsLD {chr_name}: exists.")
        return ld_out

    for m in (running, done):
        if m.exists(): m.unlink()
    running.touch()

    # Preflight: verify POS rows match Beagle markers and look sane
    try:
        import gzip
        beagle_rows = 0
        with gzip.open(str(beagle), "rt") as bf:
            header = bf.readline().strip().split() if True else []
            if len(header) >= 6 and (len(header) - 3) % 3 == 0:
                pass
            for _ in bf:
                beagle_rows += 1
        pos_rows = sum(1 for _ in open(pos, "r"))
        if pos_rows != beagle_rows:
            raise RuntimeError(f"[preflight] POS rows ({pos_rows}) != Beagle markers ({beagle_rows}) for {chr_name}. Regenerating POS from Beagle may help.")
    except Exception as e:
        print(str(e))
    pos_rows = sum(1 for _ in open(pos, 'r'))
    cmd = ngsld_cmd(beagle, pos, ld_out, args.nind, args.max_kb_dist, args.min_maf, args.threads, pos_rows)
    print(f"[run] {chr_name} ngsLD -> {ld_out.name}")
    ret = run_cmd(cmd, log)
    if ret != 0 or not ld_out.exists():
        running.unlink(missing_ok=True)
        raise RuntimeError(f"ngsLD failed for {chr_name}. See log: {log}")
    running.unlink(missing_ok=True); done.touch()
    return ld_out


def run_prune_graph_for_chr(args, chr_name: str, ld_tsv: Path) -> Path:
    pruned = args.outdir / f"{chr_name}.ld.pruned"
    log = args.outdir / "logs" / f"prune_graph_{chr_name}.log"
    running, done = step_markers(pruned, "prune")

    if pruned.exists() and not args.overwrite:
        print(f"[skip] prune_graph {chr_name}: exists.")
        return pruned

    for m in (running, done):
        if m.exists(): 
            m.unlink()
    running.touch()

    # Build and run prune_graph without attempting to reference 'beagle' here
    cmd = prune_graph_cmd(ld_tsv, pruned, args.min_dist, args.min_r2, args.threads, True)
    print(f"[run] {chr_name} prune_graph -> {pruned.name}")
    ret = run_cmd(cmd, log)
    if ret != 0 or not pruned.exists():
        running.unlink(missing_ok=True)
        raise RuntimeError(f"prune_graph failed for {chr_name}. See log: {log}")

    running.unlink(missing_ok=True)
    done.touch()
    return pruned

    for m in (running, done):
        if m.exists(): m.unlink()
    running.touch()

    # Preflight: verify POS rows match Beagle markers and look sane
    try:
        import gzip
        beagle_rows = 0
        with gzip.open(str(beagle), "rt") as bf:
            header = bf.readline().strip().split() if True else []
            if len(header) >= 6 and (len(header) - 3) % 3 == 0:
                pass
            for _ in bf:
                beagle_rows += 1
        pos_rows = sum(1 for _ in open(pos, "r"))
        if pos_rows != beagle_rows:
            raise RuntimeError(f"[preflight] POS rows ({pos_rows}) != Beagle markers ({beagle_rows}) for {chr_name}. Regenerating POS from Beagle may help.")
    except Exception as e:
        print(str(e))
    env = os.environ.copy()
    if args.threads and args.threads > 0:
        env["OMP_NUM_THREADS"] = str(args.threads)

    cmd = prune_graph_cmd(ld_tsv, pruned, args.min_dist, args.min_r2, args.threads, True)
    print(f"[run] {chr_name} prune_graph -> {pruned.name}")
    ret = run_cmd(cmd, log, env=env)
    if ret != 0 or not pruned.exists():
        running.unlink(missing_ok=True)
        raise RuntimeError(f"prune_graph failed for {chr_name}. See log: {log}")
    running.unlink(missing_ok=True); done.touch()
    return pruned

def process_chr(args, chr_name: str):
    beagle = perchrd(perchrd_dir, chr_name, ".beagle.gz")
    if not beagle.exists():
        raise FileNotFoundError(f"{beagle} not found")

    pos = ensure_pos(perchrd_dir, args.outdir, chr_name, True)

    # Count number of sites directly from POS (needed for ngsLD bookkeeping)
    n_sites = sum(1 for _ in open(pos, "r"))
    with open(args.outdir / f"{chr_name}.n_sites", "w") as w:
        w.write(str(n_sites) + "\n")
    print(f"[pos] {chr_name}: POS has n_sites={n_sites}")

    ld_tsv = run_ngsld_for_chr(args, chr_name, beagle, pos)
    pruned = run_prune_graph_for_chr(args, chr_name, ld_tsv)

    # Per-chr sites/unlinked outputs
    pruned_sites = parse_unlinked_positions(pruned)
    ids_path, tsv_path = write_unlinked_formats(args.outdir, chr_name, pruned_sites)

    mafs = perchrd(perchrd_dir, chr_name, ".mafs.gz")
    if mafs.exists():
        sites4 = write_sites_4col(args.outdir, chr_name, pruned_sites, mafs)
        print(f"[sites] {chr_name}: wrote {ids_path.name}, {tsv_path.name}, {sites4.name}")
    else:
        # Fallback: Beagle-derived alleles for pruned set
        sites4 = write_sites_4col_from_beagle(args.outdir, chr_name, pruned_sites, beagle)
        print(f"[sites] {chr_name}: wrote {ids_path.name}, {tsv_path.name}, {sites4.name} (from Beagle)")

    print(f"[done] {chr_name}")
    return (chr_name, ld_tsv, pruned, ids_path, tsv_path)

# ---------- Main ----------

def main():
    global perchrd_dir, PREFIX
    ap = argparse.ArgumentParser(description="LD prune ANGSD Beagle (original CLI preserved) with logs, defaults, resume, sites, and merging.")
    ap.add_argument("-c", "--chr-list", help="Chromosome names, one per line")
    ap.add_argument("-p", "--perchr-dir", required=False, help="Dir with per-chr ANGSD {CHR}.mafs.gz / .beagle.gz (or {prefix}.{CHR}.*) [unless -B]")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory (per-chr files & logs)")
    ap.add_argument("-B", "--single-beagle", help="Merged Beagle (plain or .gz) to split into per-chr {CHR}.beagle.gz in --outdir (not implemented here)")
    ap.add_argument("-j", "--jobs", type=int, default=1, help="Parallel jobs for BOTH ngsLD and prune_graph")
    ap.add_argument("-t", "--threads", type=int, default=1, help="Threads per job (ngsLD, prune_graph OMP, bgzip)")
    ap.add_argument("-n", "--nind", type=int, help="ngsLD --n_ind")
    ap.add_argument("-k", "--max-kb-dist", dest="max_kb_dist", type=int, default=1000, help="ngsLD --max_kb_dist")
    ap.add_argument("-m", "--min-maf", dest="min_maf", type=float, default=0.05, help="ngsLD --min_maf")
    ap.add_argument("-r", "--min-r2", dest="min_r2", type=float, default=0.1, help="threshold used in prune_graph weight-filter (column_7 > MIN_R2)")
    ap.add_argument("-d", "--min-dist", dest="min_dist", type=int, default=50000, help="threshold used in prune_graph weight-filter (column_3 > MIN_DIST)")
    ap.add_argument("-x", "--prefix", help="Preferred prefix to match {prefix}.{CHR}.* inputs; also used when writing missing .pos")
    ap.add_argument("-R", "--reverse-sort", action="store_true", help="Reverse version sort of chromosomes")
    ap.add_argument("-g", "--gl-basename", dest="gl_basename", default="merged", help="Merged Beagle basename (no .gz); written in --outdir")
    ap.add_argument("-G", "--pruned-gl-basename", dest="pruned_gl_basename", default="pruned", help="Pruned Beagle basename (no .gz); written in --outdir")
    ap.add_argument("-I", "--keep-intermediates", action="store_true", help="Keep per-chr .pos/.ld/.ld.pruned")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs and rerun steps")
    ap.add_argument("--merge", action="store_true", help="Merge per-chrom outputs (mafs and site lists)")
    ap.add_argument("--external-pruned-list", help="File with full paths to .ld.pruned files from prior runs to include in merged beagle and mafs")

    args = ap.parse_args()
    # verbose default (no flag)
    args.outdir = Path(args.outdir).resolve()
    (args.outdir / "logs").mkdir(parents=True, exist_ok=True)

    # Determine input source: per-chr dir or single merged Beagle
    if args.single_beagle:
        single = Path(args.single_beagle).resolve()
        if not single.exists():
            ap.error(f"--single-beagle file not found: {single}")
        print(f"[input] splitting merged Beagle -> per-chr in {args.outdir}")
        # Create per-chr beagles in outdir and use outdir as our working per-chr dir
        split_chrs = split_beagle_to_perchr(single, args.outdir)

        # If a matching single MAFS exists (same basename, .mafs.gz), split it too
        sname = single.name
        if sname.endswith(".beagle.gz"):
            mafs_single = single.with_name(sname[:-len(".beagle.gz")] + ".mafs.gz")
        else:
            mafs_single = single.with_suffix("").with_suffix(".mafs.gz")
        if mafs_single.exists():
            print(f"[input] found matching MAFS: {mafs_single.name} -> splitting per-chr")
            split_mafs_to_perchr(mafs_single, args.outdir)

        perchrd_dir = args.outdir
        if not split_chrs:
            ap.error("No chromosomes found when splitting the merged Beagle")
        print(f"[input] discovered chromosomes from merged Beagle: {', '.join(split_chrs)}")
    else:
        if not args.perchr_dir:
            ap.error("Provide either -p/--perchr-dir or -B/--single-beagle")
        perchrd_dir = Path(args.perchr_dir).resolve()
        if not perchrd_dir.exists():
            ap.error(f"--perchr-dir not found: {perchrd_dir}")

    # Initial TEMP discovery to help prefix inference
    PREFIX = args.prefix or ""
    chrs = discover_chrs(perchrd_dir, Path(args.chr_list) if args.chr_list else None, args.reverse_sort)
    if not PREFIX:
        inferred = infer_prefix(perchrd_dir, chrs)
        if inferred:
            PREFIX = inferred
            print(f"[prefix] inferred '{PREFIX}' from per-chr files")
        # Refresh chromosome discovery now that PREFIX may be set
        chrs = discover_chrs(perchrd_dir, Path(args.chr_list) if args.chr_list else None, args.reverse_sort)

    print("[config] perchr-dir:", perchrd_dir)
    print("[config] outdir:", args.outdir)
    print("[config] chrs:", ", ".join(chrs))
    print(f"[config] ngsLD: --max_kb_dist {args.max_kb_dist} --min_maf {args.min_maf} --n_threads {args.threads} --n_ind {args.nind}")
    print(f"[config] prune_graph: --in/--out with --weight-filter 'column_3 > {args.min_dist} && column_7 > {args.min_r2}' -n {args.threads}")
    print(f"[config] overwrite={args.overwrite}, keep_intermediates={args.keep_intermediates}")

    # TODO: implement -B split if needed (preserving CLI; not requested to change here)
    if args.single_beagle:
        print("[info] Using per-chr Beagles derived from merged input")

    # Process chromosomes (parallel)
    results = []
    if args.jobs and args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(process_chr, args, c): c for c in chrs}
            for fut in as_completed(futs):
                results.append(fut.result())
    else:
        for c in chrs:
            results.append(process_chr(args, c))

    # Merge artifacts
    # Build union of pruned sites including external pruned list
    union_sites = build_union_pruned_sites(chrs, args.outdir, Path(args.external_pruned_list) if args.external_pruned_list else None)

    # Conditional merged unlinked/site lists (current run + external appended)
    if args.merge:
        mafs_list = [perchrd(perchrd_dir, c, ".mafs.gz") for c in chrs if perchrd(perchrd_dir, c, ".mafs.gz").exists()]
        if mafs_list:
            merge_gz_mafs(mafs_list, args.outdir / "merged.mafs.gz")
            print("[merge] merged.mafs.gz")

        ids_list   = [args.outdir / f"{c}.unlinked_ids.txt" for c in chrs]
        tsv_list   = [args.outdir / f"{c}.unlinked_sites.tsv" for c in chrs]
        sites_list = [args.outdir / f"{c}.sites.tsv" for c in chrs]
        merge_plain([p for p in ids_list   if p.exists()], args.outdir / "merged.unlinked_ids.txt")
        merge_plain([p for p in tsv_list   if p.exists()], args.outdir / "merged.unlinked_sites.tsv")
        merge_plain([p for p in sites_list if p.exists()], args.outdir / "merged.sites.tsv")
        print("[merge] merged unlinked/site files (current run)")

        if args.external_pruned_list and Path(args.external_pruned_list).exists():
            ext_pairs = []
            for line in Path(args.external_pruned_list).read_text().splitlines():
                fp = line.strip()
                if not fp: continue
                for chrom, pos in parse_unlinked_positions(Path(fp)):
                    ext_pairs.append((chrom, pos))
            if ext_pairs:
                with open(args.outdir / "merged.unlinked_ids.txt", "a") as f1, open(args.outdir / "merged.unlinked_sites.tsv", "a") as f2:
                    for chrom, pos in ext_pairs:
                        f1.write(f"{chrom}_{pos}\n")
                        f2.write(f"{chrom}\t{pos}\n")
                print("[merge] appended external pruned to merged unlinked files")

    # Always write merged GLs and pruned MAFS
    full_gl = args.outdir / f"{args.gl_basename}.beagle.gz"
    pruned_gl = args.outdir / f"{args.pruned_gl_basename}.beagle.gz"
    write_merged_beagle_full(full_gl, sorted(chrs, key=chrom_key), perchrd_dir)
    write_merged_beagle_pruned(pruned_gl, sorted(chrs, key=chrom_key), perchrd_dir, union_sites)
    print(f"[merge] wrote merged GLs: {full_gl.name}, {pruned_gl.name}")

    pruned_mafs = args.outdir / f"{args.pruned_gl_basename}.mafs.gz"
    write_merged_pruned_mafs(pruned_mafs, sorted(chrs, key=chrom_key), perchrd_dir, union_sites)
    print(f"[merge] wrote merged pruned MAFS: {pruned_mafs.name}")

    # Write per-chr n_sites sidecar if snps.pos exists
    for c in chrs:
        snps = args.outdir / f"{c}.snps.pos"
        if snps.exists():
            n = sum(1 for _ in open(snps, "r"))
            with open(args.outdir / f"{c}.n_sites", "w") as w:
                w.write(str(n) + "\n")
    # Merge n_sites into a single TSV
    with open(args.outdir / "merged.n_sites.tsv", "w") as w:
        w.write("chr\tn_sites\n")
        for c in chrs:
            ns = args.outdir / f"{c}.n_sites"
            if ns.exists():
                w.write(f"{c}\t{open(ns).read().strip()}\n")

    # Cleanup
    if not args.keep_intermediates:
        for c in chrs:
            for suffix in (".pos",):  # keep original .ld and .ld.pruned from prune_graph
                p = args.outdir / f"{c}{suffix}"
                if p.exists(): p.unlink()

    print("All chromosomes completed.")

if __name__ == "__main__":
    main()
