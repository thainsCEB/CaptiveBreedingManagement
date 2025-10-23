#!/usr/bin/env python3
from pathlib import Path  # added by patch
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

def merge_q(ids, out_suffix, outdir, merged_q=None, taxa_labels=None):
    """
    Merge per-sample .qopt files into a single merged Q file (NO header).

    Returns:
      rows: List[Tuple[str, List[float]]] where each item is (sample_id, values[K])
      K: int
      merged_q: output file path
    """
    import os, glob, re, sys

    if merged_q is None:
        merged_q = os.path.join(outdir, "merged.Q")
    os.makedirs(outdir, exist_ok=True)

    def _is_float(x):
        try:
            float(x); return True
        except Exception:
            return False

    def parse_qopt_first_row(qfile):
        """Return list[float] from the first usable data row in qfile (skip header-like lines)."""
        with open(qfile, "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if not lines:
            return None
        def _is_header(tokens):
            return (tokens and tokens[0].upper() == "ID") or all(t.upper().startswith("K") for t in tokens)
        data_lines = [ln for ln in lines if not _is_header(ln.split())]
        if not data_lines:
            return None
        for ln in data_lines:
            toks = ln.split()
            if toks and not _is_float(toks[0]) and any(_is_float(x) for x in toks[1:]):
                toks = toks[1:]
            nums = [float(x) for x in toks if _is_float(x)]
            if nums:
                return nums
        return None

    def find_qopt_candidates_for_id(sid):
        suf = out_suffix or ""
        patterns = [
            os.path.join(outdir, f"{sid}*{suf}*.qopt"),
            os.path.join(outdir, f"{sid}*.{suf}.qopt"),
            os.path.join(outdir, f"{sid}{suf}.qopt"),
            os.path.join(outdir, f"{sid}*.qopt"),
            os.path.join(outdir, "**", f"{sid}*{suf}*.qopt"),
            os.path.join(outdir, "**", f"{sid}*.{suf}.qopt"),
            os.path.join(outdir, "**", f"{sid}{suf}.qopt"),
            os.path.join(outdir, "**", f"{sid}*.qopt"),
        ]
        hits = []
        for p in patterns:
            hits.extend(glob.glob(p, recursive=True))
            if p.lower().endswith(".qopt"):
                hits.extend(glob.glob(p[:-5] + ".QOPT", recursive=True))
        seen, uniq = set(), []
        for h in hits:
            if h not in seen and os.path.isfile(h):
                seen.add(h); uniq.append(h)
        return uniq

    def find_all_qopt(outdir):
        pats = [
            os.path.join(outdir, "**", "*.qopt"),
            os.path.join(outdir, "**", "*.QOPT"),
            os.path.join(outdir, "*.qopt"),
            os.path.join(outdir, "*.QOPT"),
        ]
        hits = []
        for p in pats:
            hits.extend(glob.glob(p, recursive=True))
        return [h for h in hits if os.path.isfile(h)]

    # Build mapping list[(sid, qfile)]
    mapping = []
    matched_any = False
    for sid in ids:
        cand = find_qopt_candidates_for_id(sid)
        if cand:
            matched_any = True
            mapping.append((sid, cand[0]))

    if not matched_any:
        # Fallback: use all qopt files; infer sid from basename up to first '.'
        all_q = sorted(find_all_qopt(outdir))
        for q in all_q:
            base = os.path.basename(q)
            sid = base.split('.')[0] if '.' in base else os.path.splitext(base)[0]
            mapping.append((sid, q))

    # Determine K
    K = None
    for _, qf in mapping:
        vals = parse_qopt_first_row(qf)
        if vals:
            K = len(vals)
            break
    if K is None and out_suffix:
        m = re.search(r'[\._-]K(\d+)[\._-]?', str(out_suffix))
        if m:
            try:
                K = int(m.group(1))
            except Exception:
                K = None
    if K is None:
        # Ensure merged file exists (empty), then raise
        open(merged_q, "w").close()
        inv_lines = []
        for root, _, files in os.walk(outdir):
            for fn in files:
                inv_lines.append(os.path.join(root, fn))
                if len(inv_lines) >= 30: break
            if len(inv_lines) >= 30: break
        inv = "\n".join(inv_lines) if inv_lines else "(no files found)"
        msg = (
            "No usable .qopt files found to determine K.\n"
            f"Checked outdir='{outdir}' with out_suffix='{out_suffix}'.\n"
            "Example of files present under outdir:\n"
            f"{inv}"
        )
        raise FileNotFoundError(msg)

    # Build rows (sid, vals) and write merged file (NO header)
    rows = []
    with open(merged_q, "w") as out:
        for sid, qf in mapping:
            vals = parse_qopt_first_row(qf)
            if not vals or len(vals) != K:
                continue
            rows.append((sid, vals))
            out.write(sid + "\t" + "\t".join(f"{v:.6f}" for v in vals) + "\n")

    # Safety: strip accidental header-like first line
    try:
        with open(merged_q, "r") as _f:
            _lines = _f.readlines()
        if _lines:
            first = _lines[0].strip().split()
            if (first and first[0].upper() == "ID") or all(tok.upper().startswith("K") for tok in first):
                with open(merged_q, "w") as _f:
                    _f.writelines(_lines[1:])
    except Exception:
        pass

    return rows, K, merged_q

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
        # Create labeled versions of every *.qopt if taxa_map is provided (robust mapping)
        if args.taxa_map:
            _tlabs, _tmap, _tpairs = None, None, None
            try:
                _tlabs = read_taxa_map(args.taxa_map, K_expected=K)
            except Exception:
                pass
            try:
                _labs2, _tmap, _tpairs = read_taxa_map2(args.taxa_map)
                if _tlabs is None:
                    _tlabs = _labs2
            except Exception:
                _tmap, _tpairs = None, None
            if _tlabs:
                label_all_qopts_in_dir(args.outdir, _tlabs, taxa_map=_tmap, taxa_pairs=_tpairs, sep=" ")

    

# Create labeled versions of every *.qopt if taxa_map is provided (robust mapping)
# Create labeled versions of every *.qopt if taxa_map is provided (robust mapping)
# === Added by patch: label qopt files with taxa/population names ==================
def _infer_has_header_qopt(lines):
    """Return True if first non-empty, non-comment line looks like a header (K1 K2 ...)."""
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        toks = re.split(r"\s+", s)
        # header if tokens like K1, K2 ... and not numeric
        nonnum = any(not re.match(r"^[0-9.+-eE]+$", t) for t in toks)
        k_style = all(re.match(r"^K\d+$", t) for t in toks)
        return nonnum and k_style
    return False

def label_qopt_file(path, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    """
    Read a fastNGSadmix .qopt file and write a labeled copy replacing K1..Kk with taxa_labels.
    Output: <basename>.labeled.qopt in the same directory.
    """
    p = Path(path)
    data = p.read_text().splitlines()
    has_header = _infer_has_header_qopt(data)
    out = []
    replaced = False
    wrote_header = False
    for ln in data:
        s = ln.strip()
        if not s or s.startswith("#"):
            out.append(ln)
            continue
        if has_header and not replaced:
            toks = re.split(r"\s+", s)
            # Replace header in place
            out.append("\t".join(taxa_labels))
            replaced = True
            continue
        else:
            # If no header, we'll prepend it once before the first data row
            if not wrote_header:
                out.append("\t".join(taxa_labels))
                wrote_header = True
            out.append(ln)
    # Edge case: empty file
    if not out and taxa_labels:
        out = ["\t".join(taxa_labels)]
    outp = p.with_suffix(p.suffix + ".labeled")
    Path(outp).write_text("\n".join(out) + "\n")
    return str(outp)

def label_all_qopts_in_dir(outdir, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    """
    Label all *.qopt files in outdir and its immediate subdirs (common when per-sample outputs are nested).
    """
    outdir = Path(outdir)
    made = []
    for root, dirs, files in os.walk(outdir):
        for fn in files:
            if fn.endswith(".qopt"):
                q = Path(root) / fn
                try:
                    outp = label_qopt_file(q, taxa_labels, taxa_map=taxa_map, taxa_pairs=taxa_pairs, sep=sep)
                    print(f"[labeled] {q} -> {outp}")
                    made.append(outp)
                except Exception as e:
                    print(f"[warn] failed to label {q}: {e}")
        # Do not recurse deeply beyond 1 level by default; fastNGSadmix usually flat
        # If deep recursion is needed, rem# === qopt labeling: robust helpers (added) =======================================
import re as _re
from pathlib import Path as _Path
import os as _os

def _tokenize_ws(s: str):
    if s is None: return []
    s = s.replace("\u00A0"," ").replace("\u2007"," ").replace("\u202F"," ")
    s = _re.sub(r"\s+", " ", s.strip())
    return s.split(" ") if s else []

def _looks_like_k_header(tokens):
    if not tokens: return False
    nonnum = any(not _re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(_re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and kstyle

def _is_non_numeric_header(tokens):
    if not tokens: return False
    nonnum = any(not _re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(_re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and not kstyle

def _is_k_header_line(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#"): return False
    return _re.fullmatch(r'(?:\s*K\d+\s+)+K\d+\s*', s) is not None

def read_taxa_map2(path: str):
    labels, mapping, multicol = [], {}, False
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"): continue
            parts = _tokenize_ws(s)
            if len(parts) == 1 and not multicol:
                labels.append(parts[0])
            elif len(parts) >= 2:
                multicol = True
                old, new = parts[0], parts[-1]
                mapping[old] = new
                labels.append(new)
    if multicol and not mapping: raise SystemExit(f"[error] No valid (old,new) pairs in taxa map: {path}")
    if not multicol and not labels: raise SystemExit(f"[error] No taxa labels found in: {path}")
    return labels, (mapping if multicol else None), (list(mapping.items()) if multicol else None)

def label_qopt_file(path, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    pth = _Path(path)
    text = pth.read_text(encoding="utf-8")
    lines = text.splitlines()

    # Remove K-headers early
    lines = [ln for ln in lines if not _is_k_header_line(ln) and not _looks_like_k_header(_tokenize_ws(ln.strip()))]

    # Whole-file replacements using mapping (sed-style word boundaries)
    if taxa_pairs:
        def repl_line(ln):
            out = ln
            for old, new in taxa_pairs:
                pat = r'(?<!\S)' + _re.escape(old) + r'(?!\S)'
                out = _re.sub(pat, new, out)
            return out
        lines = [repl_line(ln) for ln in lines]

    # Remove K-headers again (post-replacement)
    lines = [ln for ln in lines if not _is_k_header_line(ln) and not _looks_like_k_header(_tokenize_ws(ln.strip()))]

    # Detect existing header (first non-empty non-comment, non-numeric-only line)
    first_idx = None
    for idx, ln in enumerate(lines):
        s = ln.strip()
        if not s or s.startswith("#"): continue
        first_idx = idx; break

    existing_header_idx, existing_tokens = None, None
    if first_idx is not None:
        toks = _tokenize_ws(lines[first_idx].strip())
        if _is_non_numeric_header(toks):
            existing_header_idx, existing_tokens = first_idx, toks

    # Build replacement header labels
    rep_labels = list(taxa_labels)
    if taxa_map is not None:
        if existing_tokens:
            rep_labels = [taxa_map.get(tok, tok) for tok in existing_tokens]
        else:
            rep_labels = list(taxa_labels)

    # Compose output: single header + body
    body = [ln for i, ln in enumerate(lines) if i != existing_header_idx]
    out_lines = [sep.join(rep_labels)] + body

    # Trim any second header line immediately following
    trimmed, header_seen = [], False
    for ln in out_lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            trimmed.append(ln); continue
        toks = _tokenize_ws(s)
        if not header_seen:
            trimmed.append(ln); header_seen = True
        else:
            if _is_k_header_line(ln) or _is_non_numeric_header(toks):
                continue
            trimmed.append(ln)
    out_lines = trimmed

    outp = pth
    outp.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return str(outp)

def label_all_qopts_in_dir(outdir, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    outdir = _Path(outdir)
    made = []
    for root, _, files in _os.walk(outdir):
        for fn in files:
            if fn.endswith(".qopt"):
                q = _Path(root) / fn
                try:
                    outp = label_qopt_file(q, taxa_labels, taxa_map=taxa_map, taxa_pairs=taxa_pairs, sep=sep)
                    print(f"[labeled] {q} -> {outp}")
                    made.append(outp)
                except Exception as e:
                    print(f"[warn] failed to label {q}: {e}")
    return made
# ================================================================================
# remove this break.
        # break
    return made
# ==================================================================================

if __name__ == "__main__":
    main()


# === qopt labeling: robust helpers (added) =======================================
import re as _re
from pathlib import Path as _Path
import os as _os

def _tokenize_ws(s: str):
    if s is None: return []
    s = s.replace("\u00A0"," ").replace("\u2007"," ").replace("\u202F"," ")
    s = _re.sub(r"\s+", " ", s.strip())
    return s.split(" ") if s else []

def _looks_like_k_header(tokens):
    if not tokens: return False
    nonnum = any(not _re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(_re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and kstyle

def _is_non_numeric_header(tokens):
    if not tokens: return False
    nonnum = any(not _re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(_re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and not kstyle

def _is_k_header_line(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#"): return False
    return _re.fullmatch(r'(?:\s*K\d+\s+)+K\d+\s*', s) is not None

def read_taxa_map2(path: str):
    labels, mapping, multicol = [], {}, False
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"): continue
            parts = _tokenize_ws(s)
            if len(parts) == 1 and not multicol:
                labels.append(parts[0])
            elif len(parts) >= 2:
                multicol = True
                old, new = parts[0], parts[-1]
                mapping[old] = new
                labels.append(new)
    if multicol and not mapping: raise SystemExit(f"[error] No valid (old,new) pairs in taxa map: {path}")
    if not multicol and not labels: raise SystemExit(f"[error] No taxa labels found in: {path}")
    return labels, (mapping if multicol else None), (list(mapping.items()) if multicol else None)

def label_qopt_file(path, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    pth = _Path(path)
    text = pth.read_text(encoding="utf-8")
    lines = text.splitlines()

    # Remove K-headers early
    lines = [ln for ln in lines if not _is_k_header_line(ln) and not _looks_like_k_header(_tokenize_ws(ln.strip()))]

    # Whole-file replacements using mapping (sed-style word boundaries)
    if taxa_pairs:
        def repl_line(ln):
            out = ln
            for old, new in taxa_pairs:
                pat = r'(?<!\S)' + _re.escape(old) + r'(?!\S)'
                out = _re.sub(pat, new, out)
            return out
        lines = [repl_line(ln) for ln in lines]

    # Remove K-headers again (post-replacement)
    lines = [ln for ln in lines if not _is_k_header_line(ln) and not _looks_like_k_header(_tokenize_ws(ln.strip()))]

    # Detect existing header (first non-empty non-comment, non-numeric-only line)
    first_idx = None
    for idx, ln in enumerate(lines):
        s = ln.strip()
        if not s or s.startswith("#"): continue
        first_idx = idx; break

    existing_header_idx, existing_tokens = None, None
    if first_idx is not None:
        toks = _tokenize_ws(lines[first_idx].strip())
        if _is_non_numeric_header(toks):
            existing_header_idx, existing_tokens = first_idx, toks

    # Build replacement header labels
    rep_labels = list(taxa_labels)
    if taxa_map is not None:
        if existing_tokens:
            rep_labels = [taxa_map.get(tok, tok) for tok in existing_tokens]
        else:
            rep_labels = list(taxa_labels)

    # Compose output: single header + body
    body = [ln for i, ln in enumerate(lines) if i != existing_header_idx]
    out_lines = [sep.join(rep_labels)] + body

    # Trim any second header line immediately following
    trimmed, header_seen = [], False
    for ln in out_lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            trimmed.append(ln); continue
        toks = _tokenize_ws(s)
        if not header_seen:
            trimmed.append(ln); header_seen = True
        else:
            if _is_k_header_line(ln) or _is_non_numeric_header(toks):
                continue
            trimmed.append(ln)
    out_lines = trimmed

    outp = pth
    outp.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return str(outp)

def label_all_qopts_in_dir(outdir, taxa_labels, taxa_map=None, taxa_pairs=None, sep=" "):
    outdir = _Path(outdir)
    made = []
    for root, _, files in _os.walk(outdir):
        for fn in files:
            if fn.endswith(".qopt"):
                q = _Path(root) / fn
                try:
                    outp = label_qopt_file(q, taxa_labels, taxa_map=taxa_map, taxa_pairs=taxa_pairs, sep=sep)
                    print(f"[labeled] {q} -> {outp}")
                    made.append(outp)
                except Exception as e:
                    print(f"[warn] failed to label {q}: {e}")
    return made
# ================================================================================
