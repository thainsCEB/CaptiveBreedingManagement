#!/usr/bin/env python3
"""
subset_beagle.py — subset/remove individuals in ANGSD .beagle(.gz) files

Behavior
--------
- Keep/remove by 0-based index (accepts 'indX' or 'X').
- If input has a header (starts with 'marker'), the output header is preserved
  for kept individuals (no renaming/renumbering). Removed individuals' header
  triplets and their 3 GL columns are deleted entirely.
- If input lacks a header, the output does NOT invent one.

Examples
--------
# keep two individuals
python beagle_edit.py -i in.beagle.gz -o out.beagle.gz --keep ind0,ind3

# remove a list from file
python beagle_edit.py -i in.beagle.gz -o out.beagle.gz --remove-file rm.txt
"""

import argparse
import gzip
import io
import os
import sys
from typing import List, Tuple, Optional

# ---------- I/O helpers ----------
def _open_maybe_gz(path: str, mode: str):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def _writable_path(p: Optional[str]):
    if p is None:
        return None
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    return _open_maybe_gz(p, "wb")

def _read_txt_list(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

# ---------- Beagle parsing ----------
def parse_header_info(t: io.TextIOBase) -> Tuple[bool, Optional[List[str]], int]:
    """
    Returns (has_header, header_tokens_or_None, n_individuals).
    Leaves stream positioned after header if present; else at first data line.
    """
    pos = t.tell()
    first = t.readline()
    if not first:
        raise SystemExit("Empty beagle file.")
    toks = first.strip().split()
    if toks and toks[0].lower() == "marker":
        ncols = len(toks)
        if (ncols - 3) % 3 != 0:
            raise SystemExit(f"Header has unexpected column count: {ncols}")
        n = (ncols - 3) // 3
        return True, toks, n
    # no header: rewind and peek data
    t.seek(pos)
    dl = t.readline()
    if not dl:
        raise SystemExit("Beagle has no data.")
    dt = dl.strip().split()
    if (len(dt) - 3) % 3 != 0:
        raise SystemExit(f"Data line has unexpected column count: {len(dt)}")
    n = (len(dt) - 3) // 3
    t.seek(pos)
    return False, None, n

# ---------- Selection helpers ----------
def _to_indices_list(raw: Optional[str], file_path: Optional[str]) -> Optional[List[int]]:
    if raw is None and file_path is None:
        return None
    items: List[str] = []
    if raw:
        items.extend([x.strip() for x in raw.split(",") if x.strip()])
    if file_path:
        items.extend(_read_txt_list(file_path))
    idxs: List[int] = []
    for x in items:
        s = x
        if s.lower().startswith("ind"):
            s = s[3:]
        try:
            i = int(s)
            if i < 0:
                raise ValueError
            idxs.append(i)
        except ValueError:
            raise SystemExit(f"Bad individual identifier: {x} (use indX or X, 0-based)")
    # de-dup preserving order
    out, seen = [], set()
    for i in idxs:
        if i not in seen:
            out.append(i)
            seen.add(i)
    return out

def compute_keep_indices(n: int, keep_idxs: Optional[List[int]], remove_idxs: Optional[List[int]]) -> List[int]:
    all_idxs = list(range(n))
    keep_set = set(all_idxs) if not keep_idxs else set(keep_idxs)
    if remove_idxs:
        keep_set -= set(remove_idxs)
    kept = [i for i in all_idxs if i in keep_set]
    if not kept:
        raise SystemExit("After keep/remove filtering, no individuals remain.")
    return kept

# ---------- Header writer (preserve tokens) ----------
def write_filtered_header_preserving_tokens(out: io.TextIOBase, header_toks: List[str], keep_idx: List[int]) -> None:
    fixed = header_toks[:3]
    kept_triplets: List[str] = []
    for i in keep_idx:
        st = 3 + 3*i
        kept_triplets.extend(header_toks[st:st+3])
    out.write("\t".join(fixed + kept_triplets) + "\n")

# ---------- Filter stream ----------
def stream_filter(in_path: str, out_path: str, keep_idx: List[int]) -> None:
    with _open_maybe_gz(in_path, "rb") as raw_in, _writable_path(out_path) as raw_out:
        tin = io.TextIOWrapper(raw_in, encoding="utf-8", newline="")
        tout = io.TextIOWrapper(raw_out, encoding="utf-8", newline="")
        has_header, header_toks, _ = parse_header_info(tin)
        # Header: preserve tokens for kept indices only (if header exists)
        if has_header and header_toks is not None:
            write_filtered_header_preserving_tokens(tout, header_toks, keep_idx)
        spans = [(3 + 3*i, 3 + 3*i + 3) for i in keep_idx]
        for ln in tin:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            toks = ln.split()
            fixed = toks[:3]
            kept = []
            for st, en in spans:
                kept.extend(toks[st:en])
            tout.write("\t".join(fixed + kept) + "\n")
        tout.flush()

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Subset/remove individuals in ANGSD .beagle(.gz) files (preserve header tokens for kept samples).")
    ap.add_argument("-i", "--in-beagle", required=True, help="Input .beagle or .beagle.gz")
    ap.add_argument("-o", "--out-beagle", required=True, help="Output .beagle(.gz)")

    # keep/remove (indX or X; 0-based)
    ap.add_argument("--keep", help="Comma-separated individuals to KEEP (indX or X).")
    ap.add_argument("--keep-file", help="File with individuals to KEEP (one per line).")
    ap.add_argument("--remove", help="Comma-separated individuals to REMOVE (indX or X).")
    ap.add_argument("--remove-file", help="File with individuals to REMOVE (one per line).")

    args = ap.parse_args()

    # Inspect input
    with _open_maybe_gz(args.in_beagle, "rb") as raw:
        tin = io.TextIOWrapper(raw, encoding="utf-8", newline="")
        _, _, n = parse_header_info(tin)

    keep_idxs = _to_indices_list(args.keep, args.keep_file)
    remove_idxs = _to_indices_list(args.remove, args.remove_file)

    kept = compute_keep_indices(n, keep_idxs, remove_idxs)
    stream_filter(args.in_beagle, args.out_beagle, kept)

if __name__ == "__main__":
    main()
