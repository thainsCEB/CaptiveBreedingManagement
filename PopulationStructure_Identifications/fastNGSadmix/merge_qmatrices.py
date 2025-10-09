#!/usr/bin/env python3
"""
merge_qmatrices.py — Merge a reference Q matrix (NGSadmix/ADMIXTURE/etc.) with a target Q matrix (fastNGSadmix)
using a single metadata file, optional headerless mode, and a cluster-key mapping (inline or 2-col file).

Features
- ONE metadata file with columns: ID, Source (Reference/Target or ref/tgt), Group (optional).
- Preserve final row order exactly as metadata order.
- --cluster-key accepts inline "K1=Name,..." OR a two-column TSV/CSV mapping (Kindex, Label).
- --lenient-counts trims per-Source to smaller of metadata vs Q if counts mismatch.
- --precision controls decimal places in output.
- --no-header treats metadata as having no header: col1=ID, col2=Source, col3=Group (optional).

Author: ChatGPT
"""

import argparse
import sys
import os
from typing import Dict, List
import pandas as pd
import numpy as np
import re

def read_q_matrix(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=r"\\s+", header=None, comment=None, engine="python")
    except Exception as e:
        raise ValueError(f"Failed to read Q matrix '{path}': {e}")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="raise")
    return df

def normalize_source(val: str) -> str:
    v = str(val).strip().lower()
    if v in ("reference", "ref", "r"):
        return "Reference"
    if v in ("target", "tgt", "t"):
        return "Target"
    raise ValueError(f"Invalid Source value '{val}'. Expected Reference/Target (or ref/tgt).")

def read_metadata(path: str, no_header: bool = False) -> pd.DataFrame:
    """
    Read a forgiving metadata file (TSV/CSV/whitespace). Accepts ragged rows.
    """
    if no_header:
        last_err = None
        for sep in ["\\t", ",", r"\\s+"]:
            try:
                df = pd.read_csv(path, sep=sep, engine="python", header=None, dtype=str, comment="#", skip_blank_lines=True)
                if df.shape[0] == 0:
                    continue
                break
            except Exception as e:
                last_err = e
                df = None
        if df is None:
            raise ValueError(f"Failed to read metadata '{path}' without header: {last_err}")
        if df.shape[1] < 2:
            raise ValueError("Metadata without header requires at least 2 columns: ID, Source (Group optional as col3).")
        out = pd.DataFrame()
        out["ID"] = df.iloc[:, 0].astype(str)
        out["Source"] = df.iloc[:, 1].apply(normalize_source)
        out["Group"] = df.iloc[:, 2].astype(str) if df.shape[1] >= 3 else np.nan
        return out

    last_err = None
    for sep in ["\\t", ",", r"\\s+"]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python", dtype=str, comment="#", skip_blank_lines=True)
            if df.shape[0] == 0:
                continue
            break
        except Exception as e:
            last_err = e
            df = None
    if df is None:
        # Final fallback: manual tolerant parsing
        try:
            rows = []
            with open(path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if "\\t" in line:
                        parts = line.split("\\t")
                    elif "," in line:
                        parts = [x.strip() for x in line.split(",")]
                    else:
                        parts = [x for x in re.split(r"\\s+", line) if x != ""]
                    rows.append(parts)
            if not rows:
                raise ValueError("No rows found in metadata.")
            maxlen = max(len(r) for r in rows)
            rows = [r + [None] * (maxlen - len(r)) for r in rows]
            df = pd.DataFrame(rows)
        except Exception as e:
            raise ValueError(f"Failed to read metadata '{path}': {last_err or e}")

    # Header detection
    if not any(str(c).lower() in ("id","source","group") for c in df.columns):
        head = [str(x).lower() for x in df.iloc[0].tolist()]
        if any(tok in head for tok in ("id","source","group")):
            df.columns = df.iloc[0].tolist()
            df = df.iloc[1:].reset_index(drop=True)
        else:
            df.columns = [f"col{i+1}" for i in range(df.shape[1])]

    cols_lower = [str(c).lower() for c in df.columns]
    id_col = df.columns[cols_lower.index("id")] if "id" in cols_lower else df.columns[0]
    source_col = df.columns[cols_lower.index("source")] if "source" in cols_lower else (df.columns[1] if df.shape[1] >= 2 else df.columns[0])
    group_col = df.columns[cols_lower.index("group")] if "group" in cols_lower else (df.columns[2] if df.shape[1] >= 3 else None)

    out = pd.DataFrame()
    out["ID"] = df[id_col].astype(str)
    out["Source"] = df[source_col].apply(normalize_source)
    out["Group"] = df[group_col].astype(str) if group_col is not None else np.nan
    return out

def _parse_cluster_inline(s: str, k: int) -> List[str]:
    mapping: Dict[int, str] = {}
    parts = re.split(r"[,\\s]+", s.strip())
    for p in parts:
        if p == "":
            continue
        if "=" not in p:
            raise ValueError(f"Bad cluster-key token '{p}'. Expected form K1=Name or 1=Name.")
        left, right = p.split("=", 1)
        m = re.fullmatch(r"[Kk]?(\\d+)", left.strip())
        if not m:
            raise ValueError(f"Bad cluster-key left side '{left}'. Use K1, k2, or just 1.")
        idx = int(m.group(1))
        if not (1 <= idx <= k):
            raise ValueError(f"Cluster index {idx} out of range 1..{k}.")
        mapping[idx] = right.strip() if right.strip() else f"K{idx}"
    return [mapping.get(i+1, f"K{i+1}") for i in range(k)]

def _parse_cluster_file(path: str, k: int) -> List[str]:
    try:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    except Exception as e:
        raise ValueError(f"Failed to read cluster-key file '{path}': {e}")
    if df.shape[1] < 2:
        raise ValueError("Cluster-key file must have at least two columns: <Kindex> <Label>.")
    df = df.iloc[:, :2].copy()
    df.columns = ["K", "Label"]
    mapping: Dict[int, str] = {}
    for _, row in df.iterrows():
        raw = str(row["K"]).strip()
        m = re.fullmatch(r"[Kk]?(\\d+)", raw)
        if not m:
            raise ValueError(f"Invalid cluster index '{raw}' in cluster-key file; expected like 'K1' or '1'.")
        idx = int(m.group(1))
        if not (1 <= idx <= k):
            raise ValueError(f"Cluster index {idx} in cluster-key file out of range 1..{k}.")
        mapping[idx] = str(row["Label"]).strip()
    return [mapping.get(i+1, f"K{i+1}") for i in range(k)]

def parse_cluster_key(s: str | None, k: int) -> List[str]:
    if not s or s.strip() == "":
        return [f"K{i+1}" for i in range(k)]
    if os.path.exists(s) and os.path.isfile(s):
        return _parse_cluster_file(s, k)
    return _parse_cluster_inline(s, k)

def assemble_from_metadata(meta: pd.DataFrame, ref_q: pd.DataFrame, tgt_q: pd.DataFrame, cluster_names: List[str], lenient: bool = False) -> pd.DataFrame:
    n_ref_meta = int((meta["Source"] == "Reference").sum())
    n_tgt_meta = int((meta["Source"] == "Target").sum())
    n_ref_q = int(ref_q.shape[0])
    n_tgt_q = int(tgt_q.shape[0])
    if (n_ref_meta != n_ref_q or n_tgt_meta != n_tgt_q) and not lenient:
        raise ValueError(f"Count mismatch: metadata Reference={n_ref_meta}, Target={n_tgt_meta}; Q Reference={n_ref_q}, Target={n_tgt_q}. Use --lenient-counts to trim.")
    if lenient:
        ref_meta_all = meta[meta["Source"] == "Reference"].copy()
        tgt_meta_all = meta[meta["Source"] == "Target"].copy()
        n_ref = min(n_ref_meta, n_ref_q)
        n_tgt = min(n_tgt_meta, n_tgt_q)
        ref_meta = ref_meta_all.iloc[:n_ref].copy()
        tgt_meta = tgt_meta_all.iloc[:n_tgt].copy()
        ref_q = ref_q.iloc[:n_ref, :].reset_index(drop=True)
        tgt_q = tgt_q.iloc[:n_tgt, :].reset_index(drop=True)
        ref_taken = tgt_taken = 0
        rows_meta = []
        for _, r in meta.iterrows():
            if r["Source"] == "Reference" and ref_taken < n_ref:
                rows_meta.append(r); ref_taken += 1
            elif r["Source"] == "Target" and tgt_taken < n_tgt:
                rows_meta.append(r); tgt_taken += 1
        meta = pd.DataFrame(rows_meta).reset_index(drop=True)
    else:
        ref_meta = meta[meta["Source"] == "Reference"].copy()
        tgt_meta = meta[meta["Source"] == "Target"].copy()

    k = ref_q.shape[1]
    if tgt_q.shape[1] != k:
        raise ValueError(f"K mismatch: ref has {k} columns, target has {tgt_q.shape[1]}.")

    ref_block = ref_q.copy(); ref_block.columns = cluster_names
    tgt_block = tgt_q.copy(); tgt_block.columns = cluster_names

    ref_block.insert(0, "ID", ref_meta["ID"].values)
    ref_block.insert(1, "Group", ref_meta["Group"].values)
    ref_block.insert(2, "Source", "Reference")

    tgt_block.insert(0, "ID", tgt_meta["ID"].values)
    tgt_block.insert(1, "Group", tgt_meta["Group"].values)
    tgt_block.insert(2, "Source", "Target")

    rows = []
    ref_pointer = tgt_pointer = 0
    for _, row in meta.iterrows():
        if row["Source"] == "Reference":
            rows.append(ref_block.iloc[[ref_pointer]]); ref_pointer += 1
        else:
            rows.append(tgt_block.iloc[[tgt_pointer]]); tgt_pointer += 1
    merged = pd.concat(rows, axis=0, ignore_index=True)
    return merged

def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge reference and target Q matrices using metadata and a cluster-key mapping.")
    p.add_argument("-R", "--ref-q", required=True, help="Reference Q matrix file (whitespace-delimited, no header).")
    p.add_argument("-T", "--tgt-q", required=True, help="Target Q matrix file (whitespace-delimited, no header).")
    p.add_argument("-m", "--metadata", required=True, help="Metadata TSV/CSV with columns: ID, Source, Group.")
    p.add_argument("--no-header", action="store_true", help="Treat metadata as having NO header; columns are: 1=ID, 2=Source, 3=Group(optional).")
    p.add_argument("-k", "--cluster-key", default=None, help='Inline "K1=A,..." OR path to two-column TSV/CSV mapping (Kindex, Label).')
    p.add_argument("--lenient-counts", action="store_true", help="Trim to the min(metadata, Q) per Source instead of erroring on count mismatches.")
    p.add_argument("-p", "--precision", type=int, default=4, help="Number of decimal places for Q values in output (default: 4).")
    p.add_argument("-o", "--out", required=True, help="Output TSV.")
    return p.parse_args(argv)

def main(argv: List[str]) -> None:
    args = parse_args(argv)
    ref_q = read_q_matrix(args.ref_q)
    tgt_q = read_q_matrix(args.tgt_q)
    k = ref_q.shape[1]
    if tgt_q.shape[1] != k:
        raise SystemExit(f"K mismatch: ref has {k} columns, target has {tgt_q.shape[1]}.")
    meta = read_metadata(args.metadata, no_header=args.no_header)
    cluster_names = parse_cluster_key(args.cluster_key, k)
    merged = assemble_from_metadata(meta, ref_q, tgt_q, cluster_names, lenient=args.lenient_counts)
    float_fmt = f"% .{args.precision}f".replace(" %", "%")
    outdir = os.path.dirname(os.path.abspath(args.out)) or "."
    os.makedirs(outdir, exist_ok=True)
    merged.to_csv(args.out, sep="\\t", index=False, float_format=float_fmt)
    sys.stderr.write(f"Wrote {merged.shape[0]} rows with K={k} to {args.out} using precision={args.precision}\\n")

if __name__ == "__main__":
    main(sys.argv[1:])
