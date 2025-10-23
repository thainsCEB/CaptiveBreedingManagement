#!/usr/bin/env python3
"""
assign_cluster: (updated) Adds short-hand flags for convenience.

Author: Taylor Hains
Date: 2025-10-13
"""
import argparse, sys, csv, os, re
from pathlib import Path
from typing import List, Optional, Tuple

def _parse_cluster_names(arg: Optional[str], K: int) -> List[str]:
    """
    Accepts:
      - "A,B,C" (comma-separated, length K)
      - "1=A,pop2=B,K3=C" (key=value pairs with 1..K or popN/KNN)
      - path to file:
          * one name per line (length K)
          * two columns: key \t name where key is index (1..K or popN/KNN)
    Returns list[str] length K.
    """
    if not arg:
        return [f"C{i+1}" for i in range(K)]
    # File?
    if os.path.exists(arg):
        mapping = {}
        with open(arg) as f:
            for raw in f:
                s = raw.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) == 1:
                    mapping[len(mapping)+1] = parts[0]
                else:
                    key, name = parts[0].strip(), parts[1].strip()
                    m = re.match(r"^(?:pop|k)?(\d+)$", key, flags=re.IGNORECASE)
                    if not m:
                        raise SystemExit(f"[ERROR] Bad key in cluster-names file: '{key}'")
                    mapping[int(m.group(1))] = name
        return [mapping.get(i+1, f"C{i+1}") for i in range(K)]
    # key=value pairs?
    if "=" in arg:
        mapping = {}
        for tok in re.split(r"\s*,\s*", arg.strip()):
            if not tok: 
                continue
            if "=" not in tok:
                raise SystemExit(f"[ERROR] Expected key=value in --cluster-names token '{tok}'")
            key, name = tok.split("=",1)
            m = re.match(r"^(?:pop|k)?(\d+)$", key.strip(), flags=re.IGNORECASE)
            if not m:
                raise SystemExit(f"[ERROR] Bad key '{key}' in --cluster-names (use 1..K or popN/KNN)")
            mapping[int(m.group(1))] = name.strip()
        return [mapping.get(i+1, f"C{i+1}") for i in range(K)]
    # Otherwise comma-separated list
    names = [s.strip() for s in arg.split(",") if s.strip()]
    if len(names) != K:
        raise SystemExit(f"[ERROR] --cluster-names expects {K} items, got {len(names)}")
    return names
import numpy as np
import pandas as pd

def read_table(path: str) -> pd.DataFrame:
    """
    Read with a permissive parser:
      - Treat any run of whitespace as a delimiter.
      - Also accept CSV/TSV by letting pandas split on whitespace after on-the-fly pre-clean.
    """
    # Fast path: whitespace-delimited
    try:
        df = pd.read_table(path, sep=r"\s+", header=None, engine="python", comment="#")
        return df
    except Exception:
        # Fallback: standard CSV (commas/tabs)
        return pd.read_csv(path, header=None, comment="#")

def coerce_numeric_frame(df: pd.DataFrame) -> Tuple[Optional[pd.Series], pd.DataFrame, List[int]]:
    """
    Decide if first column is IDs (non-numeric across most rows). If so, peel it off.
    Convert remaining to numeric, drop all-NaN columns (trailing delimiters), return (ids, q, dropped_cols).
    """
    nrows, ncols = df.shape
    ids = None
    first_col_nonnum_ratio = None

    # Heuristic: if >90% of values in col0 are non-numeric, treat as IDs
    col0 = df.iloc[:,0].astype(str)
    def is_floatish(x: str) -> bool:
        try:
            float(x)
            return True
        except Exception:
            return False
    first_col_nonnum_ratio = 1.0 - (col0.map(is_floatish).mean())

    if first_col_nonnum_ratio > 0.9:
        ids = col0
        num = df.iloc[:,1:].copy()
    else:
        num = df.copy()

    # Coerce to numeric
    for c in num.columns:
        num[c] = pd.to_numeric(num[c], errors="coerce")

    # Drop columns that are entirely NaN (e.g., from trailing whitespace/delimiters)
    before_cols = num.shape[1]
    keep_mask = ~num.isna().all(axis=0)
    dropped_cols = [int(c) for c in num.columns[~keep_mask]]
    num = num.loc[:, keep_mask]

    # If any row is entirely NaN, that means parsing failed badly
    if (num.notna().sum(axis=1) == 0).any():
        bad = int((num.notna().sum(axis=1) == 0).sum())
        raise SystemExit(f"[error] {bad} rows have no numeric values after parsing. Check the Q file formatting.")

    # Fill any remaining NaNs with 0 and renormalize rows to sum~1 (if sums>0)
    arr = num.to_numpy(dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    row_sums = arr.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        norm = np.where(row_sums[:,None] > 0, arr / row_sums[:,None], arr)
    num = pd.DataFrame(norm, index=num.index, columns=range(norm.shape[1]))

    return (ids, num, dropped_cols)

def read_inputs(q_path: str, q_has_ids: bool, ids_path: Optional[str]) -> Tuple[List[str], pd.DataFrame]:
    raw = read_table(q_path)
    ids_auto, q_num, dropped = coerce_numeric_frame(raw)

    if dropped:
        sys.stderr.write(f"[warn] Dropped {len(dropped)} empty column(s) from Q (likely trailing delimiters): {dropped}\n")

    if q_has_ids:
        # If user said IDs are in first column but auto-detect disagrees, try to recover
        if ids_auto is None:
            # Force treat first column as IDs (even if numeric-looking)
            ids = raw.iloc[:,0].astype(str).tolist()
            q_only = raw.iloc[:,1:].apply(pd.to_numeric, errors="coerce")
            # Drop all-NaN cols and normalize
            keep = ~q_only.isna().all(axis=0)
            if (~keep).any():
                dropped2 = [int(i) for i, k in zip(q_only.columns, keep) if not k]
                sys.stderr.write(f"[warn] Dropped {len(dropped2)} empty column(s) after forcing ID column: {dropped2}\n")
            q_only = q_only.loc[:, keep].fillna(0.0)
            arr = q_only.to_numpy(float)
            rs = arr.sum(axis=1)
            arr = np.where(rs[:,None] > 0, arr/rs[:,None], arr)
            q = pd.DataFrame(arr, columns=range(arr.shape[1]))
        else:
            ids = ids_auto.astype(str).tolist()
            q = q_num
    else:
        if ids_path:
            ids = [ln.strip().split()[0] for ln in open(ids_path) if ln.strip()]
            if len(ids) != q_num.shape[0]:
                raise SystemExit(f"[error] --ids count ({len(ids)}) != Q rows ({q_num.shape[0]}).")
            q = q_num
        else:
            # If IDs not provided and we auto-detected IDs, use them. Otherwise synthesize IDs.
            if ids_auto is not None:
                sys.stderr.write("[info] Auto-detected ID column in Q; using it.\n")
                ids = ids_auto.astype(str).tolist()
                q = q_num
            else:
                ids = [f"ind{i}" for i in range(q_num.shape[0])]
                q = q_num

    return ids, q

def main():
    ap = argparse.ArgumentParser(description="Assign clusters to individuals from a Q matrix, robust to messy delimiters.")
    ap.add_argument('-q', "--q", required=True, help="Q matrix file (ADMIXTURE/NGSadmix/fastNGSadmix style).")
    ap.add_argument('-I', "--q-has-ids", action="store_true", help="Treat first column of --q as sample IDs.")
    ap.add_argument('-i', "--ids", help="Text file with sample IDs (one per line) matching Q row order (used if --q-has-ids not given or auto-detect fails).")
    ap.add_argument('-c', "--cluster-names", help="Comma-separated names for clusters (length K). Default: C1..CK")
    ap.add_argument('-n', "--min-prop", type=float, default=0.0, help="Minimum top proportion to assign (default 0).")
    ap.add_argument('-m', "--min-delta", type=float, default=0.0, help="Minimum margin (top - second) to assign (default 0).")
    ap.add_argument('-a', "--ambig-label", default="Admixed", help="Label for ambiguous assignments.")
    ap.add_argument('-s', "--slice-k", type=int, default=None,
                    help="Force K by slicing to first K columns (useful if your Q has an accidental extra column).")
    ap.add_argument('-o', "--out", required=True, help="Output TSV with assignments and proportions.")
    ap.add_argument('-u', "--summary-out", help="Optional TSV with counts per assigned label.")
    args = ap.parse_args()

    sample_ids, q = read_inputs(args.q, args.q_has_ids, args.ids)

    # Optionally force K (e.g., if you KNOW it should be 3 but file has a junk 4th col)
    if args.slice_k is not None:
        if args.slice_k < 1 or args.slice_k > q.shape[1]:
            raise SystemExit(f"[error] --slice-k must be between 1 and current K={q.shape[1]}.")
        if args.slice_k != q.shape[1]:
            sys.stderr.write(f"[warn] Slicing Q from K={q.shape[1]} to K={args.slice_k}\n")
            q = q.iloc[:, :args.slice_k]

    K = q.shape[1]
    if args.cluster_names:
        names = [x.strip() for x in args.cluster_names.split(",") if x.strip() != ""]
        if len(names) != K:
            raise SystemExit(f"[error] --cluster-names length ({len(names)}) != K ({K}). "
                             f"Use --slice-k to match, or fix the names.")
    else:
        names = [f"C{i+1}" for i in range(K)]

    arr = q.to_numpy(float)
    # Safety: re-normalize tiny drift
    rs = arr.sum(axis=1)
    arr = np.where(rs[:,None] > 0, arr/rs[:,None], arr)

    top_idx = arr.argmax(axis=1)
    if K > 1:
        sorted_idx = np.argsort(-arr, axis=1)
        second_val = arr[np.arange(arr.shape[0]), sorted_idx[:, 1]]
    else:
        second_val = np.zeros(arr.shape[0])
    top_val = arr[np.arange(arr.shape[0]), top_idx]
    delta = top_val - second_val

    assignable = (top_val >= args.min_prop) & (delta >= args.min_delta)
    labels = np.where(assignable, np.array(names)[top_idx], args.ambig_label)

    out = pd.DataFrame({
        "sample_id": sample_ids,
        "assigned_cluster": labels,
        "max_prop": np.round(top_val, 6),
        "delta_next": np.round(delta, 6),
    })
    for i, cname in enumerate(names):
        out[f"q_{cname}"] = np.round(arr[:, i], 6)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    if args.summary_out:
        summary = (
            out["assigned_cluster"]
            .value_counts(dropna=False)
            .rename_axis("assigned_cluster")
            .reset_index(name="count")
        )
        Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary_out, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    sys.stderr.write(f"[info] Parsed K={K}; wrote {out.shape[0]} rows to {args.out}\n")
    sys.stderr.write(f"[info] Assignable={int(assignable.sum())} / {out.shape[0]} "
                     f"(min_prop={args.min_prop}, min_delta={args.min_delta})\n")

if __name__ == "__main__":
    main()