#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# KING-only relationship classification with subtype inference.
# Interface:
#   -i/--input  : required input file (TSV/CSV/whitespace)
#   -o/--out    : output TSV (default relationship_calls.tsv)
#   -n/--named  : positional mode where IDs are col3 & col4 (1-based), KING is col5 (or col7 fallback)
#
# Non-named mode:
#   - IDs are ALWAYS taken from column 1 and 2 (1-based).
#   - KING is detected by header (king/phi/phi_hat/kinship/...) or the first numeric column after column 2.
#
# Output columns only (in this order):
#   ID1  ID2  r_ab  theta  king  k0  k1  k2  ibs0  final_degree  degree_basis  relationship_inferred

import argparse, sys, re
from typing import Optional, Tuple, List
import pandas as pd

DESIRED_COLUMNS = [
    "ID1","ID2","r_ab","theta","king","k0","k1","k2","ibs0",
    "final_degree","degree_basis","relationship_inferred"
]

def _normalize(s: str) -> str:
    return re.sub(r"[-._\s]+", "", s.lower())

def _guess_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    norm_cols = {_normalize(c): c for c in df.columns}
    for cand in candidates:
        nc = _normalize(cand)
        if nc in norm_cols:
            return norm_cols[nc]
    for c in df.columns:
        nc = _normalize(c)
        if any(_normalize(cand) in nc for cand in candidates):
            return c
    return None

def _maybe_remap_ngsrelate_k012(df: pd.DataFrame) -> pd.DataFrame:
    # Map ngsRelate J7/J8/J9 -> k2/k1/k0 if explicit k0-2 not present.
    cols = { _normalize(c): c for c in df.columns }
    j7 = cols.get("j7"); j8 = cols.get("j8"); j9 = cols.get("j9")
    has_k0 = "k0" in cols or "kingk0" in cols
    has_k1 = "k1" in cols or "kingk1" in cols
    has_k2 = "k2" in cols or "kingk2" in cols
    out = df.copy()
    # User mapping: K0==J9; K1==J8; K2==J7
    if j9 is not None and not has_k0:
        out["k0"] = out[j9]
    if j8 is not None and not has_k1:
        out["k1"] = out[j8]
    if j7 is not None and not has_k2:
        out["k2"] = out[j7]
    return out

def _classify_by_king(phi: float) -> Tuple[str,str]:
    if pd.isna(phi): return ("Unknown","king")
    if phi > 0.354: return ("Monozygotic/Duplicate","king")
    if phi > 0.1875: return ("1st degree","king")
    if phi > 0.09375: return ("2nd degree","king")
    if phi > 0.046875: return ("3rd degree","king")
    if phi > 0.0234375: return ("4th degree","king")
    if phi > 0.01171875: return ("5th degree","king")
    return ("Unrelated","king")

def _infer_subtype(phi,k0,k1,k2,ibs0,degree_label) -> str:
    # Heuristic subtype calls using k0/k1/k2/ibs0 (if present).
    if degree_label == "Monozygotic/Duplicate": return "Duplicates/Clonal"
    if pd.isna(phi): return "Unspecified"
    k0v = float(k0) if k0 is not None else float('nan')
    k1v = float(k1) if k1 is not None else float('nan')
    k2v = float(k2) if k2 is not None else float('nan')
    r0v = float(ibs0) if ibs0 is not None else float('nan')
    if degree_label == "1st degree":
        po_like = ((pd.isna(k2v) or k2v<=0.05) and (pd.isna(k1v) or k1v>=0.85)
                   and (pd.isna(k0v) or k0v<=0.15) and (pd.isna(r0v) or r0v<=0.01))
        if po_like: return "PO"
        fs_like = ((not pd.isna(k2v) and k2v>=0.05) or (not pd.isna(r0v) and r0v>0.01))
        if fs_like: return "FS"
        return "1st-degree (unspecified)"
    if degree_label == "2nd degree":
        hs_like = ((pd.isna(k2v) or k2v<=0.03) and (pd.isna(k1v) or 0.3<=k1v<=0.7) and (pd.isna(k0v) or 0.3<=k0v<=0.7))
        gp_av_like = ((not pd.isna(k0v) and k0v>=0.65) and (pd.isna(k2v) or k2v<=0.03))
        if gp_av_like and not hs_like: return "GP/Avuncular"
        if hs_like: return "HS"
        return "2nd-degree (unspecified)"
    if degree_label == "3rd degree":
        cousin_like = ((pd.isna(k2v) or k2v<=0.02) and (pd.isna(k1v) or 0.1<=k1v<=0.35) and (pd.isna(k0v) or k0v>=0.6))
        if cousin_like: return "Cousins"
        return "3rd-degree (unspecified)"
    if degree_label in {"4th degree","5th degree"}:
        return degree_label+" (unspecified)"
    return "Unspecified"

def _read_table(path:str)->pd.DataFrame:
    for sep in ["\t",",",";",r"\s+"]:
        try:
            df=pd.read_csv(path,sep=sep,engine="python")
            if df.shape[1]>1: return df
        except Exception:
            continue
    return pd.read_csv(path,sep=None,engine="python")

def _select_named_columns(df: pd.DataFrame):
    # Named mode: IDs are in columns 3 and 4 (1-based), KING in column 5 (or 7 fallback).
    if df.shape[1] < 5:
        sys.exit("[ERROR] Named mode expects at least 5 columns (IDs in 3&4, KING in 5).")
    id1 = df.iloc[:, 2].astype(str)
    id2 = df.iloc[:, 3].astype(str)
    king5 = pd.to_numeric(df.iloc[:, 4], errors="coerce")
    if king5.notna().any():
        king = king5
    elif df.shape[1] >= 7:
        king7 = pd.to_numeric(df.iloc[:, 6], errors="coerce")
        if king7.notna().any():
            king = king7
        else:
            sys.exit("[ERROR] Named mode: KING not found in column 5 or 7 (numeric parse failed).")
    else:
        sys.exit("[ERROR] Named mode: KING not found in column 5, and no column 7 present.")
    return id1, id2, king

def main():
    ap=argparse.ArgumentParser(description="KING-only classification + subtype inference. Outputs ONLY core columns.")
    ap.add_argument("-i","--input", required=True, help="Input TSV/CSV/whitespace-delimited file")
    ap.add_argument("-o","--out", default="relationship_calls.tsv", help="Output TSV")
    ap.add_argument("-n","--named", action="store_true",
                    help="Positional mode: IDs=col3&col4; KING=col5 (fallback col7)")
    args=ap.parse_args()

    df=_read_table(args.input)
    if df.empty: sys.exit("[ERROR] Empty input")

    df=_maybe_remap_ngsrelate_k012(df)

    # IDs and KING selection
    if args.named:
        ID1, ID2, KING = _select_named_columns(df)
    else:
        if df.shape[1] < 2:
            sys.exit("[ERROR] Need at least two columns for IDs in non-named mode.")
        ID1 = df.iloc[:, 0].astype(str)
        ID2 = df.iloc[:, 1].astype(str)
        king_col = _guess_col(df, ["king","phi","phi_hat","kinship","king_phi","king_kinship"])
        if king_col is None:
            king_idx = None
            for i, c in enumerate(df.columns):
                if i in (0,1):
                    continue
                try:
                    pd.to_numeric(df[c])
                    king_idx = i
                    break
                except Exception:
                    continue
            if king_idx is None:
                sys.exit("[ERROR] Could not determine KING column automatically. Try -n named mode if positions are known.")
            KING = pd.to_numeric(df.iloc[:, king_idx], errors="coerce")
        else:
            KING = pd.to_numeric(df[king_col], errors="coerce")

    # Optional metrics by header autodetect (not positional)
    r_ab_col=_guess_col(df,["r_ab","rab","rel_ab"])
    theta_col=_guess_col(df,["theta","kinship_theta","ngsrelate_theta","king_theta"])
    k0_col=_guess_col(df,["k0","kingk0"])
    k1_col=_guess_col(df,["k1","kingk1"])
    k2_col=_guess_col(df,["k2","kingk2"])
    ibs0_col=_guess_col(df,["ibs0","r0","ibs_0"])

    out_cols={
      "ID1": ID1,
      "ID2": ID2,
      "r_ab": pd.to_numeric(df[r_ab_col],errors="coerce") if r_ab_col else pd.Series([float('nan')]*len(df)),
      "theta": pd.to_numeric(df[theta_col],errors="coerce") if theta_col else pd.Series([float('nan')]*len(df)),
      "king": KING,
      "k0": pd.to_numeric(df[k0_col],errors="coerce") if k0_col else pd.Series([float('nan')]*len(df)),
      "k1": pd.to_numeric(df[k1_col],errors="coerce") if k1_col else pd.Series([float('nan')]*len(df)),
      "k2": pd.to_numeric(df[k2_col],errors="coerce") if k2_col else pd.Series([float('nan')]*len(df)),
      "ibs0": pd.to_numeric(df[ibs0_col],errors="coerce") if ibs0_col else pd.Series([float('nan')]*len(df)),
    }
    out=pd.DataFrame(out_cols)

    # KING-only degree
    degs,bases=[],[]
    for phi in out["king"]:
        d,b=_classify_by_king(phi); degs.append(d); bases.append(b)
    out["final_degree"]=degs; out["degree_basis"]=bases

    # Subtype inference
    out["relationship_inferred"]=[_infer_subtype(phi,k0,k1,k2,r0,d)
                                  for phi,k0,k1,k2,r0,d in zip(out["king"],out["k0"],out["k1"],out["k2"],out["ibs0"],out["final_degree"])]

    # Only desired columns
    out = out.reindex(columns=DESIRED_COLUMNS)
    out.to_csv(args.out, sep="\t", index=False, float_format="%.6f", na_rep="")
    sys.stderr.write(f"[OK] Wrote KING-only classifications with subtypes to {args.out}\n")

if __name__=="__main__": main()
