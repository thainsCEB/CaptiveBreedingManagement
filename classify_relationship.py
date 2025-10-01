#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Classify relationships by combining ngsRelate r_ab, theta (kinship), and KING kinship.
Optional:
  - Generate GRM from r_ab (-g/--grm)
  - Infer specific relationship subtype (PO/FS/HS/Grandparent-Avuncular/Cousin-like) (-R/--infer-subtype)
    using k0,k1,k2 and R0 (IBS0).

Lean output (tab-delimited):
  ID1 ID2 r_ab theta king k0 k1 k2 ibs0 degree_rab degree_theta degree_king final_degree degree_basis [relationship_inferred]
"""

import argparse, sys, math, csv
from typing import Tuple, List, Dict, Optional
import pandas as pd
from collections import defaultdict

# ---- Degree labels & bins ----
DEGREE_LABELS = {1:"1st degree",2:"2nd degree",3:"3rd degree",4:"4th degree",5:"5th degree",6:"Unrelated"}

RAB_BINS = [
    (1, 0.375,  1.01),
    (2, 0.1875, 0.375),
    (3, 0.09375,0.1875),
    (4, 0.046875,0.09375),
    (5, 0.0234375,0.046875)
]
RAB_MIN_RELATED = 0.015

THETA_HALFR_BINS = [
    (1, 0.1875,  1.01),
    (2, 0.09375, 0.1875),
    (3, 0.046875,0.09375),
    (4, 0.0234375,0.046875),
    (5, 0.01171875,0.0234375),
]
THETA_HALFR_MIN_RELATED = 0.0075

KING_STYLE_BINS = [
    (1, 0.177,   1.01),
    (2, 0.0884,  0.177),
    (3, 0.0442,  0.0884),
    (4, 0.0221,  0.0442),
    (5, 0.01105, 0.0221),
]
KING_STYLE_MIN_RELATED = 0.0075

RAB_CENTERS = [(1,0.5),(2,0.25),(3,0.125),(4,0.0625),(5,0.03125)]
KIN_CENTERS = [(d, r/2.0) for d,r in RAB_CENTERS]

# ---- Helpers ----
def val_or_none(x):
    try:
        if pd.isna(x): return None
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def classify_by_bins(value: float, bins: List[tuple], min_related: float) -> int:
    if value is None or value < min_related:
        return 6
    for deg, lo, hi in bins:
        if lo <= value < hi:
            return deg
    return 6

def continuous_degree(value: float, centers: List[tuple]) -> float:
    centers_sorted = sorted(centers, key=lambda t: t[1], reverse=True)
    if value is None: return 6.0
    if value >= centers_sorted[0][1]:
        return float(centers_sorted[0][0])
    if value <= centers_sorted[-1][1]:
        low_deg, low_val = centers_sorted[-1]
        return min(6.0, low_deg + (low_val - value)/max(low_val,1e-9)*(6-low_deg))
    for i in range(len(centers_sorted)-1):
        d1,v1 = centers_sorted[i]; d2,v2 = centers_sorted[i+1]
        if v1 >= value >= v2:
            t = (v1 - value)/max((v1 - v2),1e-12)
            return d1*(1-t) + d2*t
    return 6.0

def vote_and_combine(deg_calls: List[int], cont_degs: List[float],
                     strict: bool=False, lenient: bool=False) -> str:
    calls = [d for d in deg_calls if d != 6]
    if not calls:
        meanc = sum(cont_degs)/len(cont_degs) if cont_degs else 6.0
        rd = int(round(meanc))
        return DEGREE_LABELS[rd] if 1 <= rd <= 6 else DEGREE_LABELS[6]
    from collections import Counter
    c = Counter(calls)
    top_deg, top_n = c.most_common(1)[0]
    ties = [k for k,v in c.items() if v==top_n]
    if len(ties) == 1:
        span = max(calls) - min(calls)
        ambiguous = (span > 0 if strict else span > (2 if lenient else 1))
        return "Ambiguous" if ambiguous else DEGREE_LABELS[top_deg]
    meanc = sum(cont_degs)/len(cont_degs) if cont_degs else float(ties[0])
    tie_deg = min(max(int(round(meanc)),1),6)
    span = max(calls) - min(calls)
    ambiguous = (True if strict else (span > (2 if lenient else 1)))
    return "Ambiguous" if ambiguous else DEGREE_LABELS[tie_deg]

# ---- Relationship inference (uses k0,k1,k2,R0) ----
PO_K2_MAX = 0.05
PO_R0_MAX = 0.005
FS_K2_BAND = (0.10, 0.40)
HS_K0_BAND = (0.40, 0.62)
AVGP_K0_BAND = (0.62, 0.90)

def infer_subtype(final_degree: str,
                  k0: Optional[float], k1: Optional[float],
                  k2: Optional[float], r0: Optional[float]) -> str:
    if final_degree == "1st degree":
        if (k2 is not None and k2 <= PO_K2_MAX) and (r0 is not None and r0 <= PO_R0_MAX):
            return "PO"
        if k2 is not None and FS_K2_BAND[0] <= k2 <= FS_K2_BAND[1]:
            return "FS"
        return "1st-degree (PO/FS)"
    if final_degree == "2nd degree":
        if k0 is not None and HS_K0_BAND[0] <= k0 <= HS_K0_BAND[1]:
            return "HS"
        if k0 is not None and AVGP_K0_BAND[0] <= k0 <= AVGP_K0_BAND[1]:
            return "Grandparent/Avuncular"
        return "2nd-degree (HS/GP/Avuncular)"
    if final_degree == "3rd degree":
        return "Cousin-like"
    return ""

# ---- Main ----
def main():
    ap = argparse.ArgumentParser(description="Classify relationships; optionally infer subtype and emit GRM.")
    ap.add_argument("-i","--input", required=True, help="Input TSV/CSV (delimiter auto-detected).")
    ap.add_argument("-o","--out", default=None, help="Write table output (default: stdout).")
    ap.add_argument("-T","--theta-mode", choices=["half-r","king-style"], default="half-r",
                    help="Theta thresholds (default half-r).")
    ap.add_argument("-K","--king-mode", choices=["king-style","half-r"], default="king-style",
                    help="KING thresholds (default king-style).")
    ap.add_argument("-p","--primary-metric", choices=["consensus","r_ab","theta","king"], default="consensus",
                    help="Metric driving final_degree.")
    ap.add_argument("-s","--strict", action="store_true", help="Consensus: any disagreement → Ambiguous.")
    ap.add_argument("-l","--lenient", action="store_true", help="Consensus: only span>2 → Ambiguous.")
    ap.add_argument("-m","--min-metrics", type=int, default=1, help="Consensus: require ≥N metrics.")
    ap.add_argument("-n","--named", action="store_true", help="Use 3rd/4th columns as ID1/ID2.")
    ap.add_argument("-g","--grm", default=None, help="Write GRM from r_ab (square matrix with IDs as header & index).")
    ap.add_argument("-R","--infer-subtype", action="store_true",
                    help="Add relationship_inferred (PO/FS/HS/Grandparent-Avuncular/Cousin-like).")
    args = ap.parse_args()

    # Detect delimiter
    with open(args.input, "rb") as fh:
        sample = fh.read(1024)
    try:
        dialect = csv.Sniffer().sniff(sample.decode("utf-8", errors="ignore"))
        sep = dialect.delimiter
    except Exception:
        text = sample.decode("utf-8", errors="ignore")
        sep = "\t" if "\t" in text else ("," if "," in text else r"\s+")

    df = pd.read_csv(args.input, sep=sep, engine="python", dtype=str)

    # IDs
    if args.named:
        if df.shape[1] < 4:
            sys.stderr.write("ERROR: --named requires ≥4 columns.\n"); sys.exit(2)
        c_id1, c_id2 = df.columns[2], df.columns[3]
    else:
        lc = {c.lower(): c for c in df.columns}
        def col(opts):
            for nm in opts:
                if nm in lc: return lc[nm]
            return None
        c_id1 = col(["id1","ind1","sample1","i1"]) or (df.columns[0] if df.shape[1]>=2 else None)
        c_id2 = col(["id2","ind2","sample2","i2"]) or (df.columns[1] if df.shape[1]>=2 else None)
        if c_id1 is None or c_id2 is None:
            sys.stderr.write("ERROR: cannot determine ID columns.\n"); sys.exit(2)

    # Metrics (case-insensitive)
    cols = {c.lower(): c for c in df.columns}
    def col(name_opts):
        for nm in name_opts:
            if nm in cols: return cols[nm]
        return None
    c_rab   = col(["r_ab","rab","relatedness","r"])
    c_theta = col(["theta","kinship","phi"])
    c_king  = col(["king","king_kinship","king_phi"])
    # ngsRelate aliases: K0==J9, K1==J8, K2==J7
    c_k0    = col(["k0","ibd0","p_ibd0","j9"])
    c_k1    = col(["k1","ibd1","p_ibd1","j8"])
    c_k2    = col(["k2","ibd2","p_ibd2","j7"])
    # IBS0 a.k.a. R0/RO
    c_ibs0  = col(["r0","ro","ibs0","p_ibs0","ibs_0","king_r0","king_ro"])

    # Threshold modes
    theta_bins = THETA_HALFR_BINS if args.theta_mode=="half-r" else KING_STYLE_BINS
    theta_min  = THETA_HALFR_MIN_RELATED if args.theta_mode=="half-r" else KING_STYLE_MIN_RELATED
    king_bins  = KING_STYLE_BINS if args.king_mode=="king-style" else THETA_HALFR_BINS
    king_min   = KING_STYLE_MIN_RELATED if args.king_mode=="king-style" else THETA_HALFR_MIN_RELATED

    out_rows = []
    r_ab_pairs: Dict[Tuple[str,str], List[float]] = defaultdict(list)
    id_set = set()

    for _, row in df.iterrows():
        id1 = str(row[c_id1]); id2 = str(row[c_id2])
        id_set.update([id1, id2])

        rab  = val_or_none(row[c_rab])   if c_rab   else None
        th   = val_or_none(row[c_theta]) if c_theta else None
        king = val_or_none(row[c_king])  if c_king  else None
        k0   = val_or_none(row[c_k0])    if c_k0    else None
        k1   = val_or_none(row[c_k1])    if c_k1    else None
        k2   = val_or_none(row[c_k2])    if c_k2    else None
        r0   = val_or_none(row[c_ibs0])  if c_ibs0  else None

        # For GRM (use r_ab only)
        if rab is not None and id1 != id2:
            a,b = (id1,id2) if id1 <= id2 else (id2,id1)
            r_ab_pairs[(a,b)].append(rab)

        # Per-metric degree calls
        deg_rab   = classify_by_bins(rab,  RAB_BINS,  RAB_MIN_RELATED) if rab  is not None else ""
        deg_theta = classify_by_bins(th,   theta_bins, theta_min)       if th   is not None else ""
        deg_king  = classify_by_bins(king, king_bins,  king_min)        if king is not None else ""

        # Final degree
        degree_basis = args.primary_metric
        if args.primary_metric == "r_ab":
            final = "Ambiguous" if rab  is None else DEGREE_LABELS[deg_rab]
        elif args.primary_metric == "theta":
            final = "Ambiguous" if th   is None else DEGREE_LABELS[deg_theta]
        elif args.primary_metric == "king":
            final = "Ambiguous" if king is None else DEGREE_LABELS[deg_king]
        else:
            present = sum([rab is not None, th is not None, king is not None])
            if present < args.min_metrics:
                final = "Ambiguous"
            else:
                disc_calls = [d for d in [
                    deg_rab if deg_rab != "" else 6,
                    deg_theta if deg_theta != "" else 6,
                    deg_king if deg_king != "" else 6] if d != 6]
                cont_vals  = [
                    continuous_degree(rab,  RAB_CENTERS)  if rab  is not None else 6.0,
                    continuous_degree(th,   KIN_CENTERS)   if th   is not None else 6.0,
                    continuous_degree(king, KIN_CENTERS)   if king is not None else 6.0]
                final = vote_and_combine(disc_calls, cont_vals, strict=args.strict, lenient=args.lenient)
            degree_basis = "consensus"

        row_out = {
            "ID1": id1, "ID2": id2,
            "r_ab": "" if rab  is None else f"{rab:.6f}",
            "theta": "" if th  is None else f"{th:.6f}",
            "king": "" if king is None else f"{king:.6f}",
            "k0": "" if k0 is None else f"{k0:.6f}",
            "k1": "" if k1 is None else f"{k1:.6f}",
            "k2": "" if k2 is None else f"{k2:.6f}",
            "ibs0": "" if r0 is None else f"{r0:.6f}",
            "degree_rab":   "" if deg_rab   == "" else DEGREE_LABELS[deg_rab],
            "degree_theta": "" if deg_theta == "" else DEGREE_LABELS[deg_theta],
            "degree_king":  "" if deg_king  == "" else DEGREE_LABELS[deg_king],
            "final_degree": final,
            "degree_basis": degree_basis
        }

        if args.infer_subtype:
            row_out["relationship_inferred"] = infer_subtype(final, k0, k1, k2, r0)

        out_rows.append(row_out)

    # Write table
    cols_out = [
        "ID1","ID2","r_ab","theta","king","k0","k1","k2","ibs0",
        "degree_rab","degree_theta","degree_king",
        "final_degree","degree_basis"
    ]
    if args.infer_subtype:
        cols_out.append("relationship_inferred")

    out_df = pd.DataFrame(out_rows, columns=cols_out)
    if args.out:
        out_df.to_csv(args.out, sep="\t", index=False)
    else:
        out_df.to_csv(sys.stdout, sep="\t", index=False)

    # ---- Write GRM (square matrix with same ordered IDs as header and index) ----
    if args.grm:
        ids = sorted(id_set)
        grm = pd.DataFrame(float("nan"), index=ids, columns=ids, dtype=float)  # square, same order
        # diagonal = 1.0
        for i in ids:
            grm.loc[i, i] = 1.0
        # off-diagonals from averaged r_ab
        for (a, b), vals in r_ab_pairs.items():
            if vals:
                v = sum(vals) / len(vals)
                grm.loc[a, b] = v
                grm.loc[b, a] = v
        grm.index.name = None  # blank top-left cell
        grm.to_csv(args.grm, sep="\t", index=True, header=True)

if __name__ == "__main__":
    main()
