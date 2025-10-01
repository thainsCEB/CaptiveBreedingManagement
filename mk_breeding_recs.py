#!/usr/bin/env python3
"""
Breeding recommendations from genomics using:
  - Mean Kinship (MK) minimization
  - Founder representation balancing (optional)
  - Inbreeding avoidance (F_offspring = φ_ij)
  - Taxon/population purity constraints (same-group default)
  - Pair labels: BREED / HOLD / DO NOT BREED
  - Individual-level "do not pair" reasons

INPUT (auto-detected if you use --input/-i):
  1) Square kinship/GRM matrix (IDs as both header and index)
  2) NGSRelate long pairwise table (columns ID1/ID2/... + metric)

You can also force a specific input by using one of:
  --matrix/-m, --grm/-g, or --ngsrelate/-n

Common transforms:
  - GRM with relationship r ≈ 2φ  -> use --transform half
  - NGSRelate r_ab (≈ 2φ)         -> use --ngsrelate-transform half
  - NGSRelate theta (≈ φ)         -> default as_kinship (no scaling)

Author: ChatGPT
"""
import argparse, sys
import pandas as pd, numpy as np
from pathlib import Path
from typing import Optional, Dict, Tuple

# -------------------- IO helpers --------------------

def _read_table(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path, sep="\t")

def _coerce_sex(x: str) -> Optional[str]:
    if pd.isna(x): return None
    s = str(x).strip().lower()
    if s in {"m","male"}: return "M"
    if s in {"f","female"}: return "F"
    return None

def _normalize_labels(idx_or_cols):
    # Cast to string, strip, keep empty as-is
    return [str(x).strip() for x in idx_or_cols]

def _set_index_if_id_col(df: pd.DataFrame) -> pd.DataFrame:
    """If first column looks like IDs (mostly non-numeric), use it as index."""
    first = df.iloc[:, 0]
    if not np.issubdtype(first.dtype, np.number):
        return df.set_index(df.columns[0])
    # If first col numeric but header row includes same labels as columns (e.g., 0..n-1),
    # we still leave as-is; later normalization will try to align names as strings.
    return df

def _try_make_square(df: pd.DataFrame, label_hint: str) -> pd.DataFrame:
    """
    Make a square symmetric numeric matrix with matching index/columns.
    Tries several fixes that commonly cause "Matrix not square" errors.
    """
    # 1) If first col is ID-like, use as index
    if df.shape[1] >= 2:
        df = _set_index_if_id_col(df)

    # 2) Normalize labels to strings
    df.index = _normalize_labels(df.index)
    df.columns = _normalize_labels(df.columns)

    # 3) If columns missing, but first row looks like IDs (e.g., file had an extra header row),
    # try using the first row as header and re-parse.
    missing = set(df.index) - set(df.columns)
    if missing:
        # Heuristic: if many column names look like integers 0..N-1, but row labels are real IDs,
        # we might have read a no-header file. Try re-reading with header=0 and set index.
        # We already used pandas auto-sep, so instead try: treat first row as header by hand.
        # If df has an unnamed first column 'Unnamed: 0' and that equals index, handle that too.
        # Last resort: align by intersection only (warn).
        pass

    # 4) Align to intersection; if that drops too much, we’ll warn
    inter = df.columns.intersection(df.index)
    if len(inter) == 0:
        # Try one more thing: if the index equals the first column values but columns are 0..N-1,
        # assume the FIRST row contained column IDs. Promote that row to header.
        try:
            first_row = df.iloc[0, :].values
            # New columns from first row, then drop it and set index from current index column
            df2 = df.iloc[1:, :].copy()
            df2.columns = _normalize_labels(first_row)
            # If left-most column duplicates the index names, set it as index
            if not np.issubdtype(df2.iloc[:,0].dtype, np.number):
                df2 = df2.set_index(df2.columns[0])
            df2.index = _normalize_labels(df2.index)
            df2.columns = _normalize_labels(df2.columns)
            inter = df2.columns.intersection(df2.index)
            if len(inter) > 0:
                df = df2
        except Exception:
            pass

    inter = df.columns.intersection(df.index)
    if len(inter) == 0:
        sys.exit(f"[ERROR] Could not reconcile {label_hint} labels between rows and columns. "
                 f"Ensure the first column is IDs and header row has the same IDs.")

    # Reindex square on intersection and ensure numeric
    M = df.loc[inter, inter].copy()
    try:
        M = M.astype(float)
    except Exception:
        # Try to coerce, converting non-numeric to NaN then fill later
        M = M.apply(pd.to_numeric, errors="coerce")

    # Fill any diagonal NAs with 0.5; off-diagonal NAs with median off-diagonal
    vals = M.values
    if np.isnan(vals).any():
        np.fill_diagonal(vals, np.where(np.isfinite(np.diag(vals)), np.diag(vals), 0.5))
        off = vals[np.triu_indices_from(vals, k=1)]
        med = np.nanmedian(off) if np.isfinite(off).any() else 0.0
        vals[np.isnan(vals)] = med
        M = pd.DataFrame(vals, index=M.index, columns=M.columns)

    # Symmetrize (average with transpose)
    M = (M + M.T) / 2.0
    return M

def _read_square_matrix(path: str, label_hint: str = "matrix/GRM") -> pd.DataFrame:
    raw = _read_table(path)
    return _try_make_square(raw, label_hint)

# -------------------- NGSRelate handling --------------------

def _metric_to_phi(series: pd.Series, metric: str, transform: str) -> pd.Series:
    x = series.astype(float)
    if metric.lower() == "theta":
        phi = x; default = "as_kinship"
    elif metric.lower() in {"r_ab","rab","r"}:
        phi = x/2.0; default = "half"
    elif metric.upper() in {"J7","J8","J9","R0"}:
        phi = x; default = "none"
    else:
        phi = x; default = "none"
    mode = (transform or default).lower()
    if mode == "as_kinship": pass
    elif mode == "half": phi *= 0.5
    elif mode == "double": phi *= 2.0
    elif mode == "none": pass
    else: sys.exit(f"[ERROR] Unknown --ngsrelate-transform: {transform}")
    return phi

def _build_phi_from_ngsrelate(path: str, metric: str, transform: str) -> pd.DataFrame:
    df = _read_table(path)
    # Detect ID columns
    id1_cols = [c for c in df.columns if str(c).strip().lower() in {"id1","ind1","sample1","i","a"}]
    id2_cols = [c for c in df.columns if str(c).strip().lower() in {"id2","ind2","sample2","j","b"}]
    if not id1_cols or not id2_cols:
        id1, id2 = df.columns[:2]
    else:
        id1, id2 = id1_cols[0], id2_cols[0]
    # Find metric column
    lc = {c.lower(): c for c in df.columns}
    metric_aliases = {
        "theta": ["theta","coancestry","kinship","phi","phi_ij"],
        "r_ab": ["r_ab","rab","relatedness","r"],
        "j7": ["j7","k2"],
        "j8": ["j8","k1"],
        "j9": ["j9","k0"],
        "r0": ["r0","ibs0"],
    }
    wanted = metric.lower()
    col_metric = None
    for cand in metric_aliases.get(wanted, [wanted]):
        if cand in lc: col_metric = lc[cand]; break
    if col_metric is None:
        if metric in df.columns: col_metric = metric
        else: sys.exit(f"[ERROR] Metric '{metric}' not found. Columns: {list(df.columns)}")

    A = df[[id1, id2, col_metric]].copy()
    A.columns = ["ID1","ID2","VAL"]
    A["VAL"] = _metric_to_phi(A["VAL"], metric, transform)

    ids = pd.Index(pd.unique(pd.concat([A["ID1"], A["ID2"]], ignore_index=True)), name="ID")
    K = pd.DataFrame(np.nan, index=ids, columns=ids, dtype=float)
    for r in A.itertuples(index=False):
        K.at[r.ID1, r.ID2] = r.VAL
        K.at[r.ID2, r.ID1] = r.VAL

    # Set diagonals to 0.5*(1+F_i); if F unknown, use 0.5
    np.fill_diagonal(K.values, 0.5)
    # Fill missing with median off-diagonal
    off = K.values[np.triu_indices_from(K.values, k=1)]
    med = np.nanmedian(off) if np.isfinite(off).any() else 0.0
    K = K.fillna(med)
    return K

def _apply_transform_matrix(K: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode == "none": return K
    if mode == "half": return K*0.5
    if mode == "double": return K*2.0
    if mode == "auto":
        d = np.diag(K.values)
        return K if 0.35 <= np.median(d) <= 0.65 else K
    sys.exit(f"[ERROR] Unknown --transform: {mode}")

# -------------------- Scoring helpers --------------------

def compute_mean_kinship(K: pd.DataFrame) -> pd.Series:
    return K.mean(axis=1)

def founder_penalties(ids, founder_map: Dict[str,str], weight: float, target_fr: Optional[float]=None):
    if not founder_map: return {}, {}
    groups = pd.Series({i: founder_map.get(i, None) for i in ids}, dtype="object").dropna()
    if groups.empty: return {}, {}
    freqs = groups.value_counts()/len(groups)
    if target_fr is None: target_fr = 1.0/len(freqs)
    over = {g: max(0.0, f - target_fr) for g, f in freqs.items()}
    max_over = max(over.values()) if over else 0.0
    if max_over > 0: over = {g: v/max_over for g, v in over.items()}
    indiv = {i: weight*over.get(groups.get(i, None), 0.0) for i in ids}
    return indiv, {g: over.get(g, 0.0) for g in freqs.index}

# -------------------- Input auto-detection --------------------

def _infer_input_kind(df: pd.DataFrame) -> str:
    """
    Decide if df looks like:
      - 'ngsrelate' long table (>=2 ID cols + metric col)
      - 'square' matrix/GRM
    Heuristics: if square-ish and headers overlap index -> 'square'. Otherwise 'ngsrelate'.
    """
    # Square-ish?
    if abs(df.shape[0] - df.shape[1]) <= 2:
        # If first col is non-numeric and appears as column header anywhere => likely square with ID column
        first = df.iloc[:, 0]
        if not np.issubdtype(first.dtype, np.number):
            return "square"
    # If it has at least two columns and looks like pair list => ngsrelate
    return "ngsrelate"

# -------------------- Main --------------------

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Inputs (auto or explicit)
    ap.add_argument("-i","--input", help="Auto-detect format from this file (matrix/GRM or NGSRelate long table).")
    ap.add_argument("-m","--matrix", help="Force square matrix format (IDs in header+index).")
    ap.add_argument("-g","--grm", help="Force square GRM/relationship matrix (IDs in header+index).")
    ap.add_argument("-n","--ngsrelate", help="Force NGSRelate long table.")
    ap.add_argument("-M","--ngsrelate-metric", default="theta",
                    choices=["theta","r_ab","J7","J8","J9","R0","rab","r0"],
                    help="NGSRelate metric to use")
    ap.add_argument("-T","--ngsrelate-transform", default="as_kinship",
                    choices=["as_kinship","half","double","none"],
                    help="Convert chosen metric to φ_ij")
    ap.add_argument("-X","--transform", choices=["none","auto","half","double"], default="none",
                    help="Transform for square matrices/GRM to φ_ij")

    # Covariates / constraints
    ap.add_argument("-I","--inb", help="Table with columns: ID, F")
    ap.add_argument("-F","--founders", help="2-col: ID, FounderGroup")
    ap.add_argument("-S","--sexes", help="2-col: ID, Sex in {M,F}")
    ap.add_argument("-A","--ages", help="2/3-col: ID, Age")
    ap.add_argument("-a","--min-age", type=float, default=None)
    ap.add_argument("-b","--max-age", type=float, default=None)
    ap.add_argument("-Z","--same-sex-ok", action="store_true")

    ap.add_argument("-Y","--taxa", help="2-col: ID, Group (taxon/population)")
    ap.add_argument("-P","--purity-mode", choices=["strict","penalize","off"], default="strict")
    ap.add_argument("-W","--w-purity", type=float, default=2.0)

    ap.add_argument("-Q","--Fmax", type=float, default=0.0625)
    ap.add_argument("-u","--w-mk", type=float, default=1.0)
    ap.add_argument("-v","--w-pair", type=float, default=2.0)
    ap.add_argument("-w","--w-fr", type=float, default=1.0)
    ap.add_argument("-r","--target-fr", type=float, default=None)

    ap.add_argument("-B","--breed-quantile", type=float, default=0.15)
    ap.add_argument("-H","--hold-quantile", type=float, default=0.60)

    ap.add_argument("-p","--prefix", default="breeding_recs")
    ap.add_argument("-L","--whitelist", help="ID list to force-include (one per line)")
    ap.add_argument("-D","--blacklist", help="ID list to exclude (one per line)")
    ap.add_argument("-d","--seed", type=int, default=1337)
    args = ap.parse_args()
    np.random.seed(args.seed)

    # --------- Decide which input we're using ---------
    src_kind = None
    if args.matrix: src_kind, src_path = "square", args.matrix
    elif args.grm:  src_kind, src_path = "square", args.grm
    elif args.ngsrelate: src_kind, src_path = "ngsrelate", args.ngsrelate
    elif args.input:
        # Peek and infer
        peek = _read_table(args.input)
        k = _infer_input_kind(peek)
        src_kind, src_path = ("square", args.input) if k == "square" else ("ngsrelate", args.input)
    else:
        sys.exit("[ERROR] Provide one of --input/-i, --matrix/-m, --grm/-g, or --ngsrelate/-n")

    # --------- Build φ matrix ---------
    if src_kind == "square":
        K = _read_square_matrix(src_path, "matrix/GRM")
        K = _apply_transform_matrix(K, args.transform)
    else:
        K = _build_phi_from_ngsrelate(src_path, args.ngsrelate_metric, args.ngsrelate_transform)

    ids = list(K.index)

    # --------- Optional covariates ---------
    def _read_id_list(p):
        s = Path(p).read_text().strip().splitlines()
        return set([x.strip().split()[0] for x in s if x.strip()])

    inb = None
    if args.inb:
        t = _read_table(args.inb)
        idc = [c for c in t.columns if str(c).lower() in {"id","sample","ind"}]
        fc  = [c for c in t.columns if str(c).lower() in {"f","inbreeding","inb"}]
        if not idc or not fc: sys.exit("[ERROR] --inb must have columns: ID, F")
        inb = t[[idc[0], fc[0]]].copy(); inb.columns = ["ID","F"]; inb = inb.set_index("ID").reindex(ids)

    sexes = None
    if args.sexes:
        t = _read_table(args.sexes)
        idc = [c for c in t.columns if str(c).lower() in {"id","sample","ind"}]
        sc  = [c for c in t.columns if str(c).lower() in {"sex","m_f","mf"}]
        if not idc or not sc: sys.exit("[ERROR] --sexes must have columns: ID, Sex")
        t = t[[idc[0], sc[0]]].copy(); t.columns = ["ID","Sex"]; t["Sex"] = t["Sex"].map(_coerce_sex)
        sexes = t.set_index("ID").reindex(ids)["Sex"]

    age_ok = pd.Series(True, index=ids)
    if args.ages:
        t = _read_table(args.ages)
        idc = [c for c in t.columns if str(c).lower() in {"id","sample","ind"}]
        ac  = [c for c in t.columns if "age" in str(c).lower()]
        if not idc or not ac: sys.exit("[ERROR] --ages must have columns: ID, Age")
        t = t[[idc[0], ac[0]]].copy(); t.columns = ["ID","Age"]
        ages = t.set_index("ID").reindex(ids)["Age"].astype(float)
        if args.min_age is not None: age_ok &= ages >= args.min_age
        if args.max_age is not None: age_ok &= ages <= args.max_age

    founder_map = {}
    if args.founders:
        t = _read_table(args.founders)
        idc = [c for c in t.columns if str(c).lower() in {"id","sample","ind"}]
        gc  = [c for c in t.columns if ("founder" in str(c).lower()) or ("group" in str(c).lower())]
        if not idc or not gc: sys.exit("[ERROR] --founders must have columns: ID, FounderGroup")
        t = t[[idc[0], gc[0]]].copy(); t.columns = ["ID","FounderGroup"]
        founder_map = {r.ID: str(r.FounderGroup) for r in t.itertuples(index=False)}

    taxa_map = {}
    if args.taxa:
        t = _read_table(args.taxa)
        idc = [c for c in t.columns if str(c).lower() in {"id","sample","ind"}]
        gc  = [c for c in t.columns if ("group" in str(c).lower()) or ("tax" in str(c).lower()) or ("pop" in str(c).lower())]
        if not idc or not gc: sys.exit("[ERROR] --taxa must have columns: ID, Group")
        t = t[[idc[0], gc[0]]].copy(); t.columns = ["ID","Group"]
        taxa_map = {r.ID: str(r.Group) for r in t.itertuples(index=False)}

    whitelist = _read_id_list(args.whitelist) if args.whitelist else None
    blacklist = _read_id_list(args.blacklist) if args.blacklist else set()

    # --------- Filter candidates ---------
    candidates = [i for i in ids if i not in blacklist and age_ok.get(i, True)]
    if not candidates: sys.exit("[ERROR] No candidates remain after filtering.")
    K = K.loc[candidates, candidates]

    # --------- MK & penalties ---------
    MK = compute_mean_kinship(K)
    indiv_fr_pen, _ = founder_penalties(candidates, founder_map, weight=args.w_fr, target_fr=args.target_fr)

    def group_of(x): return taxa_map.get(x, np.nan)
    def sex_of(x):   return None if sexes is None else sexes.get(x, None)

    # --------- Score pairs ---------
    rows = []
    for i, a in enumerate(candidates):
        for b in candidates[i+1:]:
            # sex filter
            if (not args.same_sex_ok) and (sexes is not None):
                s1, s2 = sex_of(a), sex_of(b)
                if s1 and s2 and s1 == s2:
                    continue

            phi = float(K.at[a, b])
            F_off = phi
            okF = F_off <= args.Fmax

            g1, g2 = group_of(a), group_of(b)
            same_group = (str(g1) == str(g2)) if not (pd.isna(g1) or pd.isna(g2)) else True  # unknown => allow
            purity_ok = True
            purity_penalty = 0.0
            flags = []

            if args.purity_mode == "strict":
                purity_ok = same_group
                if not same_group:
                    flags.append("cross-group")
            elif args.purity_mode == "penalize":
                purity_ok = True
                if not same_group:
                    purity_penalty = args.w_purity
                    flags.append("cross-group")

            mk_term = args.w_mk * ((MK[a] + MK[b]) / 2.0)
            pair_term = args.w_pair * F_off
            fr_term = (indiv_fr_pen.get(a, 0.0) + indiv_fr_pen.get(b, 0.0)) / 2.0
            score = mk_term + pair_term + fr_term + purity_penalty

            rows.append({
                "ID1": a, "sex1": sex_of(a), "group1": g1, "MK1": MK[a],
                "ID2": b, "sex2": sex_of(b), "group2": g2, "MK2": MK[b],
                "pair_kinship": phi, "F_offspring_pred": F_off,
                "founder_penalty": fr_term, "purity_penalty": purity_penalty,
                "score": score, "okF": okF, "purity_ok": bool(purity_ok),
                "flags": ";".join(flags) if flags else ""
            })

    recs = pd.DataFrame(rows)
    if recs.empty: sys.exit("[ERROR] No pair candidates generated. Check filters/inputs.")

    # whitelist note
    recs["whitelist_pair"] = recs.apply(lambda r: (r.ID1 in whitelist) or (r.ID2 in whitelist), axis=1) if whitelist else False

    # Sort and rank
    jitter = np.random.uniform(0, 1e-9, size=len(recs))
    recs = recs.sort_values(by=["okF","purity_ok","score"], ascending=[False, False, True]).assign(_j=jitter)
    recs = recs.sort_values(by=["okF","purity_ok","score","_j"], ascending=[False, False, True, True]).drop(columns=["_j"]).reset_index(drop=True)
    recs.insert(0, "rank", np.arange(1, len(recs)+1))

    # Pair-level recommendation labels
    elig = recs["okF"] & recs["purity_ok"]
    recs["recommendation"] = "DO NOT BREED"
    if elig.any():
        eligible_scores = recs.loc[elig, "score"]
        pct = eligible_scores.rank(pct=True, method="average")
        breed_cut = args.breed_quantile
        hold_cut  = args.hold_quantile
        lbl = []
        for idx in recs.index:
            if not elig.loc[idx]:
                lbl.append("DO NOT BREED")
            else:
                pr = pct.loc[idx]
                if pr <= breed_cut: lbl.append("BREED")
                elif pr <= hold_cut: lbl.append("HOLD")
                else: lbl.append("DO NOT BREED")
        recs["recommendation"] = lbl

    # Outputs
    prefix = args.prefix
    recs.to_csv(f"{prefix}.recommended_pairs.csv", index=False)

    mk_df = MK.sort_values(ascending=True).reset_index()
    mk_df.columns = ["ID","MK"]
    mk_df.to_csv(f"{prefix}.priority_by_MK.csv", index=False)
    mk_df.sort_values("MK", ascending=False).to_csv(f"{prefix}.suppress_by_MK.csv", index=False)

    if founder_map:
        grp = pd.DataFrame({"ID": list(founder_map.keys()), "FounderGroup": list(founder_map.values())})
        grp = grp.merge(mk_df, how="left", on="ID")
        grp.to_csv(f"{prefix}.founder_groups_with_MK.csv", index=False)

    # Individual "do not pair"
    reasons = {i: [] for i in candidates}
    for i in candidates:
        if i in blacklist: reasons[i].append("blacklisted")
        rows_i = recs[((recs["ID1"] == i) | (recs["ID2"] == i)) & recs["okF"] & recs["purity_ok"]]
        if rows_i.empty: reasons[i].append("no safe partners under constraints")
    top_decile_cut = mk_df["MK"].quantile(0.90) if len(mk_df) >= 10 else mk_df["MK"].max()
    for i, mkv in mk_df.itertuples(index=False):
        if mkv >= top_decile_cut: reasons[i].append("high MK (overrepresented)")
    dnpi = [{"ID": i,
             "Group": taxa_map.get(i, np.nan),
             "Sex": (None if sexes is None else sexes.get(i, None)),
             "MK": MK[i],
             "Reasons": "; ".join(sorted(set(rs)))}
            for i, rs in reasons.items() if rs]
    pd.DataFrame(dnpi).sort_values(["Reasons","MK"], ascending=[True, False]) \
        .to_csv(f"{prefix}.do_not_pair_individuals.csv", index=False)

    # Console
    print(f"[OK] Wrote: {prefix}.recommended_pairs.csv")
    print(f"[OK] Wrote: {prefix}.priority_by_MK.csv (low MK ⇒ prioritize)")
    print(f"[OK] Wrote: {prefix}.suppress_by_MK.csv (high MK ⇒ consider hold/contracept)")
    print(f"[OK] Wrote: {prefix}.do_not_pair_individuals.csv")
    if founder_map:
        print(f"[OK] Wrote: {prefix}.founder_groups_with_MK.csv")
    print("\nTop 10 BREED recommendations:")
    cols = [c for c in ["rank","ID1","sex1","group1","MK1","ID2","sex2","group2","MK2","pair_kinship","F_offspring_pred","score","recommendation","flags"] if c in recs.columns]
    print(recs.loc[recs["recommendation"]=="BREED"].head(10)[cols].to_string(index=False))

if __name__ == "__main__":
    main()
