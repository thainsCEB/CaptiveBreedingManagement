#!/usr/bin/env python3
"""
Classify mitochondrial taxa from a Newick tree using reference individuals,
merge onto a summary table, and produce classification reports.

Outputs
-------
1) Annotated summary (TSV): adds Mito_taxon and Taxon_status, and standardizes Nuclear_taxon/Nuclear_status.
2) Hybrids list:         <OUT>.hybrids.tsv
3) Taxon counts:         <OUT>.taxon_counts.tsv
4) Status counts:        <OUT>.status_counts.tsv
5) Hybrid count:         <OUT>.hybrids.count.txt
6) Combined summary:     <OUT>.counts_summary.tsv  (unifies taxon counts, hybrid count, and status counts)

Notes
-----
- Taxon_status values:
    * "Pure"            (when Nuclear_taxon == Mito_taxon)
    * "hybrid (...)"    (if Nuclear_status indicates hybrid)
    * "mito-nuclear discordant (...)" (different mito vs nuclear on non-hybrids)
    * "uncertain" or "mito_only" for edge cases
- Classification methods:
    * "nearest"       : by patristic distance to nearest reference
    * "lca-majority"  : by majority reference taxon in the smallest clade containing the query + refs
- Optional --name-map to translate tree labels to sample_id.
- Optional --distance-cache to speed repeated runs.

Requirements
------------
pip install biopython pandas
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from typing import Dict, Optional

import pandas as pd

try:
    from Bio import Phylo
    HAVE_BIOPYTHON = True
except Exception:
    HAVE_BIOPYTHON = False


def read_table(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)


def load_name_map(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    nm = read_table(path)
    need = {"tree_label", "sample_id"}
    if not need.issubset(nm.columns):
        raise SystemExit(f"--name-map needs columns {sorted(need)}; got {list(nm.columns)}")
    return dict(zip(nm["tree_label"].astype(str), nm["sample_id"].astype(str)))


def build_parent_map(root):
    parent = {}
    for clade in root.find_clades(order="level"):
        for child in clade.clades:
            parent[child] = clade
    return parent


def lca_majority_assign(tree, leaf_to_taxon: Dict[str, str], query: str) -> Optional[str]:
    node = next(tree.find_clades(name=query), None)
    if node is None:
        return None
    parent = build_parent_map(tree.root)
    cur = node
    while cur is not None:
        leaves = [t.name for t in cur.get_terminals() if t.name]
        taxa = [leaf_to_taxon[l] for l in leaves if l in leaf_to_taxon]
        if taxa:
            most_common = Counter(taxa).most_common()
            top_taxon, top_n = most_common[0]
            if len(most_common) > 1 and most_common[1][1] == top_n:
                return None  # tie -> uncertain
            return top_taxon
        cur = parent.get(cur, None)
    return None


def nearest_reference_assign(tree, leaf_to_taxon: Dict[str, str], query: str,
                             max_distance: Optional[float],
                             dist_cache: Dict[tuple[str, str], float]) -> Optional[str]:
    leaves = {t.name for t in tree.get_terminals() if t.name}
    if query not in leaves:
        return None
    best_taxon = None
    best_dist = None
    for leaf, taxon in leaf_to_taxon.items():
        if leaf == query:
            continue
        key = (query, leaf) if query < leaf else (leaf, query)
        d = dist_cache.get(key)
        if d is None:
            try:
                d = tree.distance(query, leaf)
            except Exception:
                continue
            dist_cache[key] = d
        if best_dist is None or d < best_dist:
            best_dist = d
            best_taxon = taxon
    if best_dist is None:
        return None
    if max_distance is not None and best_dist > max_distance:
        return "uncertain"
    return best_taxon


def main():
    ap = argparse.ArgumentParser(description="Classify mito taxa from Newick + refs; merge and report.")
    ap.add_argument("--tree", required=True, help="Newick tree file (leaf names are individuals)")
    ap.add_argument("--refs", required=True, help="Reference TSV/CSV with columns: sample_id,taxon")
    ap.add_argument("--summary", required=True, help="Existing summary TSV/CSV with at least sample_id")
    ap.add_argument("--out", required=True, help="Output annotated summary TSV path")
    ap.add_argument("--name-map", default=None, help="Optional TSV/CSV with columns: tree_label,sample_id")
    ap.add_argument("--ref-id-col", default="sample_id", help="Column in refs for sample IDs (default: sample_id)")
    ap.add_argument("--ref-taxon-col", default="taxon", help="Column in refs for taxon (default: taxon)")
    ap.add_argument("--method", choices=["nearest", "lca-majority"], default="nearest")
    ap.add_argument("--max-distance", type=float, default=None, help="If set with nearest, far hits -> 'uncertain'")
    ap.add_argument("--distance-cache", default=None, help="JSON file to cache pairwise distances")
    ap.add_argument("--report-prefix", default=None, help="Prefix for report TSVs. Default: use --out")
    args = ap.parse_args()

    if not HAVE_BIOPYTHON:
        raise SystemExit("Biopython is required. Install with: pip install biopython")

    # Load inputs
    tree = Phylo.read(args.tree, "newick")
    refs = read_table(args.refs)
    if args.ref_id_col not in refs.columns or args.ref_taxon_col not in refs.columns:
        raise SystemExit(f"Reference file must have '{args.ref_id_col}' and '{args.ref_taxon_col}'. Got: {list(refs.columns)}")

    name_map = load_name_map(args.name_map)
    inv_name_map = {v: k for k, v in name_map.items()}

    # Reference leaf->taxon map
    leaf_to_taxon: Dict[str, str] = {}
    for _, row in refs.iterrows():
        sid = str(row[args.ref_id_col])
        tax = str(row[args.ref_taxon_col])
        leaf = inv_name_map.get(sid, sid)
        leaf_to_taxon[leaf] = tax

    leaves = [t.name for t in tree.get_terminals() if t.name]

    # Distance cache
    dist_cache: Dict[tuple[str, str], float] = {}
    if args.distance_cache and os.path.exists(args.distance_cache):
        try:
            with open(args.distance_cache, "r") as fh:
                raw = json.load(fh)
            for k, v in raw.items():
                a, b = k.split("||", 1)
                dist_cache[(a, b)] = v
        except Exception:
            pass

    # Assign mito taxa
    mito_assignments: Dict[str, str] = {}
    for leaf in leaves:
        if leaf in leaf_to_taxon:
            mito_assignments[leaf] = leaf_to_taxon[leaf]
            continue
        if args.method == "nearest":
            tax = nearest_reference_assign(tree, leaf_to_taxon, leaf, args.max_distance, dist_cache)
        else:
            tax = lca_majority_assign(tree, leaf_to_taxon, leaf)
        mito_assignments[leaf] = tax if tax is not None else "uncertain"

    if args.distance_cache:
        serial = {"||".join(k): v for k, v in dist_cache.items()}
        with open(args.distance_cache, "w") as fh:
            json.dump(serial, fh, indent=2)

    # Load summary and standardize columns
    summ = read_table(args.summary)
    if "sample_id" not in summ.columns:
        # case-insensitive rescue
        for c in list(summ.columns):
            if c.lower() == "sample_id":
                summ = summ.rename(columns={c: "sample_id"})
                break
        else:
            raise SystemExit("summary must have a 'sample_id' column")

    # Standardize Nuclear columns from legacy names
    if "Nuclear_taxon" not in summ.columns and "Putative_taxon" in summ.columns:
        summ = summ.rename(columns={"Putative_taxon": "Nuclear_taxon"})
    if "Nuclear_status" not in summ.columns and "Status" in summ.columns:
        summ = summ.rename(columns={"Status": "Nuclear_status"})

    # Map mito assignments to sample IDs (via name_map if provided)
    mito_to_sample = { (name_map.get(leaf, leaf)): tax for leaf, tax in mito_assignments.items() }
    summ["Mito_taxon"] = summ["sample_id"].astype(str).map(mito_to_sample).fillna("absent_in_tree")

    # Build Taxon_status
    in_nuc_col = "Nuclear_taxon" if "Nuclear_taxon" in summ.columns else None
    in_status_col = "Nuclear_status" if "Nuclear_status" in summ.columns else None

    def compute_status(row) -> str:
        mito = str(row["Mito_taxon"])
        if mito in ("uncertain", "absent_in_tree", "None", "nan"):
            return "uncertain"
        nuc = str(row[in_nuc_col]) if in_nuc_col else None
        status_val = str(row[in_status_col]).lower() if in_status_col else ""
        if nuc is None or nuc == "None" or nuc == "nan":
            return "mito_only"
        if nuc == mito:
            return "Pure"
        if "hybrid" in status_val:
            return f"hybrid (nuclear: {nuc}, mito: {mito})"
        return f"mito-nuclear discordant (nuclear: {nuc}, mito: {mito})"

    summ["Taxon_status"] = summ.apply(compute_status, axis=1)

    # Write annotated summary
    outdir = os.path.dirname(os.path.abspath(args.out)) or "."
    os.makedirs(outdir, exist_ok=True)
    summ.to_csv(args.out, sep="\t", index=False)

    # === Reporting ===
    rep_prefix = args.report_prefix or args.out

    # (A) Taxon counts (by mito)
    taxon_counts = (
        summ["Mito_taxon"]
        .fillna("uncertain")
        .value_counts(dropna=False)
        .rename_axis("taxon")
        .reset_index(name="count")
    )
    taxon_counts.to_csv(f"{rep_prefix}.taxon_counts.tsv", sep="\t", index=False)

    # (B) Hybrids table and count
    is_hybrid = summ["Taxon_status"].str.startswith("hybrid", na=False) if "Taxon_status" in summ.columns else pd.Series(False, index=summ.index)
    hybrid_cols = ["sample_id"]
    if "Nuclear_taxon" in summ.columns:
        hybrid_cols.append("Nuclear_taxon")
    hybrid_cols += ["Mito_taxon", "Taxon_status"]
    hybrids = summ.loc[is_hybrid, hybrid_cols].copy()
    hybrids.to_csv(f"{rep_prefix}.hybrids.tsv", sep="\t", index=False)
    hybrid_count = int(len(hybrids))
    with open(f"{rep_prefix}.hybrids.count.txt", "w") as fh:
        fh.write(str(hybrid_count) + "\n")

    # (C) Status counts (bucketed)
    def status_bucket(s: str) -> str:
        s = str(s)
        if s == "Pure":
            return "Pure"
        if s.startswith("hybrid"):
            return "hybrid"
        if s.startswith("mito-nuclear discordant"):
            return "mito-nuclear discordant"
        if s == "uncertain":
            return "uncertain"
        if s == "mito_only":
            return "mito_only"
        return "other"

    status_counts = (
        summ["Taxon_status"].map(status_bucket)
        .value_counts(dropna=False)
        .rename_axis("status")
        .reset_index(name="count")
    )
    status_counts.to_csv(f"{rep_prefix}.status_counts.tsv", sep="\t", index=False)

    # (D) Combined counts summary table
    combined_rows = []
    for _, r in taxon_counts.iterrows():
        combined_rows.append({"section": "taxon_counts", "label": str(r["taxon"]), "count": int(r["count"])})

    combined_rows.append({"section": "hybrid_count", "label": "hybrid", "count": hybrid_count})

    for _, r in status_counts.iterrows():
        combined_rows.append({"section": "status_counts", "label": str(r["status"]), "count": int(r["count"])})

    combined_df = pd.DataFrame(combined_rows, columns=["section", "label", "count"])
    combined_df.to_csv(f"{rep_prefix}.counts_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
