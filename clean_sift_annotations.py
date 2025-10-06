#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clean SIFT annotations in VCF and table outputs.

- VCF: remove the INFO key (default: SIFTINFO) if its value mentions UTR or NONCODING.
       If the entire INFO field becomes empty, write '.'.

- Table (.xlsx/.tsv/.csv): by default, BLANK SIFT-related columns on rows where
       any of those columns include UTR or NONCODING (case-insensitive).
       Optionally, drop such rows instead (--table-mode drop-rows).

Example:
  python clean_sift_annotations.py \
      --vcf-in input.vcf.gz \
      --vcf-out input.clean.vcf.gz \
      --table-in SIFT_annot.xlsx \
      --table-out SIFT_annot_clean.xlsx

Author: ChatGPT
"""
import argparse
import gzip
import os
import re
import sys
from typing import List, Optional

def open_by_ext(path: str, mode: str):
    """
    Open a file with transparent gzip handling based on extension.
    mode like 'rt' or 'wt' for text; 'rb' or 'wb' for binary.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding=None if "b" in mode else "utf-8")

def clean_vcf_info_field(info: str, info_key: str, pattern: re.Pattern) -> str:
    """
    Remove entries of `info_key` from INFO if its value matches pattern.
    Keep other INFO keys intact. If resulting INFO is empty, return '.'.
    INFO format is semicolon-separated; key=value (or flag keys without '=')
    """
    if info == "." or not info:
        return info

    parts = info.split(";")
    kept: List[str] = []
    for p in parts:
        if not p:
            continue
        if "=" in p:
            k, v = p.split("=", 1)
        else:
            k, v = p, None

        if k == info_key and v is not None:
            # Drop this key ONLY if its value matches the pattern
            if pattern.search(v):
                continue  # skip it
            else:
                kept.append(p)
        else:
            kept.append(p)

    if not kept:
        return "."
    return ";".join(kept)

def clean_table_grep_like(table_in: str, table_out: str, pattern: str):
    import gzip, io, os, sys, re
    patt = re.compile(pattern)
    # Decide opener based on extension
    def _open(path, mode):
        return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode, encoding='utf-8', newline='')
    # Read line-wise and write out lines that DO NOT match pattern (grep -v)
    with _open(table_in, 'rt') as fin, _open(table_out, 'wt') as fout:
        removed = 0
        kept = 0
        for line in fin:
            if patt.search(line):
                removed += 1
                continue
            fout.write(line)
            kept += 1
    print(f"[INFO] Grep-like table clean: kept={kept} removed={removed} pattern={pattern}", file=sys.stderr)

def clean_vcf(in_path: str, out_path: str, info_key: str, regex: str, drop_low_confidence: bool = True, lowconf_pattern: str = "(?:low\\s*confidence)"):
    patt = re.compile(regex, flags=re.IGNORECASE)
    lowp = re.compile(lowconf_pattern, flags=re.IGNORECASE) if drop_low_confidence else None
    with open_by_ext(in_path, "rt") as fin, open_by_ext(out_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                # malformed line, write as-is
                fout.write(line)
                continue
            info = fields[7]
            # First remove UTR/NONCODING SIFT entries
            info = clean_vcf_info_field(info, info_key, patt)
            # Then remove low-confidence SIFT entries (do not drop the variant)
            if lowp is not None and info not in (".", ""):
                parts = info.split(";")
                kept = []
                for p in parts:
                    if "=" in p:
                        k, v = p.split("=", 1)
                    else:
                        k, v = p, None
                    if k == info_key and v is not None and lowp.search(v):
                        continue
                    kept.append(p)
                info = "." if not kept else ";".join(kept)
            fields[7] = info
            fout.write("\t".join(fields) + "\n")

def _read_table_auto(path: str, no_header: bool = False):
    """
    Read .xlsx/.csv/.tsv/.txt into (df, kind, sep)
      kind: 'xlsx' or 'delim'
      sep:  delimiter for 'delim' (',' or '\t')
    """
    import pandas as pd
    lower = path.lower()
    if lower.endswith(".xlsx"):
        df = pd.read_excel(path, header=None if no_header else 0)
        return df, "xlsx", None
    # default to delimiter-based
    # Choose tab if .tsv/.txt else comma for .csv
    if lower.endswith(".tsv") or lower.endswith(".txt"):
        df = pd.read_csv(path, sep="\t", dtype=str, header=None if no_header else 0, keep_default_na=False)
        return df, "delim", "\t"
    # csv fallback
    df = pd.read_csv(path, sep=",", dtype=str, header=None if no_header else 0, keep_default_na=False)
    return df, "delim", ","

def _write_table_auto(df, out_path: str, kind: str, sep: Optional[str]):
    lower = out_path.lower()
    if kind == "xlsx" or lower.endswith(".xlsx"):
        # ensure xlsx writer
        df.to_excel(out_path, index=False)
    else:
        # delimiter-based
        if lower.endswith(".tsv") or (sep == "\t" and not lower.endswith(".csv")):
            df.to_csv(out_path, sep="\t", index=False)
        else:
            df.to_csv(out_path, sep=",", index=False)

def clean_table(table_in: str,
                table_out: str,
                info_key: str,
                regex: str,
                table_mode: str = "blank",
                sift_columns: Optional[List[str]] = None,
                debug: bool = False,
                no_header: bool = False,
                match_cols: Optional[List[int]] = None,
                blank_cols: Optional[List[int]] = None,
                fix_scaffold_typos: bool = False,
                drop_all_na_annotations: bool = False,
                drop_empty_lines: bool = True,
                drop_low_confidence: bool = True,
                lowconf_pattern: str = "(?:low\\s*confidence)"):
    """
    Clean the SIFT table:
      - Identify SIFT-related columns (names containing 'SIFT' or matching `info_key`),
        or use `sift_columns` if provided.
      - Find rows where ANY of those columns match regex (case-insensitive).
      - If table_mode == 'blank': set those SIFT columns to empty string on matched rows.
        If table_mode == 'drop-rows': drop matched rows entirely.
    """
    import pandas as pd

    df, kind, sep = _read_table_auto(table_in, no_header=no_header)

    import pandas as pd, re as _re

    # Helper: NA-like
    def _na_like_series(s):
        # Treat true NaNs as NA-like *before* casting to string
        na_mask = s.isna()
        s2 = s.astype(str).str.strip()
        # Also consider common NA tokens and literal doubled-quotes
        tokens = {"", ".", "NA", "N/A", "NULL", "None", "nan", "NaN", '""'}
        tok_mask = s2.str.upper().isin({t.upper() for t in tokens})
        return na_mask | tok_mask

    # Drop rows that are entirely NA-like across ALL columns (handles '""' lines)
    if drop_empty_lines:
        all_na_mask = None
        for c in df.columns:
            col_na = _na_like_series(df[c])
            all_na_mask = col_na if all_na_mask is None else (all_na_mask & col_na)
        # Additional explicit grep-like rule: drop rows that are exactly '""' (possibly with whitespace)
        single_quote_mask = None
        if df.shape[1] >= 1:
            # Build a mask that's True only if ALL columns are empty AND
            # any column literally equals '""' (after strip), which captures raw lines like '""'
            # when parsed into a single column.
            literal_mask = df.apply(lambda row: any(str(v).strip() == '""' for v in row), axis=1)
            single_quote_mask = literal_mask
            all_na_mask = all_na_mask | single_quote_mask
        if all_na_mask is not None:
            n_empty = int(all_na_mask.sum())
            if debug and n_empty:
                print(f"[DEBUG] Dropping {n_empty} completely empty/NA-like or literal \"\" row(s).", file=sys.stderr)
            if n_empty:
                df = df.loc[~all_na_mask].copy()

    # Optionally drop low-confidence rows up-front (any column contains lowconf pattern)
    if drop_low_confidence:
        lowp = _re.compile(lowconf_pattern, _re.IGNORECASE)
        low_mask = None
        for c in df.columns:
            s = df[c].astype(str)
            m = s.str.contains(lowp, regex=True, na=False)
            low_mask = m if low_mask is None else (low_mask | m)
        n_low = int(low_mask.sum()) if low_mask is not None else 0
        if debug and n_low:
            print(f"[DEBUG] Dropping {n_low} row(s) with low-confidence predictions.", file=sys.stderr)
        if n_low:
            df = df.loc[~low_mask].copy()

    import pandas as pd

    # Assign synthetic headers if requested
    if no_header:
        df.columns = [f"c{i+1}" for i in range(df.shape[1])]
        if debug:
            print(f"[DEBUG] Assigned synthetic headers: {list(df.columns)}", file=sys.stderr)

    # Optional scaffold typo fix in first column
    if fix_scaffold_typos and df.shape[1] >= 1:
        first_col = df.columns[0]
        before_sample = df[first_col].head(3).tolist()
        # Fix patterns: "ffolds_123" -> "Scaffolds_123"; "scaffolds_123" -> "Scaffolds_123"
        df[first_col] = df[first_col].astype(str).str.replace(r'^[Ff]?folds_(\d+)$', r'Scaffolds_\1', regex=True)
        df[first_col] = df[first_col].str.replace(r'^[sS]caffolds_(\d+)$', r'Scaffolds_\1', regex=True)
        after_sample = df[first_col].head(3).tolist()
        if debug:
            print(f"[DEBUG] First column before fix (sample): {before_sample}", file=sys.stderr)
            print(f"[DEBUG] First column after fix  (sample): {after_sample}", file=sys.stderr)

    if match_cols:
        # Convert 1-based indices to names and validate
        total = df.shape[1]
        sel = []
        for idx in match_cols:
            if 1 <= idx <= total:
                sel.append(df.columns[idx-1])
            else:
                print(f"[WARN] --match-cols index {idx} out of range (1..{total})", file=sys.stderr)
        cols = sel
        if not cols:
            print(f"[WARN] No valid --match-cols provided; nothing to clean.", file=sys.stderr)
            _write_table_auto(df, table_out, kind, sep)
            return
    elif sift_columns:
        cols = [c for c in sift_columns if c in df.columns]
        if not cols:
            print(f"[WARN] None of the provided --sift-columns exist in the table. No changes made.", file=sys.stderr)
            _write_table_auto(df, table_out, kind, sep)
            return
    else:
        # auto-detect SIFT columns: contain 'sift' (case-insensitive) OR exact match to info_key
        cols = [c for c in df.columns if ("sift" in c.lower()) or (c == info_key)]
        # also consider generic annotation columns that might hold NONCODING/UTR flags
        for extra in ("Annotation", "Consequence", "Region", "Effect"):
            if extra in df.columns and extra not in cols:
                cols.append(extra)
        if not cols:
            print(f"[WARN] No SIFT-like columns found; nothing to clean. Writing table unchanged.", file=sys.stderr)
            _write_table_auto(df, table_out, kind, sep)
            return
    if debug:
        print(f"[DEBUG] Candidate columns used for matching: {cols}", file=sys.stderr)

    patt = re.compile(regex, flags=re.IGNORECASE)

    # Build mask: any of the chosen columns contains the pattern
    def col_matches(series):
        # Ensure string; handle NA safely
        s = series.astype(str)
        return s.str.contains(patt, na=False, regex=True)

    mask = None
    for c in cols:
        m = col_matches(df[c])
        if debug:
            try:
                cnt = int(m.sum())
            except Exception:
                cnt = -1
            print(f"[DEBUG] Matches in column '{c}': {cnt}", file=sys.stderr)
        mask = m if mask is None else (mask | m)

    # Determine which columns to blank when matched
    if blank_cols:
        total = df.shape[1]
        blank_names = []
        for idx in blank_cols:
            if 1 <= idx <= total:
                blank_names.append(df.columns[idx-1])
            else:
                print(f"[WARN] --blank-cols index {idx} out of range (1..{total})", file=sys.stderr)
        if not blank_names:
            blank_names = cols[:]  # fallback
    else:
        blank_names = cols[:]

    if debug:
        print(f"[DEBUG] Columns that will be blanked on match: {blank_names}", file=sys.stderr)


    affected = int(mask.sum()) if mask is not None else 0
    print(f"[INFO] Table rows affected: {affected}", file=sys.stderr)

    if affected == 0:
        # Even if no direct matches, we may still drop all-NA rows if requested (rare but possible)
        df2 = df.copy()
    else:
        if table_mode == "drop-rows":
            df2 = df.loc[~mask].copy()
        else:
            df2 = df.copy()
            for c in blank_names:
                df2.loc[mask, c] = ""  # blank out

    # Optionally drop rows with all NA/empty values across chosen annotation columns
    if drop_all_na_annotations:
        # Decide which columns to check: prefer explicit blank_names; else cols
        check_cols = blank_names if len(blank_names) > 0 else cols
        # Build a boolean frame of "is NA-like"
        def _is_na_like(s):
            na_mask = s.isna()
            s2 = s.astype(str).str.strip()
            tokens = {"", ".", "NA", "N/A", "NULL", "None", "nan", "NaN", '""'}
            tok_mask = s2.str.upper().isin({t.upper() for t in tokens})
            return na_mask | tok_mask
        na_mask = s.isna()
        s2 = s.astype(str).str.strip()
        tokens = {"", ".", "NA", "N/A", "NULL", "None", "nan", "NaN", '""'}
        tok_mask = s2.str.upper().isin({t.upper() for t in tokens})
        return na_mask | tok_mask
        na_matrix = None
        for c in check_cols:
            series_na = _is_na_like(df2[c])
            na_matrix = series_na if na_matrix is None else (na_matrix & series_na)
        to_drop = na_matrix if na_matrix is not None else None
        if to_drop is not None:
            n_drop = int(to_drop.sum())
            if debug:
                print(f"[DEBUG] Rows to drop because all selected annotation columns are NA/empty: {n_drop}", file=sys.stderr)
            if n_drop > 0:
                df2 = df2.loc[~to_drop].copy()

    
    # Final sweep: drop any row that still contains a literal '""' token in any cell
    try:
        qq_any_mask = df2.apply(lambda row: row.astype(str).str.strip().eq('""').any(), axis=1)
        n_qq = int(qq_any_mask.sum())
        if debug and n_qq:
            print(f"[DEBUG] Dropping {n_qq} row(s) containing literal \"\" tokens before write.", file=sys.stderr)
        if n_qq:
            df2 = df2.loc[~qq_any_mask].copy()
    except Exception as _e:
        if debug:
            print(f"[DEBUG] Final '\"\"' sweep skipped due to error: {_e}", file=sys.stderr)

    _write_table_auto(df2, table_out, kind, sep)

def main():
    p = argparse.ArgumentParser(
        description="Clean SIFT annotations in VCF and tables by removing/blanking entries marked as UTR or NONCODING."
    )
    p.add_argument("--vcf-in", help="Input VCF(.vcf or .vcf.gz) annotated with SIFTINFO in INFO (optional).")
    p.add_argument("--vcf-out", help="Output VCF path. Defaults to <vcf-in>.clean.vcf[.gz]")
    p.add_argument("--table-in", help="SIFT table (.xlsx/.tsv/.csv) to clean (optional).")
    p.add_argument("--table-grep", action="store_true", default=True,
                   help="Use literal grep-like filtering for tables (default: True).")
    p.add_argument("--table-grep-pattern", default=r"UTR|Low confidence",
                   help="Pattern for grep -v style exclusion (default: 'UTR|Low confidence').")
    p.add_argument("--no-table-grep", action="store_false", dest="table_grep",
                   help="Disable grep-like filtering and use advanced table cleaning instead.")
    p.add_argument("--table-out", help="Output table path. Defaults to add _clean before extension.")
    p.add_argument("--info-key", default="SIFTINFO", help="INFO key for SIFT annotation in VCF (default: SIFTINFO).")
    p.add_argument("--pattern", default=r"(?:UTR|NONCODING)", help="Regex to match entries to remove/blank (case-insensitive).")
    p.add_argument("--table-mode", choices=["blank","drop-rows"], default="blank",
                   help="How to clean table rows that match the pattern in SIFT columns (default: blank).")
    p.add_argument("--sift-columns", nargs="+",
                   help="Explicit list of column names in the table to treat as SIFT columns (optional).")
    p.add_argument("--no-header", action="store_true",
                   help="Treat the input table as headerless; synthetic headers c1..cN will be assigned.")
    p.add_argument("--match-cols", nargs="+", type=int,
                   help="1-based indices of columns to search for regex matches (used when no/unknown headers).")
    p.add_argument("--blank-cols", nargs="+", type=int,
                   help="1-based indices of columns to blank when a row matches (default: the same as --match-cols or detected SIFT columns).")
    p.add_argument("--fix-scaffold-typos", action="store_true",
                   help="Fix common scaffold typos in first column (e.g., ffolds_123 -> Scaffolds_123).")
    p.add_argument("--drop-all-na-annotations", action="store_true",
                   help="After cleaning, drop rows where all selected annotation columns are NA/empty/\".\".")
    p.add_argument("--drop-empty-lines", action="store_true", default=True,
                   help="Drop rows that are entirely empty/NA-like (default: True).")
    p.add_argument("--keep-empty-lines", action="store_false", dest="drop_empty_lines",
                   help="Keep rows that look empty/NA-like (overrides default).")
    p.add_argument("--drop-low-confidence", action="store_true", default=True,
                   help="Drop rows (table) and remove SIFTINFO entries (VCF) that include low-confidence predictions (default: True).")
    p.add_argument("--lowconf-pattern", default=r"(?:low\s*confidence)",
                   help="Regex used to detect low-confidence predictions (case-insensitive). Default matches \"low confidence\".")
    p.add_argument("--debug", action="store_true", help="Verbose reporting of detected columns and match counts.")
    args = p.parse_args()

    if not args.vcf_in and not args.table_in:
        print("[ERROR] Provide at least --vcf-in or --table-in", file=sys.stderr)
        sys.exit(2)

    # VCF cleaning
    if args.vcf_in:
        vcf_out = args.vcf_out
        if not vcf_out:
            base = args.vcf_in
            if base.endswith(".vcf.gz"):
                vcf_out = base.replace(".vcf.gz", ".clean.vcf.gz")
            elif base.endswith(".vcf"):
                vcf_out = base.replace(".vcf", ".clean.vcf")
            else:
                vcf_out = base + ".clean.vcf"
        print(f"[INFO] Cleaning VCF: {args.vcf_in} -> {vcf_out}", file=sys.stderr)
        clean_vcf(args.vcf_in, vcf_out, args.info_key, args.pattern, args.drop_low_confidence, args.lowconf_pattern)

    # Table cleaning
    if args.table_in:
        tbl_out = args.table_out
        if not tbl_out:
            root, ext = os.path.splitext(args.table_in)
            tbl_out = f"{root}_clean{ext or '.tsv'}"
        print(f"[INFO] Cleaning table: {args.table_in} -> {tbl_out}", file=sys.stderr)
        if args.table_grep:
            clean_table_grep_like(args.table_in, tbl_out, args.table_grep_pattern)
        else:
            clean_table(args.table_in, tbl_out, args.info_key, args.pattern, args.table_mode, args.sift_columns, args.debug, args.no_header, args.match_cols, args.blank_cols, args.fix_scaffold_typos, args.drop_all_na_annotations, args.drop_empty_lines, args.drop_low_confidence, args.lowconf_pattern)

if __name__ == "__main__":
    main()
