#!/usr/bin/env python3

"""
split_beagle.py

Split an ANGSD-style .beagle[.gz] file into per-individual BEAGLEs.

Header behavior:
- For each split file, write a header:
    marker  allele1  allele2  IndN  IndN  IndN
  where N is the 0-based index of that individual in BEAGLE order.

Data rows:
- Each data row keeps the original three GL values (triplet) for that individual.

Filenames:
- If --samples-list/-s is provided, those names are used for *filenames only* (in order).
- If not provided, filenames default to Ind0, Ind1, ...
- Header column names remain IndN repeated 3x regardless.

CLI:
  split_beagle.py -i cohort.beagle.gz -x -s samples.txt -o split_beagles
"""

import argparse
import gzip
import os
import re
import sys
from typing import IO, List, Optional


def smart_open(path: str, mode: str = 'rt') -> IO:
    """Open plain or gz file transparently."""
    if path.endswith('.gz'):
        return gzip.open(path, mode)  # type: ignore
    return open(path, mode)  # type: ignore


def infer_n_individuals_from_data(first_data_line: str) -> int:
    """Infer number of individuals N from a data line with 3 + 3*N columns."""
    parts = first_data_line.rstrip('\n').split()
    if len(parts) < 6:
        raise ValueError("Data row has too few columns to be a valid BEAGLE line.")
    n_cols = len(parts)
    if (n_cols - 3) % 3 != 0:
        raise ValueError(f"Columns count {n_cols} is not 3 + 3*N; cannot infer N.")
    return (n_cols - 3) // 3


def load_samples_list(path: str) -> List[str]:
    """Load one sample name per line, stripping whitespace; ignore empty lines."""
    names = []
    with open(path, 'rt') as fh:
        for line in fh:
            name = line.strip()
            if name:
                names.append(name)
    if not names:
        raise ValueError("Provided --samples-list/-s is empty.")
    return names


def sanitize_filename(name: str) -> str:
    """Make a safe filename component from a sample name."""
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name).strip('_') or "Ind"


def split_beagle_by_individual(in_beagle: str, out_dir: str, samples_list: Optional[str] = None) -> None:
    """
    Split a BEAGLE into per-individual BEAGLEs.

    Output files are named using the provided samples list (if any) or IndN otherwise.
    Header in each output: marker allele1 allele2 IndN IndN IndN
    """
    os.makedirs(out_dir, exist_ok=True)

    with smart_open(in_beagle, 'rt') as fh:
        # Read first non-empty line
        first_line = fh.readline()
        while first_line and not first_line.strip():
            first_line = fh.readline()
        if not first_line:
            raise ValueError("Input BEAGLE is empty.")

        # Find first data line (skip any header-like lines starting with 'marker')
        data_line = first_line
        while data_line and (not data_line.strip() or data_line.split()[0].lower() == 'marker'):
            data_line = fh.readline()
        if not data_line:
            raise ValueError("No data lines found in BEAGLE.")

        n_ind = infer_n_individuals_from_data(data_line)

        # Build names for filenames
        if samples_list:
            file_names = [sanitize_filename(nm) for nm in load_samples_list(samples_list)]
            if len(file_names) != n_ind:
                raise ValueError(f"--samples-list/-s has {len(file_names)} names but BEAGLE has {n_ind} individuals.")
        else:
            file_names = [f"Ind{i}" for i in range(n_ind)]

        # Prepare writers
        writers: List[IO] = []
        out_paths: List[str] = []
        try:
            for i in range(n_ind):
                ind_id = f"Ind{i}"
                fname = file_names[i]
                out_path = os.path.join(out_dir, f"{fname}.beagle.gz")
                out_paths.append(out_path)
                w = gzip.open(out_path, 'wt')
                writers.append(w)
                # Header: marker allele1 allele2 IndN IndN IndN
                w.write(f"marker\tallele1\tallele2\t{ind_id}\t{ind_id}\t{ind_id}\n")
            # Stream through file from the beginning
            fh.seek(0)
            for line in fh:
                if not line.strip():
                    continue
                # Skip header-ish lines
                if line.split()[0].lower() == 'marker':
                    continue
                parts = line.rstrip('\n').split()
                if len(parts) < 3 + 3 * n_ind:
                    raise ValueError("Encountered a line with too few columns for expected N.")
                fixed = parts[:3]
                gls = parts[3:]
                for i in range(n_ind):
                    triplet = gls[3*i:3*(i+1)]
                    writers[i].write('\t'.join(fixed + triplet) + '\n')
        finally:
            for w in writers:
                try:
                    w.close()
                except Exception:
                    pass

    # Report
    sys.stdout.write(f"Split {in_beagle} into {n_ind} files in {out_dir}:\n")
    for fn in out_paths:
        sys.stdout.write(f"  - {fn}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Split ANGSD-style BEAGLE into per-individual files."
    )
    ap.add_argument("-i", "--in-beagle", required=True, help="Input .beagle or .beagle.gz")
    ap.add_argument("-o", "--out-dir", default="out_beagle", help="Output directory")
    ap.add_argument("-x", "--split", action="store_true", help="Split multi-individual BEAGLE into per-individual files")
    ap.add_argument("-s", "--samples-list", help="Text file with sample names (used for filenames only)")
    args = ap.parse_args()

    if args.split:
        split_beagle_by_individual(args.in_beagle, args.out_dir, args.samples_list)
    else:
        ap.error("No action specified. Add --split/-x to split per individual.")


if __name__ == "__main__":
    main()
