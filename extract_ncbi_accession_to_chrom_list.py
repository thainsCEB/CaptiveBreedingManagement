#!/usr/bin/env python3
# extract_ncbi_accession_to_chrom_map.py
# Description: Extract accession->chromosome mapping from an NCBI assembly report
#              (preferred) or FASTA file, with flexible category-based filtering
#              (mt, autosomes, sex, autosomes_sex, unlocalized, unplaced, placed_*).
# Author: Taylor Hains
# Date: 2025-11-13

import argparse
import gzip
import re
import sys
from typing import Optional, List

# Category choices for --keep / --remove
CATEGORY_CHOICES = {
    "mt",
    "autosomes",
    "sex",
    "autosomes_sex",
    "unlocalized",
    "unplaced",
    "placed_plus_mt",
    "placed_no_mt",
}

SEX_NAMES = {"X", "Y", "Z", "W"}


def _open_maybe_gzip(path: str):
    if path == "-" or path is None:
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def parse_assembly_report(path: str):
    rows = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "sequence_name": parts[0],
                    "sequence_role": parts[1],
                    "assigned_molecule": parts[2],
                    "assigned_molecule_location_type": parts[3],
                    "genbank_accn": parts[4],
                    "refseq_accn": parts[5],
                }
            )
    return rows


def classify_from_report_row(r: dict) -> str:
    role = r["sequence_role"].lower()
    mol = (
        r["assigned_molecule"]
        .upper()
        .replace("CHR", "")
        .replace("CHROMOSOME", "")
        .strip()
    )
    if "mitochondr" in role or mol in {"MT", "M"}:
        return "mt"
    if "unlocalized" in role:
        return "unlocalized"
    if "unplaced" in role:
        return "unplaced"
    if role == "chromosome":
        if mol in SEX_NAMES:
            return "sex"
        if re.fullmatch(r"[0-9]+", mol) or re.fullmatch(r"[IVXLCDM]+", mol):
            return "autosomes"
        return "other_chromosome"
    return "other"


def ucsc_name(mol: str) -> str:
    m = mol.upper()
    return "chrM" if m in {"MT", "M"} else f"chr{m}"


def best_accession(refseq_accn: str, genbank_accn: str) -> Optional[str]:
    if refseq_accn and refseq_accn.lower() != "na":
        return refseq_accn
    if genbank_accn and genbank_accn.lower() != "na":
        return genbank_accn
    return None


FASTA_ACC_RE = re.compile(r"^>(\S+)")
RE_CHROM = re.compile(r"\bchromosome\s+([A-Za-z0-9]+)\b", re.IGNORECASE)
RE_MT = re.compile(r"\bmitochondr", re.IGNORECASE)
RE_UNLOCALIZED = re.compile(r"\bunlocalized\b", re.IGNORECASE)
RE_UNPLACED = re.compile(r"\bunplaced\b", re.IGNORECASE)


def parse_fasta_headers(path: str):
    res = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = FASTA_ACC_RE.match(line)
            if not m:
                continue
            acc = m.group(1)
            desc = line[1:].strip()
            category = "other"
            molecule = None

            if RE_MT.search(desc):
                category = "mt"
                molecule = "MT"
            elif RE_UNLOCALIZED.search(desc):
                category = "unlocalized"
            elif RE_UNPLACED.search(desc):
                category = "unplaced"

            chrom_m = RE_CHROM.search(desc)
            if chrom_m:
                molecule = chrom_m.group(1)
                up = molecule.upper()
                if up in SEX_NAMES:
                    category = "sex"
                elif up in {"MT", "M"}:
                    category = "mt"
                elif re.fullmatch(r"[0-9]+", molecule) or re.fullmatch(
                    r"[IVXLCDM]+", molecule
                ):
                    category = "autosomes"
                else:
                    category = "other_chromosome"

            res.append(
                {
                    "accession": acc,
                    "molecule": molecule,
                    "category": category,
                    "sequence_name": acc,
                }
            )
    return res


def matches_category(rec_cat: str, token: str, is_placed_chrom: bool) -> bool:
    """
    Decide if a record with category `rec_cat` matches the filter `token`.

    token is one of CATEGORY_CHOICES:
      - mt
      - autosomes
      - sex
      - autosomes_sex  (NEW: both autosomes and sex chromosomes)
      - unlocalized
      - unplaced
      - placed_plus_mt
      - placed_no_mt
    """
    if token == "mt":
        return rec_cat == "mt"
    if token == "autosomes":
        return rec_cat == "autosomes"
    if token == "sex":
        return rec_cat == "sex"
    if token == "autosomes_sex":
        return rec_cat in {"autosomes", "sex"}
    if token == "unlocalized":
        return rec_cat == "unlocalized"
    if token == "unplaced":
        return rec_cat == "unplaced"
    if token == "placed_plus_mt":
        return is_placed_chrom or rec_cat == "mt"
    if token == "placed_no_mt":
        return is_placed_chrom and rec_cat != "mt"
    return False


def parse_category_list(arg: Optional[str]) -> List[str]:
    if not arg:
        return []
    toks = [t.strip().lower() for t in arg.split(",") if t.strip()]
    for t in toks:
        if t not in CATEGORY_CHOICES:
            raise ValueError(
                f"Unknown category '{t}'. Choose from: {', '.join(sorted(CATEGORY_CHOICES))}"
            )
    return toks


def build_map_from_assembly_report(path, style, keep, remove):
    rows = parse_assembly_report(path)
    keep = keep or []
    remove = remove or []
    out = []
    counts = {
        "mt": 0,
        "autosomes": 0,
        "sex": 0,
        "unlocalized": 0,
        "unplaced": 0,
        "other_chromosome": 0,
        "other": 0,
    }
    for r in rows:
        cat = classify_from_report_row(r)
        counts[cat] = counts.get(cat, 0) + 1
        is_placed = r["sequence_role"].lower() == "chromosome"

        include = True if not keep else any(
            matches_category(cat, k, is_placed) for k in keep
        )
        if include and remove and any(
            matches_category(cat, k, is_placed) for k in remove
        ):
            include = False
        if not include:
            continue

        acc = best_accession(r["refseq_accn"], r["genbank_accn"])
        if not acc:
            continue

        assigned = (r["assigned_molecule"] or r["sequence_name"]).strip()
        if style == "raw":
            chrom = assigned
        elif style == "ucsc":
            if r["sequence_role"].lower() == "chromosome" or cat == "mt":
                chrom = ucsc_name(assigned)
            else:
                chrom = r["sequence_name"] or acc
        else:
            raise ValueError("Unknown --style")

        out.append((acc, chrom, cat))
    return out, counts


def build_map_from_fasta(path, style, keep, remove):
    entries = parse_fasta_headers(path)
    keep = keep or []
    remove = remove or []
    out = []
    counts = {
        "mt": 0,
        "autosomes": 0,
        "sex": 0,
        "unlocalized": 0,
        "unplaced": 0,
        "other_chromosome": 0,
        "other": 0,
    }
    for e in entries:
        cat = e["category"]
        counts[cat] = counts.get(cat, 0) + 1
        is_placed = cat in {"autosomes", "sex", "other_chromosome"}

        include = True if not keep else any(
            matches_category(cat, k, is_placed) for k in keep
        )
        if include and remove and any(
            matches_category(cat, k, is_placed) for k in remove
        ):
            include = False
        if not include:
            continue

        acc = e["accession"]
        mol = e["molecule"] or e["sequence_name"]
        if style == "raw":
            chrom = mol
        elif style == "ucsc":
            if cat in {"autosomes", "sex", "mt", "other_chromosome"} and e["molecule"]:
                chrom = ucsc_name(e["molecule"])
            else:
                chrom = e["sequence_name"]
        else:
            raise ValueError("Unknown --style")

        out.append((acc, chrom, cat))
    return out, counts


def main():
    p = argparse.ArgumentParser(
        description=(
            "Extract accession->chromosome mapping from an NCBI assembly report "
            "(preferred) or FASTA (default)."
        )
    )
    p.add_argument(
        "-a",
        "--assembly-report",
        help="NCBI assembly_report.txt (plain or .gz). If provided, overrides FASTA.",
    )
    p.add_argument(
        "-f",
        "--fasta",
        default="-",
        help="Genome FASTA (.gz OK). Default: '-' (stdin) if no assembly report.",
    )
    p.add_argument(
        "-s",
        "--style",
        choices=["raw", "ucsc"],
        default="ucsc",
        help="Naming style (default: ucsc).",
    )
    p.add_argument(
        "-k",
        "--keep",
        default=None,
        help=(
            "Comma-separated categories to KEEP "
            f"({', '.join(sorted(CATEGORY_CHOICES))})."
        ),
    )
    p.add_argument(
        "-r",
        "--remove",
        default=None,
        help=(
            "Comma-separated categories to REMOVE "
            f"({', '.join(sorted(CATEGORY_CHOICES))})."
        ),
    )
    p.add_argument(
        "-l",
        "--list-categories",
        action="store_true",
        help="Print category counts to stderr before writing TSV.",
    )
    p.add_argument(
        "-o",
        "--out-tsv",
        default=None,
        help="Write TSV (default: stdout).",
    )
    args = p.parse_args()

    keep_list = parse_category_list(args.keep)
    remove_list = parse_category_list(args.remove)

    if args.assembly_report:
        records, counts = build_map_from_assembly_report(
            args.assembly_report, args.style, keep_list, remove_list
        )
    else:
        records, counts = build_map_from_fasta(
            args.fasta, args.style, keep_list, remove_list
        )

    if args.list_categories:
        sys.stderr.write("Category counts (pre-filter):\n")
        for k in sorted(counts.keys()):
            sys.stderr.write(f"  {k:16s} {counts[k]:8d}\n")
        sys.stderr.write("\n")

    seen = set()
    lines = []
    for acc, name, _ in records:
        if (acc, name) in seen:
            continue
        seen.add((acc, name))
        lines.append(f"{acc}\t{name}")

    out = "\n".join(lines) + ("\n" if lines else "")
    if args.out_tsv:
        with open(args.out_tsv, "w", encoding="utf-8") as oh:
            oh.write(out)
    else:
        sys.stdout.write(out)


if __name__ == "__main__":
    main()
