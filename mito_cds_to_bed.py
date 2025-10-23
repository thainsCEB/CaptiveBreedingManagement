#!/usr/bin/env python3
"""
mito_cds_to_bed.py
Extract mitochondrial protein-coding CDS to BED (BED12 by default) with the
gene name in the BED "name" field, ready for:
    bedtools getfasta -fi genome.fa -bed mito_13pcg.bed -s -name -split > mito_13pcg.fa

Defaults
--------
- Restricts to the canonical 13 mitochondrial protein-coding genes:
  ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6
- Accepts common synonyms (COI/COX I, COB, NAD1, mt-ATP6, etc.)
- Auto-detects MT contigs (chrM, MT, mitochondrion...), or set with --mito-ids
- BED12 so bedtools -split reconstructs spliced CDS
- NEW: --sort / -S sorts BED by (chrom, start, end, name) before writing

Override:
- Use --include-all-cds to output *all* mitochondrial CDS, not just the 13.
- Or use --genes "COX1,COX2,COX3" to specify a custom set.

Usage
-----
# Default: 13 PCGs only, sorted BED
python mito_cds_to_bed.py -a annotations.gff3 -o mito_13pcg.bed --sort

# All MT CDS unsorted
python mito_cds_to_bed.py -a annotations.gtf --include-all-cds -o mito_all_cds.bed

# Custom list + sorted BED6
python mito_cds_to_bed.py -a annotations.gff3 --genes "COX1,COX2,COX3,ND1" --bed6 -S -o mito_subset.bed

Then extract sequences:
bedtools getfasta -fi genome.fa -bed mito_13pcg.bed -s -name -split > mito_13pcg.fa
"""

import argparse
import sys
import re
from collections import defaultdict, namedtuple
from pathlib import Path
from typing import Optional, List, Tuple

Feat = namedtuple("Feat", "chrom start end strand")

# Canonical 13 PCGs (uppercase)
CANONICAL_13 = [
    "ATP6", "ATP8",
    "COX1", "COX2", "COX3",
    "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
]

def _syn_re(pat):
    return re.compile(pat, re.IGNORECASE)

SYNONYMS = [
    (_syn_re(r"^(?:mt)?atp6$"), "ATP6"),
    (_syn_re(r"^(?:mt)?atp8$"), "ATP8"),
    (_syn_re(r"^(?:mt)?(?:co1|cox1|cox[i1]|cytochromecoxidasesubuniti|coi)$"), "COX1"),
    (_syn_re(r"^(?:mt)?(?:co2|cox2|coxii|cytochromecoxidasesubunitii|coii)$"), "COX2"),
    (_syn_re(r"^(?:mt)?(?:co3|cox3|coxiii|cytochromecoxidasesubunitiii|coiii)$"), "COX3"),
    (_syn_re(r"^(?:mt)?(?:cytb|cob|cytochromeb)$"), "CYTB"),
    (_syn_re(r"^(?:mt)?(?:nd1|nad1|nadhdehydrogenasesubunit1)$"), "ND1"),
    (_syn_re(r"^(?:mt)?(?:nd2|nad2|nadhdehydrogenasesubunit2)$"), "ND2"),
    (_syn_re(r"^(?:mt)?(?:nd3|nad3|nadhdehydrogenasesubunit3)$"), "ND3"),
    (_syn_re(r"^(?:mt)?(?:nd4|nad4|nadhdehydrogenasesubunit4)$"), "ND4"),
    (_syn_re(r"^(?:mt)?(?:nd4l|nad4l|nadhdehydrogenasesubunit4l)$"), "ND4L"),
    (_syn_re(r"^(?:mt)?(?:nd5|nad5|nadhdehydrogenasesubunit5)$"), "ND5"),
    (_syn_re(r"^(?:mt)?(?:nd6|nad6|nadhdehydrogenasesubunit6)$"), "ND6"),
]

def guess_format(attrs_field: str) -> str:
    if "=" in attrs_field and ";" in attrs_field and '"' not in attrs_field:
        return "gff3"
    if '="' in attrs_field or '"' in attrs_field:
        return "gtf"
    return "gff3"

def parse_gff3_attrs(s: str) -> dict:
    d = {}
    for part in s.strip().strip(";").split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = v.strip()
        else:
            bits = part.strip().split()
            if len(bits) == 2:
                d[bits[0]] = bits[1]
    return d

def parse_gtf_attrs(s: str) -> dict:
    d = {}
    for part in s.strip().strip(";").split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        k, v = part.split(" ", 1)
        v = v.strip().strip('"')
        d[k.strip()] = v
    return d

def sanitize_name(name: str) -> str:
    if not name:
        return "unknown"
    name = re.sub(r"\s+", "_", name)
    return name.replace("/", "_").replace("\\", "_")

def normalize_gene_string(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[\s\-\._]", "", s)
    s = s.replace("(", "").replace(")", "")
    return s

def canonicalize_gene(name: str) -> Optional[str]:
    """Return canonical 13-PCG symbol if name matches a synonym; else None."""
    if not name:
        return None
    n = normalize_gene_string(name)
    quick = {
        "atp6": "ATP6", "atp8": "ATP8",
        "co1": "COX1", "cox1": "COX1", "coi": "COX1", "coxi": "COX1",
        "co2": "COX2", "cox2": "COX2", "coii": "COX2",
        "co3": "COX3", "cox3": "COX3", "coiii": "COX3",
        "cytb": "CYTB", "cob": "CYTB", "cytochromeb": "CYTB",
        "nd1": "ND1", "nad1": "ND1",
        "nd2": "ND2", "nad2": "ND2",
        "nd3": "ND3", "nad3": "ND3",
        "nd4": "ND4", "nad4": "ND4",
        "nd4l": "ND4L", "nad4l": "ND4L",
        "nd5": "ND5", "nad5": "ND5",
        "nd6": "ND6", "nad6": "ND6",
    }
    if n.startswith("mt"):
        n2 = n[2:]
        if n2 in quick:
            return quick[n2]
    if n in quick:
        return quick[n]
    for rx, canon in SYNONYMS:
        if rx.match(n):
            return canon
    return None

def autodetect_mito_contigs(seen_chroms):
    mito = set()
    patt_any = re.compile(r"(mito|mitochond)", re.IGNORECASE)
    for c in seen_chroms:
        if c in ("chrM", "chrMT", "MT", "MtDNA"):
            mito.add(c)
        elif patt_any.search(c):
            mito.add(c)
    if not mito and "M" in seen_chroms:
        mito.add("M")
    return mito

def main():
    ap = argparse.ArgumentParser(description="Extract mitochondrial protein-coding CDS (default: canonical 13 PCGs) to BED.")
    ap.add_argument("-a", "--annotation", required=True, help="GFF3 or GTF annotation file")
    ap.add_argument("-o", "--out-bed", default="mito_13pcg.bed", help="Output BED path (BED12 by default)")
    ap.add_argument("--summary-tsv", default="mito_13pcg.summary.tsv",
                    help="Optional summary TSV (gene,chrom,start,end,strand,n_blocks,length)")
    ap.add_argument("-m", "--mito-ids", default=None,
                    help="Comma-separated contig IDs representing mitochondrion (e.g., chrM,MT); auto-detected if omitted.")
    ap.add_argument("--gene-tags", default="gene_name,gene,Name,locus_tag,product",
                    help="Comma-separated attribute keys to try for gene names (priority order).")
    ap.add_argument("--bed6", action="store_true", help="Write BED6 (merged per gene) instead of BED12.")
    ap.add_argument("--force", action="store_true", help="Overwrite output files if they exist.")
    ap.add_argument("--include-all-cds", action="store_true",
                    help="Include *all* MT CDS (not restricted to canonical 13).")
    ap.add_argument("--genes", default=None,
                    help="Custom gene list (comma-separated). If set, restricts to these (case-insensitive, canonicalized if possible).")
    ap.add_argument("--sort", "-S", action="store_true",
                    help="Sort BED by (chrom, start, end, name) before writing.")
    args = ap.parse_args()

    ann_path = Path(args.annotation)
    if not ann_path.exists():
        sys.exit(f"[ERROR] Annotation not found: {ann_path}")
    if Path(args.out_bed).exists() and not args.force:
        sys.exit(f"[ERROR] {args.out_bed} exists. Use --force to overwrite.")

    # Build allowed set
    if args.include_all_cds:
        allowed_canonical = None  # allow everything
    elif args.genes:
        allowed_canonical = set()
        for g in [x.strip() for x in args.genes.split(",") if x.strip()]:
            can = canonicalize_gene(g) or g.strip().upper()
            allowed_canonical.add(can)
    else:
        allowed_canonical = set(CANONICAL_13)

    # First pass: detect format, collect mappings
    id_to_gene = {}
    tx_to_gene = {}
    seen_chroms = set()
    fmt = None
    Raw = namedtuple("Raw", "chrom ftype start end strand attrs_line")
    raw_rows = []

    with ann_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            seen_chroms.add(chrom)
            if fmt is None:
                fmt = guess_format(attrs)
            if fmt == "gff3":
                A = parse_gff3_attrs(attrs)
                fid = A.get("ID")
                if ftype.lower() == "gene" and fid:
                    name = None
                    for k in args.gene_tags.split(","):
                        if k in A and A[k]:
                            name = A[k]
                            break
                    if not name:
                        name = A.get("Name", fid)
                    id_to_gene[fid] = sanitize_name(name)
                elif ftype.lower() in ("mrna", "transcript"):
                    tid = A.get("ID")
                    name = None
                    for k in args.gene_tags.split(","):
                        if k in A and A[k]:
                            name = A[k]
                            break
                    if not name:
                        parent_gene = A.get("Parent")
                        if parent_gene and parent_gene in id_to_gene:
                            name = id_to_gene[parent_gene]
                    if tid and name:
                        tx_to_gene[tid] = sanitize_name(name)
            else:
                A = parse_gtf_attrs(attrs)
                if ftype.lower() == "gene":
                    gid = A.get("gene_id")
                    name = None
                    for k in args.gene_tags.split(","):
                        kk = k if k != "gene" else "gene_name"
                        name = A.get(kk)
                        if name:
                            break
                    if not name:
                        name = gid
                    if gid and name:
                        id_to_gene[gid] = sanitize_name(name)
                elif ftype.lower() in ("transcript",):
                    tid = A.get("transcript_id")
                    name = A.get("gene_name") or A.get("gene_id")
                    if tid and name:
                        tx_to_gene[tid] = sanitize_name(name)

            raw_rows.append(Raw(chrom, ftype, int(start), int(end), strand, attrs))

    if not fmt:
        sys.exit("[ERROR] Could not detect file format or found no features.")

    # Mito contigs
    if args.mito_ids:
        mito_contigs = set([x.strip() for x in args.mito_ids.split(",") if x.strip()])
    else:
        mito_contigs = autodetect_mito_contigs(seen_chroms)
        if not mito_contigs:
            print("[WARN] Could not auto-detect mitochondrial contig(s). "
                  "Provide via --mito-ids (e.g., chrM,MT). Proceeding with none.",
                  file=sys.stderr)

    # Collect CDS per gene on MT contigs (filtering to 13 by default)
    gene_blocks = defaultdict(list)
    gene_meta = {}
    missing_name_count = 0
    kept_counts = defaultdict(int)
    skipped_noncanonical = 0

    def pick_gene_name(attrs_line: str) -> str:
        nonlocal missing_name_count
        if fmt == "gff3":
            A = parse_gff3_attrs(attrs_line)
            for k in args.gene_tags.split(","):
                if k in A and A[k]:
                    return sanitize_name(A[k])
            parent = A.get("Parent")
            if parent:
                pid = parent.split(",")[0]
                if pid in tx_to_gene:
                    return sanitize_name(tx_to_gene[pid])
            name = A.get("Name") or A.get("gene") or A.get("gene_name") or f"unknown_{missing_name_count+1}"
            if name.startswith("unknown_"):
                missing_name_count += 1
            return sanitize_name(name)
        else:
            A = parse_gtf_attrs(attrs_line)
            name = A.get("gene_name") or A.get("gene_id")
            if not name:
                tid = A.get("transcript_id")
                name = tx_to_gene.get(tid, f"unknown_{missing_name_count+1}")
            if str(name).startswith("unknown_"):
                missing_name_count += 1
            return sanitize_name(name)

    for row in raw_rows:
        if row.ftype != "CDS":
            continue
        if mito_contigs and row.chrom not in mito_contigs:
            continue
        raw_name = pick_gene_name(row.attrs_line)
        if allowed_canonical is None:
            canon = canonicalize_gene(raw_name) or raw_name.upper()
            keep = True
        else:
            canon = canonicalize_gene(raw_name)
            keep = (canon in allowed_canonical)
            if not keep:
                skipped_noncanonical += 1
                continue
        gname = canon
        gene_blocks[gname].append(Feat(row.chrom, row.start, row.end, row.strand))
        kept_counts[gname] += 1
        if gname not in gene_meta:
            gene_meta[gname] = (row.chrom, row.strand)
        else:
            c0, s0 = gene_meta[gname]
            if c0 != row.chrom or s0 != row.strand:
                print(f"[WARN] Gene {gname} has mixed chrom/strand; keeping first location.", file=sys.stderr)

    if not gene_blocks:
        msg = "[ERROR] No mitochondrial CDS matched"
        if allowed_canonical is None:
            msg += " (searched for all MT CDS)."
        elif args.genes:
            msg += f" the provided set: {sorted(allowed_canonical)}."
        else:
            msg += " the canonical 13 PCGs."
        msg += " Try --mito-ids to specify the MT contig(s)."
        sys.exit(msg)

    # Assemble BED rows in memory
    bed_rows: List[Tuple] = []
    for gname, blocks in gene_blocks.items():
        blocks = sorted(blocks, key=lambda x: (x.start, x.end))
        chrom, strand = gene_meta[gname]
        tx_start0 = min(b.start for b in blocks) - 1
        tx_end = max(b.end for b in blocks)
        if args.bed6:
            bed_rows.append((chrom, tx_start0, tx_end, gname, 0, strand))
        else:
            sizes = []
            rel_starts = []
            for b in blocks:
                b_start0 = b.start - 1
                size = b.end - b.start + 1
                sizes.append(size)
                rel_starts.append(b_start0 - tx_start0)
            sizes_str = ",".join(str(x) for x in sizes) + ","
            rel_starts_str = ",".join(str(x) for x in rel_starts) + ","
            block_count = len(sizes)
            # Keep all 12 BED fields as a tuple for sorting/writing
            bed_rows.append((
                chrom, tx_start0, tx_end, gname, 0, strand,
                tx_start0, tx_end, "0", block_count, sizes_str, rel_starts_str
            ))

    # Sort if requested (chrom, start, end, name)
    if args.sort:
        if args.bed6:
            bed_rows.sort(key=lambda r: (str(r[0]), int(r[1]), int(r[2]), str(r[3])))
        else:
            bed_rows.sort(key=lambda r: (str(r[0]), int(r[1]), int(r[2]), str(r[3])))

    # Write BED
    out_bed = Path(args.out_bed)
    with out_bed.open("w") as out:
        if args.bed6:
            for r in bed_rows:
                out.write("\t".join(map(str, r)) + "\n")
        else:
            for r in bed_rows:
                out.write("\t".join(map(str, r)) + "\n")

    # Summary
    if args.summary_tsv:
        with Path(args.summary_tsv).open("w") as sf:
            sf.write("gene\tchrom\tstart0\tend\tstrand\tn_blocks\tlength_bp\tcds_parts_seen\n")
            for gname, blocks in sorted(gene_blocks.items()):
                blocks = sorted(blocks, key=lambda x: (x.start, x.end))
                chrom, strand = gene_meta[gname]
                tx_start0 = min(b.start for b in blocks) - 1
                tx_end = max(b.end for b in blocks)
                n_blocks = len(blocks)
                length_bp = sum(b.end - b.start + 1 for b in blocks)
                sf.write(f"{gname}\t{chrom}\t{tx_start0}\t{tx_end}\t{strand}\t{n_blocks}\t{length_bp}\t{kept_counts[gname]}\n")

    # Messages
    if allowed_canonical is None:
        scope_msg = "all MT CDS"
    elif args.genes:
        scope_msg = f"{len(allowed_canonical)} provided gene(s)"
    else:
        scope_msg = "the canonical 13 PCGs"
    print(f"[OK] Wrote {len(bed_rows)} genes ({scope_msg}) to {out_bed}")
    if skipped_noncanonical and allowed_canonical is not None:
        print(f"[INFO] Skipped {skipped_noncanonical} non-canonical MT CDS not in target set.", file=sys.stderr)
    if args.summary_tsv:
        print(f"[OK] Wrote summary to {args.summary_tsv}")
    if not args.bed6:
        print(f"Tip: bedtools getfasta -fi genome.fa -bed {out_bed} -s -name -split > mito_pcgs.fa")
    else:
        print(f"Tip: bedtools getfasta -fi genome.fa -bed {out_bed} -s -name > mito_pcgs.fa")

if __name__ == "__main__":
    main()
