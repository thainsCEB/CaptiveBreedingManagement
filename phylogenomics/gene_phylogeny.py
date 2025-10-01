#!/usr/bin/env python3
"""
Gene alignment & phylogeny pipeline (mitogenome version)
"""

import argparse
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def parse_args():
    p = argparse.ArgumentParser(
        description="Build per-gene alignments from mitogenome FASTAs, concatenate, and optionally run IQ-TREE",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--map", required=True, help="Two-column text file: sample_name <tab> path_to_mitogenome_fasta")
    p.add_argument("-o", "--out-dir", required=True, help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=4, help="Threads for MAFFT and IQ-TREE where applicable")
    p.add_argument(
        "-d", "--datatype",
        choices=["DNA", "AA"],
        default="DNA",
        help="Partition datatype written to partitions file",
    )
    p.add_argument(
        "-i", "--iqtree",
        action="store_true",
        help="Run IQ-TREE on the concatenated alignment with partitions",
    )
    p.add_argument(
        "-b", "--ufboot",
        type=int,
        default=None,
        help="Ultrafast bootstrap replicates for IQ-TREE (adds '-bb <N>'). Ignored unless --iqtree is set.",
    )
    p.add_argument(
        "-g", "--outgroup",
        default=None,
        help="Outgroup taxon name for IQ-TREE rooting (adds '-o <taxon>'). Ignored unless --iqtree is set.",
    )
    p.add_argument(
        "-p", "--prefix",
        default=None,
        help="Prefix for IQ-TREE output files (passed to '-pre'). Ignored unless --iqtree is set.",
    )
    p.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force overwrite of existing output files",
    )
    return p.parse_args()


def read_fasta(fp: Path) -> List[Tuple[str, str]]:
    records = []
    header = None
    seq_chunks: List[str] = []
    with fp.open() as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.replace(" ", ""))
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records


def write_fasta(records: List[Tuple[str, str]], fp: Path, linewrap: int = 80, force: bool = False):
    if fp.exists() and not force:
        sys.stderr.write(f"ERROR: Output file exists and --force not set: {fp}\n")
        sys.exit(1)
    with fp.open("w") as out:
        for h, s in records:
            out.write(f">{h}\n")
            if linewrap and linewrap > 0:
                for i in range(0, len(s), linewrap):
                    out.write(s[i : i + linewrap] + "\n")
            else:
                out.write(s + "\n")


def run_cmd(cmd: List[str], cwd: Optional[Path] = None):
    try:
        subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Command failed: {' '.join(cmd)}\n")
        raise e


def align_with_mafft(in_fp: Path, out_fp: Path, threads: int, force: bool):
    if out_fp.exists() and not force:
        sys.stderr.write(f"ERROR: Alignment file exists and --force not set: {out_fp}\n")
        sys.exit(1)
    cmd = ["mafft", "--auto", "--thread", str(threads), str(in_fp)]
    with out_fp.open("w") as out:
        try:
            subprocess.run(cmd, stdout=out, check=True)
        except FileNotFoundError:
            sys.stderr.write("ERROR: 'mafft' not found in PATH.\n")
            raise


def load_alignment(fp: Path) -> Dict[str, str]:
    return {h: s for h, s in read_fasta(fp)}


def concatenate_alignments(aln_files: List[Tuple[str, Path]]) -> Tuple[List[str], Dict[str, str], List[Tuple[str, int, int]]]:
    per_gene: List[Tuple[str, Dict[str, str]]] = []
    taxa = set()
    gene_lengths: List[Tuple[str, int]] = []

    for gene, fp in aln_files:
        aln = load_alignment(fp)
        per_gene.append((gene, aln))
        for t in aln.keys():
            taxa.add(t)
        lengths = {len(s) for s in aln.values()}
        if len(lengths) != 1:
            raise ValueError(f"Alignment {fp} has sequences of differing lengths: {lengths}")
        gene_len = lengths.pop() if lengths else 0
        gene_lengths.append((gene, gene_len))

    taxa_list = sorted(taxa)
    concat_map: Dict[str, List[str]] = {t: [] for t in taxa_list}
    partitions: List[Tuple[str, int, int]] = []

    pos = 1
    for (gene, aln), (_, gene_len) in zip(per_gene, gene_lengths):
        for t in taxa_list:
            seg = aln.get(t, "-" * gene_len)
            concat_map[t].append(seg)
        start = pos
        end = pos + gene_len - 1 if gene_len > 0 else pos - 1
        partitions.append((gene, start, end))
        pos = end + 1

    concat_str_map: Dict[str, str] = {t: "".join(parts) for t, parts in concat_map.items()}
    return taxa_list, concat_str_map, partitions


def write_partitions(partitions: List[Tuple[str, int, int]], datatype: str, fp: Path, force: bool = False):
    if fp.exists() and not force:
        sys.stderr.write(f"ERROR: Partition file exists and --force not set: {fp}\n")
        sys.exit(1)
    with fp.open("w") as out:
        for gene, start, end in partitions:
            out.write(f"{datatype}, {gene} = {start}-{end}\n")


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_gene_dir = out_dir / "01.per_gene"
    aligned_dir = out_dir / "02.aligned"
    per_gene_dir.mkdir(exist_ok=True)
    aligned_dir.mkdir(exist_ok=True)

    if (args.ufboot is not None or args.outgroup is not None or args.prefix is not None) and not args.iqtree:
        sys.stderr.write("Note: --ufboot/--outgroup/--prefix are ignored unless --iqtree is specified.\n")

    # Parse map file
    sample_to_fasta: List[Tuple[str, Path]] = []
    with open(args.map) as f:
        for line in f:
            if not line.strip():
                continue
            toks = line.strip().split("\t")
            if len(toks) < 2:
                sys.stderr.write(f"Invalid line in map file: {line}")
                continue
            sample, fasta_path = toks[0], Path(toks[1])
            sample_to_fasta.append((sample, fasta_path))

    # Collect per-gene sequences
    gene_to_records: Dict[str, List[Tuple[str, str]]] = {}
    for sample, fasta_fp in sample_to_fasta:
        records = read_fasta(fasta_fp)
        for hdr, seq in records:
            gene = hdr.split()[0]
            gene_to_records.setdefault(gene, []).append((sample, seq))

    # Write per-gene FASTAs
    gene_files: List[Tuple[str, Path]] = []
    for gene, recs in gene_to_records.items():
        out_fp = per_gene_dir / f"{gene}.fasta"
        write_fasta(recs, out_fp, force=args.force)
        gene_files.append((gene, out_fp))

    if not gene_files:
        sys.stderr.write("No genes found in input FASTAs.\n")
        sys.exit(1)

    # Align with MAFFT
    aln_files: List[Tuple[str, Path]] = []
    for gene, fp in gene_files:
        aln_fp = aligned_dir / f"{gene}.aln.fasta"
        align_with_mafft(fp, aln_fp, threads=args.threads, force=args.force)
        aln_files.append((gene, aln_fp))

    # Concatenate
    taxa_list, concat_map, partitions = concatenate_alignments(aln_files)
    concat_fp = out_dir / "concatenated.fasta"
    partition_fp = out_dir / "partitions.txt"

    concat_records = [(t, seq) for t, seq in concat_map.items()]
    write_fasta(concat_records, concat_fp, force=args.force)
    write_partitions(partitions, args.datatype, partition_fp, force=args.force)

    # Run IQ-TREE if requested
    if args.iqtree:
        iqtree_exe = shutil.which("iqtree")
        if iqtree_exe is None:
            sys.stderr.write("ERROR: 'iqtree' not found in PATH. Install IQ-TREE first.\n")
            sys.exit(2)

        cmd = ["iqtree", "-s", str(concat_fp), "-p", str(partition_fp), "-T", str(args.threads), "-m", "MFP"]
        if args.ufboot is not None:
            cmd += ["-bb", str(args.ufboot)]
        if args.outgroup is not None:
            cmd += ["-o", args.outgroup]
        if args.prefix is not None:
            cmd += ["-pre", str(Path(out_dir) / args.prefix)]
        run_cmd(cmd)

    print("Done.")
    print(f"Per-gene FASTAs: {per_gene_dir}")
    print(f"Aligned FASTAs: {aligned_dir}")
    print(f"Concatenated alignment: {concat_fp}")
    print(f"Partition file: {partition_fp}")
    if args.iqtree:
        print("IQ-TREE run complete. Check output files next to concatenated.fasta (or under prefix if specified).")


if __name__ == "__main__":
    main()
