#!/usr/bin/env python3
# Name: mitoForger.py
# Title: MitoForger — batch subsampling, MitoZ runs, and mitogenome post-processing
# Description:
#   - Reads a sample sheet (TSV) with 4 columns: sampleID, read1, read2, speciesName
#     * speciesName MUST be the scientific name (e.g., "Ara glaucogularis" or "Homo sapiens"),
#       and is passed verbatim to MitoZ as: --species 'speciesName'
#   - Optionally subsamples paired reads to a given percent using seqfu (interleave/deinterleave) + seqtk
#   - Runs MitoZ for each sample (native or via Singularity) with sensible defaults
#   - Post-processes MitoZ outputs:
#       * Extract protein-coding (PCG) and 'all genes' FASTA with simplified headers
#       * Rewrites complete mitogenome FASTA header to '>MT mitochondrion'
#       * Copies the .gbf file next to derived FASTAs
#       * Appends 'Potential missing genes:' lines to a global mito.log
#   - Optionally tar.gz each <sample>.result directory
#  Coverage & MDR notes:
#   - MitoZ generally performs well when the mitochondrial read proportion (MDR)
#     is ≥0.1% of total reads, ideally ≥0.3–0.5%.
#   - If MDR <0.05% or mitochondrial coverage <30–50×, assemblies may be fragmented
#     or missing genes (especially tRNAs).
#   - Avoid excessive subsampling (-p) when coverage is uncertain; ensure enough
#     mitochondrial reads remain for reliable assembly.
#   - For very low-MDR datasets, consider read filtering or reference-guided
#     assembly (e.g., NOVOPlasty, GetOrganelle) as alternatives.
#
# Author: Taylor Hains
# Date: 2025-10-23
#
# Quick start:
#   python mitoForger.py -s samples.tsv -o out -t 12 -p 20 -T
#   python mitoForger.py -s samples.tsv -o out -t 12 -S /path/mitoz.sif
#

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Optional, Tuple, List

def sh(cmd: List[str], cwd: Optional[str] = None, capture: bool = False) -> Tuple[int, str]:
    try:
        if capture:
            out = subprocess.check_output(cmd, cwd=cwd, stderr=subprocess.STDOUT, text=True)
            return 0, out
        else:
            subprocess.check_call(cmd, cwd=cwd)
            return 0, ""
    except subprocess.CalledProcessError as e:
        msg = f"Command failed ({e.returncode}): {' '.join(cmd)}\n---\n{e.output if hasattr(e, 'output') and e.output else ''}"
        raise RuntimeError(msg) from e

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def is_fraction_or_percent(x: float) -> float:
    if x <= 0:
        return 0.0
    if x > 1.0:
        return min(x / 100.0, 1.0)
    return x

def gzip_file(src: Path, dst: Optional[Path] = None):
    if dst is None:
        dst = src.with_suffix(src.suffix + ".gz")
    with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    src.unlink(missing_ok=True)
    return dst

def find_one(pattern: str, root: Path) -> Optional[Path]:
    hits = list(root.rglob(pattern))
    return hits[0] if hits else None

def simplify_gene_headers(in_fa: Path, out_fa: Path):
    with open(in_fa, "r") as fi, open(out_fa, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                m = re.match(r"^>[^;]*;([^;()]+)", line.strip())
                if m:
                    fo.write(f">{m.group(1).strip()}\n")
                else:
                    fo.write(line)
            else:
                fo.write(line)

def rewrite_mitogenome_header(in_fa: Path, out_fa: Path):
    with open(in_fa, "r") as fi, open(out_fa, "w") as fo:
        wrote_header = False
        for line in fi:
            if line.startswith(">"):
                if not wrote_header:
                    fo.write(">MT mitochondrion\n")
                    wrote_header = True
            else:
                fo.write(line)

def tar_gz_dir(src_dir: Path, tar_path: Path):
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(src_dir, arcname=src_dir.name)

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def subsample_pair(read1: Path, read2: Path, outprefix: Path, frac: float):
    if frac <= 0 or frac >= 1.0:
        raise ValueError("subsample_pair() requires 0 < frac < 1.0")
    if not which("seqtk"):
        raise RuntimeError("seqtk not found in PATH.")
    if not which("seqfu"):
        raise RuntimeError("seqfu not found in PATH (needed for interleave/deinterleave).")

    tmpdir = outprefix.parent / f".tmp_{outprefix.name}"
    ensure_dir(tmpdir)
    interleaved = tmpdir / "interleaved.fastq"
    interleaved_sub = tmpdir / "interleaved.sub.fastq"

    # Interleave
    sh(["seqfu", "ilv", "-1", str(read1), "-2", str(read2), "-o", str(interleaved)])

    # Subsample with seqtk
    with open(interleaved, "rb") as fi:
        p = subprocess.Popen(["seqtk", "sample", "-s", "100", "-", str(frac)], stdin=fi, stdout=subprocess.PIPE)
        with open(interleaved_sub, "wb") as fo:
            shutil.copyfileobj(p.stdout, fo)
        p.wait()
        if p.returncode != 0:
            raise RuntimeError("seqtk sample failed.")

    # Deinterleave
    deint_dir = tmpdir / "deint"
    ensure_dir(deint_dir)
    sh(["seqfu", "deinterleave", "-o", str(deint_dir), str(interleaved_sub)])

    r1_candidates = list(deint_dir.glob("*_1.fastq"))
    r2_candidates = list(deint_dir.glob("*_2.fastq"))
    if not r1_candidates or not r2_candidates:
        raise RuntimeError("Failed to locate deinterleaved files.")

    r1_out = outprefix.parent / f"{outprefix.name}_1.sub.fq"
    r2_out = outprefix.parent / f"{outprefix.name}_2.sub.fq"
    shutil.move(str(r1_candidates[0]), r1_out)
    shutil.move(str(r2_candidates[0]), r2_out)

    r1_gz = gzip_file(r1_out)
    r2_gz = gzip_file(r2_out)
    shutil.rmtree(tmpdir, ignore_errors=True)
    return r1_gz, r2_gz

def build_mitoz_cmd(args, sample_id: str, fq1: Path, fq2: Path, sample_outdir: Path) -> List[str]:
    cmd = []
    if args.singularity_image:
        if not which("singularity"):
            raise RuntimeError("singularity not found in PATH, but --singularity-image was provided.")
        bind_paths = {str(args.outdir.resolve())}
        bind_paths.update({str(Path(fq1).parent.resolve()), str(Path(fq2).parent.resolve())})
        binds = ",".join(sorted(bind_paths))
        cmd = ["singularity", "exec", "--bind", binds, str(args.singularity_image)]

    mitoz_bin = args.mitoz_cmd
    cmd += [
        mitoz_bin, "all",
        "--genetic_code", str(args.genetic_code),
        "--clade", args.clade,
        "--outprefix", str(sample_outdir / sample_id),
        "--workdir", str(sample_outdir / f"temp_{sample_id}"),
        "--thread_number", str(args.threads),
        "--fq1", str(fq1),
        "--fq2", str(fq2),
        "--fastq_read_length", str(args.read_length),
        "--insert_size", str(args.insert_size),
        "--requiring_taxa", args.requiring_taxa,
        "--species", args.species
    ]
    return cmd

def run_mitoz_for_sample(args, sample_id: str, species_name: Optional[str], fq1: Path, fq2: Path, sample_outdir: Path, logfh):
    species_arg = species_name if species_name else args.species  # keep spaces verbatim
    cmd = build_mitoz_cmd(args, sample_id, fq1, fq2, sample_outdir)
    try:
        sp_idx = cmd.index("--species")
        cmd[sp_idx + 1] = species_arg
    except ValueError:
        pass

    logfh.write(f"[STEP] Running mitoz for {sample_id}\n")
    logfh.write(f"[CMD] {' '.join(cmd)}\n")
    sys.stdout.flush()
    sh(cmd)

def postprocess_mitoz_outputs(sample_id: str, outdir: Path, logfh, tar_results: bool):
    root = outdir / "mitogenomes"
    sample_root = root / f"{sample_id}.result"
    if not sample_root.exists():
        logfh.write(f"[WARN] {sample_id}: result dir not found: {sample_root}\n")
        return

    pcg_fa = find_one("*.gbf.cds.fasta", sample_root)
    all_fa = find_one("*.gbf.gene.fasta", sample_root)
    mt_fa  = find_one("*.gbf.fasta", sample_root)

    out_pcg = root / f"{sample_id}.mito.genes.pcg.fasta"
    out_all = root / f"{sample_id}.mito.genes.all.fasta"
    out_mt  = root / f"{sample_id}.mitogenome.fasta"

    if pcg_fa:
        simplify_gene_headers(pcg_fa, out_pcg)
    else:
        logfh.write(f"[WARN] {sample_id}: PCG fasta not found.\n")

    if all_fa:
        simplify_gene_headers(all_fa, out_all)
    else:
        logfh.write(f"[WARN] {sample_id}: ALL-genes fasta not found.\n")

    if mt_fa:
        rewrite_mitogenome_header(mt_fa, out_mt)
    else:
        logfh.write(f"[WARN] {sample_id}: Mitogenome fasta not found.\n")

    gbf = find_one("*.gbf", sample_root)
    if gbf:
        shutil.copy2(gbf, root / f"{sample_id}.gbf")
    else:
        logfh.write(f"[WARN] {sample_id}: .gbf file not found.\n")

    summary = find_one("summary.txt", sample_root)
    if summary:
        with open(summary, "r") as sf, open(outdir / "mito.log", "a") as mlf:
            for line in sf:
                if "Potential missing genes:" in line:
                    mlf.write(f"{sample_id}\t{line.strip()}\n")

    if tar_results and sample_root.exists():
        tar_path = sample_root.with_suffix(".tar.gz")
        tar_gz_dir(sample_root, tar_path)

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="MitoForger — batch subsampling, MitoZ runs, and mitogenome post-processing."
    )
    # Core I/O
    p.add_argument("-t", "--threads", type=int, default=8, help="Threads for MitoZ.")
    p.add_argument("-s", "--sample-sheet", required=True, help="TSV: sampleID  read1  read2  speciesName (scientific).")
    p.add_argument("-p", "--percent", type=float, default=100.0, help="Percent (or fraction) of reads to keep. 100 or 1.0 skips subsampling.")
    p.add_argument("-o", "--outdir", required=True, help="Output directory (creates subdirs).")

    # MitoZ config
    p.add_argument("-m", "--mitoz-cmd", default="mitoz", help="Path to 'mitoz' executable.")
    p.add_argument("-S", "--singularity-image", default=None, help="Path to MitoZ Singularity (.sif); auto-binds OUTDIR and FASTQ parent dirs.")
    p.add_argument("-g", "--genetic-code", type=int, default=2, help="MitoZ --genetic_code.")
    p.add_argument("-c", "--clade", default="Chordata", help="MitoZ --clade.")
    p.add_argument("-r", "--requiring-taxa", default="Aves", help="MitoZ --requiring_taxa.")
    p.add_argument("-a", "--species", default="Unknown_species", help="Fallback species name if sheet omits it.")
    p.add_argument("-l", "--read-length", type=int, default=150, help="MitoZ --fastq_read_length.")
    p.add_argument("-i", "--insert-size", type=int, default=200, help="MitoZ --insert_size.")

    # Behaviors
    p.add_argument("-T", "--tar-results", action="store_true", help="Tar.gz each <sample>.result directory after post-processing.")
    return p.parse_args()

def main():
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    mito_dir = outdir / "mitogenomes"
    mito_fastq = outdir / "mito_fastq"
    ensure_dir(mito_dir)
    ensure_dir(mito_fastq)

    frac = is_fraction_or_percent(args.percent)
    do_subsample = 0 < frac < 1.0

    if do_subsample:
        if not which("seqtk"):
            sys.exit("ERROR: seqtk not found in PATH (required for subsampling).")
        if not which("seqfu"):
            sys.exit("ERROR: seqfu not found in PATH (required for interleave/deinterleave).")
    if args.singularity_image:
        if not which("singularity"):
            sys.exit("ERROR: singularity not found in PATH but --singularity-image provided.")
    else:
        if not which(args.mitoz_cmd):
            sys.exit(f"ERROR: '{args.mitoz_cmd}' not found in PATH.")

    with open(outdir / "batch_mitoz.log", "a") as logfh:
        with open(args.sample_sheet) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or row[0].strip().startswith("#"):
                    continue
                if len(row) < 3:
                    logfh.write(f"[WARN] Skipping malformed line: {row}\n")
                    continue
                sample_id, r1, r2 = row[0].strip(), row[1].strip(), row[2].strip()
                species_name = row[3].strip() if len(row) >= 4 and row[3].strip() else None  # keep spaces

                logfh.write(f"[STEP] Sample {sample_id}\n")
                if do_subsample:
                    logfh.write(f"[STEP] Subsampling {sample_id} to {frac:.4f} of reads\n")
                    outprefix = mito_fastq / sample_id
                    fq1_gz, fq2_gz = subsample_pair(Path(r1), Path(r2), outprefix, frac)
                else:
                    fq1_gz = Path(r1)
                    fq2_gz = Path(r2)

                run_mitoz_for_sample(args, sample_id, species_name, fq1_gz, fq2_gz, mito_dir, logfh)
                postprocess_mitoz_outputs(sample_id, outdir, logfh, args.tar_results)

    print("All done. See batch_mitoz.log and mito.log for details.")

if __name__ == "__main__":
    main()
