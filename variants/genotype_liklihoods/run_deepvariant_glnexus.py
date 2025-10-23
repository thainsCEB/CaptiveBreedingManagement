#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeepVariant + GLnexus batch runner
----------------------------------
Loops over samples to run DeepVariant per-BAM (producing both VCF & GVCF), then
cohort-merges all GVCFs with GLnexus, converts BCF -> VCF.GZ, and indexes with tabix.

- DeepVariant via Singularity or Docker (GPU optional).
- GLnexus via container or host (`glnexus_cli` in PATH).
- Cleans DeepVariant intermediates by default (toggle with -K).
- Supports regions/PAR BEDs, single region, stats report, and pass-through flags.
- Adds bcftools/tabix conversion and indexing for cohort output.
- NEW: `--output-prefix/-O` to control per-sample filenames **and** cohort filenames.

Author: Taylor Hains
Date: 2025-10-18
"""

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Optional

# --------------------------- helpers ---------------------------

def log(msg: str):
    print(f"[dv-glx] {msg}", flush=True)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def read_samples(path: Path) -> List[Tuple[str, Path]]:
    rows = []
    with path.open() as f:
        sniffer = csv.Sniffer()
        sample = f.read(1024)
        f.seek(0)
        delimiter = "\t"
        try:
            dialect = sniffer.sniff(sample)
            delimiter = dialect.delimiter
        except Exception:
            delimiter = "\t"
        reader = csv.reader(f, delimiter=delimiter)
        # skip header if present
        first = True
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if first and any(h.lower() in ("sample","sample_id","id") for h in row[:2]):
                first = False
                continue
            first = False
            if len(row) < 2:
                raise ValueError(f"Bad row (need 2 cols): {row}")
            sid = row[0].strip()
            bam = Path(row[1].strip())
            rows.append((sid, bam))
    return rows

def path_binds_for(paths: List[Optional[Path]]) -> List[Path]:
    binds = set()
    for p in paths:
        if p is None:
            continue
        p = p.resolve()
        if p.is_file():
            binds.add(p.parent)
        else:
            binds.add(p)
    return sorted(binds)

def run_checked(cmd: List[str], log_file: Path) -> None:
    ensure_dir(log_file.parent)
    with log_file.open("w") as lf:
        lf.write("$ " + " ".join(shlex.quote(c) for c in cmd) + "\n")
        lf.flush()
        proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(f"Command failed (exit {ret}). See log: {log_file}")

def build_singularity_cmd(image: str, binds: List[Path], inner_cmd: List[str], gpu: bool) -> List[str]:
    cmd = ["singularity", "exec"]
    if gpu:
        cmd.append("--nv")
    for b in binds:
        cmd += ["-B", f"{str(b)}:{str(b)}"]
    cmd += [image]
    cmd += inner_cmd
    return cmd

def build_docker_cmd(image: str, binds: List[Path], inner_cmd: List[str], gpu: bool, workdir: Optional[Path]=None) -> List[str]:
    cmd = ["docker", "run", "--rm"]
    if gpu:
        cmd += ["--gpus", "all"]
    for b in binds:
        cmd += ["-v", f"{str(b)}:{str(b)}"]
    if workdir is not None:
        cmd += ["-w", str(workdir)]
    cmd += [image]
    cmd += inner_cmd
    return cmd

def build_deepvariant_inner(
    model_type: str,
    reference: Path,
    bam: Path,
    vcf_out: Path,
    gvcf_out: Path,
    shards: int,
    par_bed: Optional[Path],
    regions_bed: Optional[Path],
    region_str: Optional[str],
    vcf_stats_report: bool,
    intermediate_dir: Optional[Path],
    extra: List[str]
) -> List[str]:
    inner = [
        "/opt/deepvariant/bin/run_deepvariant",
        f"--model_type={model_type}",
        f"--ref={reference}",
        f"--reads={bam}",
        f"--output_vcf={vcf_out}",
        f"--output_gvcf={gvcf_out}",
        f"--num_shards={shards}",
        "--logging_dir=" + str((vcf_out.parent / "logs").resolve())
    ]
    if vcf_stats_report:
        inner += ["--vcf_stats_report=true"]
    if par_bed:
        inner += [f"--par_regions_bed={par_bed}"]
    # Allow both a BED and a one-off region string (DeepVariant supports repeated --regions)
    if regions_bed:
        inner += [f"--regions={regions_bed}"]
    if region_str:
        inner += [f"--regions={region_str}"]
    if intermediate_dir:
        inner += ["--intermediate_results_dir", str(intermediate_dir)]
    if extra:
        inner += extra
    return inner

def build_glnexus_inner(
    gvcfs: List[Path],
    out_bcf: Path,
    config: Optional[str],
    config_yaml: Optional[Path],
    extra: List[str],
    threads: int
) -> List[str]:
    inner = ["glnexus_cli"]
    if config_yaml:
        inner += ["--config", str(config_yaml)]
    elif config:
        inner += ["--config", config]
    if threads > 0 and not any(x in ("--threads", "-t") for x in extra):
        inner += ["--threads", str(threads)]
    inner += extra
    inner += [str(p) for p in gvcfs]
    inner += ["-O", str(out_bcf)]
    return inner

def detect_cpus(default: int=1) -> int:
    try:
        import multiprocessing as mp
        return max(default, mp.cpu_count())
    except Exception:
        return default

# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser(
        prog="deepvariant_glnexus_batch.py",
        description="Run DeepVariant per-BAM and GLnexus cohort merge (Singularity/Docker/host GLnexus)."
    )
    # Core I/O
    ap.add_argument("-s","--samples", required=True, help="TSV/CSV with sample_id and BAM path.")
    ap.add_argument("-r","--reference", required=True, help="FASTA reference (indexed).")
    ap.add_argument("-o","--outdir", required=True, help="Output directory.")
    ap.add_argument("-O","--output-prefix", default=None, help="Prefix for per-sample outputs AND cohort outputs. Defaults: per-sample '<sample>.deepvariant', cohort 'cohort'.")

    # DeepVariant engine & options
    ap.add_argument("-E","--engine", choices=["singularity","docker"], required=True, help="Container engine for DeepVariant.")
    ap.add_argument("-I","--dv-image", required=True, help="DeepVariant image: path to .sif (Singularity) or Docker image name.")
    ap.add_argument("-M","--model-type", default="WGS", choices=["WGS","WES","PACBIO","ONT"], help="DeepVariant --model_type (default: WGS).")
    ap.add_argument("-p","--par-bed", default=None, help="BED of PAR regions (optional; passed to --par_regions_bed).")
    ap.add_argument("-R","--regions-bed", default=None, help="BED to limit calling (optional; DeepVariant --regions).")
    ap.add_argument("-g","--region", default=None, help="Single region string for DeepVariant --regions (e.g., chr20:10,000,000-10,010,000).")
    ap.add_argument("-S","--vcf-stats-report", action="store_true", help="Enable DeepVariant --vcf_stats_report=true.")
    ap.add_argument("-i","--intermediate-dir", default=None, help="Root directory for DeepVariant intermediate outputs. Per-sample subdirs will be created as needed.")
    ap.add_argument("-K","--keep-intermediates", action="store_true", help="Do NOT delete per-sample DeepVariant intermediate directories.")
    ap.add_argument("-n","--shards", type=int, default=0, help="DeepVariant --num_shards (default: auto nproc).")
    ap.add_argument("-u","--gpu", action="store_true", help="Enable GPU: adds --nv (Singularity) or --gpus all (Docker).")
    ap.add_argument("-x","--dv-extra", nargs="*", default=[], help="Extra args to pass to run_deepvariant.")

    # GLnexus engine
    ap.add_argument("-H","--glx-engine", choices=["container","host"], default="container", help="Run GLnexus in container or on host (PATH).")
    ap.add_argument("-G","--glx-image", default=None, help="GLnexus image (required if --glx-engine=container).")
    ap.add_argument("-C","--glnexus-config", default="deepvariantWGS", help="GLnexus preset (e.g., deepvariantWGS, deepvariantWES). Ignored if --glnexus-yaml is set.")
    ap.add_argument("-Y","--glnexus-yaml", default=None, help="Custom GLnexus YAML config (overrides --glnexus-config).")
    ap.add_argument("-X","--glnexus-extra", nargs="*", default=[], help="Extra args for glnexus_cli (e.g., --mem-gbytes 16).")

    # Threads
    ap.add_argument("-t","--threads", type=int, default=0, help="Threads for GLnexus and bcftools/tabix (default: auto nproc).")

    # Post-merge conversion
    ap.add_argument("-V","--no-bcf2vcf", action="store_true", help="Skip bcftools conversion BCF->VCF.GZ and tabix indexing.")

    # General
    ap.add_argument("-f","--force", action="store_true", help="Re-run even if outputs exist.")
    ap.add_argument("-N","--dry-run", action="store_true", help="Print planned commands; do not execute.")
    ap.add_argument("-c","--cohort-name", default="cohort", help="Legacy: Prefix for cohort outputs if --output-prefix is not set.")
    args = ap.parse_args()

    samples = read_samples(Path(args.samples))
    if not samples:
        raise SystemExit("No samples parsed from --samples.")

    reference = Path(args.reference).resolve()
    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)
    dv_logs = outdir / "logs" / "deepvariant"
    glx_logs = outdir / "logs" / "glnexus"
    bcftools_logs = outdir / "logs" / "bcftools"
    ensure_dir(dv_logs); ensure_dir(glx_logs); ensure_dir(bcftools_logs)

    # compute resources
    shards = (args.shards if args.shards > 0 else detect_cpus())
    threads = (args.threads if args.threads > 0 else detect_cpus())

    # option paths
    par_bed = Path(args.par_bed).resolve() if args.par_bed else None
    regions_bed = Path(args.regions_bed).resolve() if args.regions_bed else None
    region_str = args.region
    glx_yaml = Path(args.glnexus_yaml).resolve() if args.glnexus_yaml else None

    # Run DeepVariant per sample
    gvcfs = []
    plan = []
    per_sample_outputs = []

    for sid, bam in samples:
        bam = bam.resolve()
        s_out = outdir / sid
        ensure_dir(s_out)

        # Per-sample prefix
        sprefix = (f\"{args.output_prefix}.{sid}.dv\" if args.output_prefix else f\"{sid}.deepvariant\")
        vcf_out = s_out / f"{sprefix}.vcf.gz"
        gvcf_out = s_out / f"{sprefix}.g.vcf.gz"
        visual_out = s_out / f"{sprefix}.visual_report.html"

        # Intermediate dir assignment
        if args.intermediate_dir:
            inter_root = Path(args.intermediate_dir).resolve()
            inter_dir = inter_root / sid
        else:
            inter_dir = s_out / "intermediate_results"

        per_sample_outputs.append({
            "sample_id": sid,
            "prefix": sprefix,
            "vcf": str(vcf_out),
            "gvcf": str(gvcf_out),
            "visual_report": str(visual_out),
            "intermediate_dir": str(inter_dir),
        })

        gvcfs.append(gvcf_out)
        need_run = args.force or (not gvcf_out.exists())

        inner = build_deepvariant_inner(
            model_type=args.model_type,
            reference=reference,
            bam=bam,
            vcf_out=vcf_out,
            gvcf_out=gvcf_out,
            shards=shards,
            par_bed=par_bed,
            regions_bed=regions_bed,
            region_str=region_str,
            vcf_stats_report=args.vcf_stats_report,
            intermediate_dir=inter_dir,
            extra=args.dv_extra,
        )

        binds = path_binds_for([reference, bam, vcf_out.parent, gvcf_out.parent, inter_dir, par_bed, regions_bed])
        if args.engine == "singularity":
            cmd = build_singularity_cmd(args.dv_image, binds, inner, gpu=args.gpu)
        else:
            cmd = build_docker_cmd(args.dv_image, binds, inner, gpu=args.gpu, workdir=None)

        plan.append(("deepvariant", sid, cmd, dv_logs / f"{sid}.deepvariant.log", need_run))

        # Index per-sample outputs and preserve visual report
        idx_log = dv_logs / f"{sid}.tabix_index.log"
        idx_cmd = ["bash", "-lc", f"tabix -p vcf {shlex.quote(str(vcf_out))} && tabix -p vcf {shlex.quote(str(gvcf_out))}"]
        plan.append(("tabix", sid, idx_cmd, idx_log, True))

        # Try to copy any *visual_report.html* to the desired filename
        preserve_cmd = [
            "bash","-lc",
            "set -euo pipefail; "
            f"rep='' ; "
            f"for d in {shlex.quote(str(vcf_out.parent))} {shlex.quote(str(inter_dir))} {shlex.quote(str(vcf_out.parent / 'logs'))}; do "
            "  if [ -d \"$d\" ]; then "
            "    cand=$(ls -1 \"$d\"/*visual_report.html 2>/dev/null | head -n1 || true); "
            "    if [ -n \"$cand\" ]; then rep=\"$cand\"; break; fi; "
            "  fi; "
            "done; "
            f"if [ -n \"$rep\" ]; then cp -f \"$rep\" {shlex.quote(str(visual_out))}; fi"
        ]
        plan.append(("preserve", sid, preserve_cmd, dv_logs / f"{sid}.preserve_report.log", True))

        # Cleanup intermediates unless kept
        if not args.keep_intermediates:
            rm_cmd = ["bash", "-lc", f"rm -rf {shlex.quote(str(inter_dir))}"]
            plan.append(("cleanup", sid, rm_cmd, dv_logs / f"{sid}.cleanup.log", True))

    # GLnexus outputs (apply output-prefix if provided; else use cohort-name)
    cohort_prefix = (f\"{args.output_prefix}.cohort.dv.glnexus\" if args.output_prefix else args.cohort_name)
    out_bcf = outdir / f"{cohort_prefix}.glnexus.bcf"
    out_vcf_gz = outdir / f"{cohort_prefix}.glnexus.vcf.gz"

    # GLnexus merge
    glx_inner = build_glnexus_inner(
        gvcfs=gvcfs,
        out_bcf=out_bcf,
        config=None if glx_yaml else args.glnexus_config,
        config_yaml=glx_yaml,
        extra=args.glnexus_extra,
        threads=threads,
    )

    glx_need_run = args.force or (not out_bcf.exists())

    if args.glx_engine == "container":
        if not args.glx_image:
            raise SystemExit("--glx-image is required when --glx-engine=container")
        binds_glx = path_binds_for(gvcfs + [outdir] + ([glx_yaml] if glx_yaml else []))
        if args.engine == "singularity":
            glx_cmd = build_singularity_cmd(args.glx_image, binds_glx, glx_inner, gpu=False)
        else:
            glx_cmd = build_docker_cmd(args.glx_image, binds_glx, glx_inner, gpu=False, workdir=None)
    else:
        glx_cmd = glx_inner

    plan.append(("glnexus", cohort_prefix, glx_cmd, glx_logs / f"{cohort_prefix}.glnexus.log", glx_need_run))

    # bcftools conversion + tabix
    if not args.no_bcf2vcf:
        bcftools_need = args.force or (not out_vcf_gz.exists())
        thr = f"-@ {threads}" if threads > 0 else ""
        bcf2vcf_cmd = ["bash", "-lc", f"bcftools view {thr} -O z -o {shlex.quote(str(out_vcf_gz))} {shlex.quote(str(out_bcf))} && tabix {thr} -p vcf {shlexquote := shlex.quote(str(out_vcf_gz))} && echo {shlexquote} >/dev/null"]
        plan.append(("bcftools", "bcf2vcf", bcf2vcf_cmd, bcftools_logs / f"{cohort_prefix}.bcf2vcf.log", bcftools_need))

    # Save manifest
    manifest = {
        "engine": args.engine,
        "dv_image": args.dv_image,
        "glx_engine": args.glx_engine,
        "glx_image": args.glx_image,
        "reference": str(reference),
        "model_type": args.model_type,
        "par_bed": str(par_bed) if par_bed else None,
        "regions_bed": str(regions_bed) if regions_bed else None,
        "region": region_str,
        "shards": shards,
        "threads": threads,
        "keep_intermediates": bool(args.keep_intermediates),
        "output_prefix": args.output_prefix,
        "cohort_prefix": cohort_prefix,
        "per_sample_outputs": per_sample_outputs,
        "gvcfs": [str(p) for p in gvcfs],
        "glnexus": {
            "config": None if glx_yaml else args.glnexus_config,
            "yaml": str(glx_yaml) if glx_yaml else None,
            "out_bcf": str(out_bcf),
            "out_vcf_gz": str(out_vcf_gz) if not args.no_bcf2vcf else None,
        }
    }
    with (outdir / "manifest.deepvariant_glnexus.json").open("w") as mf:
        json.dump(manifest, mf, indent=2)

    # Execute or print
    for step, ident, cmd, log_file, need in plan:
        if not need:
            log(f"[skip] {step} {ident} already complete.")
            continue
        log(f"[run ] {step} {ident}")
        log(f"       log -> {log_file}")
        if args.dry_run:
            print("$ " + " ".join(shlex.quote(c) for c in cmd))
        else:
            run_checked(cmd, log_file)

    log("Done. Outputs:")
    for rec in per_sample_outputs:
        print("  gVCF:", rec["gvcf"])
        print("  VCF:", rec["vcf"])
        if os.path.exists(rec["visual_report"]):
            print("  Visual report:", rec["visual_report"])
    print("  Cohort BCF:", out_bcf)
    if not args.no_bcf2vcf:
        print("  Cohort VCF.GZ:", out_vcf_gz)

if __name__ == "__main__":
    main()
