#!/usr/bin/env python3
"""
ANGSD runner — fast per-chromosome parallelization with asyncio.

Per-chromosome mode (with --chr-list / -c):
  - Outputs in: <out>.perchr/
  - Per-chr prefix: <out>.perchr/<out>.<CHR>.*
  - Per-chr logs:   <out>.perchr/logs/<out>.<CHR>.log  (ANGSD stdout+stderr, streamed)
  - Master log:     <out>.perchr/<out>.progress.log    (also mirrors to stdout)
  - Concurrency:    --procs simultaneous ANGSD jobs, each with -P threads
  - OMP settings:   OMP_NUM_THREADS=-P for each child process

Whole-genome mode (no --chr-list):
  - Single run with prefix <out>

Major/minor handling
--------------------
- Default:                 -doMajorMinor 1
- With --polarize / -z:    -doMajorMinor 5   (requires --anc / -A)
- With --fixed-sites / -S: -doMajorMinor 3 -sites FILE  (mutually exclusive with --polarize)

Missingness & Depth (no silent defaults)
----------------------------------------
- --missingness / -m:
    * 0..1  → fraction missing allowed; minInd = floor((1 - frac) * N)
    * >=1   → explicit minInd integer
  If omitted → do not pass -minInd.

- --mean-cov / -M:
    * setMinDepthInd = floor(mean/3), min 1
    * setMaxDepthInd = ceil(2*mean),  min 1
  If omitted → do not pass -setMinDepthInd/-setMaxDepthInd.
"""

from __future__ import annotations

import argparse
import asyncio
import gzip
import math
import os
import sys
from pathlib import Path
from typing import List, Iterable, Tuple, Optional


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Execution / resources
    p.add_argument("-b", "--bamlist", required=True, help="BAM list file")
    p.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    p.add_argument("-A", "--anc", default=None, help="Ancestral FASTA (required if --polarize)")
    p.add_argument("-o", "--out", required=True,
                   help="Output prefix (global; per-chr adds .CHR and is placed under <out>.perchr/)")
    p.add_argument("-P", "--threads-per-job", type=int, default=10,
                   help="Threads for each ANGSD process (-P)")
    p.add_argument("-p", "--procs", type=int, default=2,
                   help="How many chromosomes to run concurrently (only when --chr-list is used)")

    # Sharding
    p.add_argument("-c", "--chr-list", default=None,
                   help="File with one chromosome name per line; if omitted, run a single whole-genome job")

    # Major/minor logic (mutually exclusive)
    g = p.add_mutually_exclusive_group()
    g.add_argument("-z", "--polarize", action="store_true",
                   help="Use ancestral polarization (-doMajorMinor 5); requires --anc")
    # --sites and --fixed are separate (can be used together)
    p.add_argument("-s", "--sites", default=None,
                   help="Sites file to restrict analysis (-sites FILE). Used with --fixed to set -doMajorMinor 3.")
    p.add_argument("-S", dest="sites", help=argparse.SUPPRESS)
    p.add_argument("-F", "--fixed", action="store_true", default=False,
               help="Toggle use of fixed major/minor info. When used together with --sites, sets -doMajorMinor 3.")

    # Missingness & depth (no defaults if not provided)
    p.add_argument("-m", "--missingness", default=None,
                   help="If 0..1 → fraction missing allowed (compute minInd). If >=1 → explicit minInd integer.")
    p.add_argument("-M", "--mean-cov", type=float, default=None,
                   help="Mean per-individual coverage -> setMinDepthInd=floor(mean/3), setMaxDepthInd=ceil(2*mean).")

    # Behavior
    p.add_argument("-f", "--overwrite", action="store_true",
                   help="Rerun a chromosome even if outputs exist")
    p.add_argument("-x", "--concat", action="store_true",
                   help="Concatenate per-chrom .mafs.gz and .beagle.gz at end (only with --chr-list)")

    p.add_argument("--genoDepth", type=int, default=5,
                   help="Minimum depth per site for genotypes (maps to --geno_minDepth). Default: 5")
    p.add_argument("--plink", action="store_true",
                   help="If set, include -doPlink 2 in ANGSD command to emit PLINK files")
    p.add_argument("--allow-improper", action="store_true",
                   help="If set, use -only_proper_pairs 0. Default behavior is 1 (only proper pairs).")
    p.add_argument("--postCutoff", type=float, default=0.95,
                   help="Posterior cutoff for genotype/posterior filtering (ANGSD -postCutoff). Default: 0.95")
    p.add_argument("-N", "--ngsrelate", action="store_true", default=False,
                   help="If set: prepare NGSrelate input by using -doGlf 3 and omitting -doGeno/geno_minDepth.")
    return p.parse_args()


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def log_progress(progress_log: Path, msg: str):
    line = msg if msg.endswith("\n") else msg + "\n"
    sys.stdout.write(line)
    sys.stdout.flush()
    with progress_log.open("a") as fh:
        fh.write(line)


def load_chromosomes(args) -> Optional[List[str]]:
    if args.chr_list:
        with open(args.chr_list) as fh:
            chrs = [line.strip() for line in fh if line.strip()]
        return chrs
    return None


def count_bams(bamlist_path: str) -> int:
    n = 0
    with open(bamlist_path) as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n

def compute_minInd(missingness_arg: Optional[str], n_bams: int) -> Optional[int]:
    if missingness_arg is None:
        return None
    try:
        val_int = int(missingness_arg)
        if val_int >= 1:
            return min(val_int, n_bams)
    except ValueError:
        pass
    try:
        frac = float(missingness_arg)
        if not (0.0 <= frac <= 1.0):
            raise ValueError
        present = math.floor((1.0 - frac) * n_bams)
        present = max(1, min(present, n_bams))
        return present
    except ValueError:
        raise SystemExit(f"--missingness must be an integer >=1 or a float in [0,1]. Got: {missingness_arg}")

def per_chr_outputs_exist(outprefix_chr: Path) -> bool:
    mafs = Path(str(outprefix_chr) + ".mafs.gz")
    return mafs.exists() and mafs.stat().st_size > 0


def depth_from_mean(mean_cov: Optional[float]) -> Optional[Tuple[int, int]]:
    if mean_cov is None:
        return None
    if mean_cov <= 0:
        raise SystemExit("--mean-cov must be > 0")
    dmin = max(1, math.floor(mean_cov / 3.0))
    dmax = max(1, math.ceil(2.0 * mean_cov))
    return dmin, dmax


def base_angsd_args(args, minInd_override: Optional[int], depth_pair: Optional[Tuple[int, int]]) -> List[str]:
    cmd = [
        "angsd",
        "-P", str(args.threads_per_job),
        "-b", args.bamlist,
        "-GL", "2",
        "-doCounts", "1",
        "-minMapQ", "30",
        "-minQ", "30",
        "-remove_bads", "1",
        "-uniqueOnly", "1",
        "-C", "50",
        "-baq", "1",
        "-minMaf", "0.05",
        "-doMaf", "2",
        "-SNP_pval", "1e-6",
        "-doglf", "2",
        "-doPost", "1",
        "-postCutoff", str(args.postCutoff),
        "-doGeno", "2",
        "-geno_minDepth", str(args.genoDepth),
        "-ref", args.ref,
        "-skipTriallelic", "1",
    ]

    # Adjust for NGSrelate mode
    if getattr(args, "ngsrelate", False):
        # Switch GLF to 3
        if "-doglf" in cmd:
            i = cmd.index("-doglf")
            if i+1 < len(cmd):
                cmd[i+1] = "3"
        else:
            cmd += ["-doglf", "3"]
        # Remove -doGeno and any geno depth flags
        while "-doGeno" in cmd:
            j = cmd.index("-doGeno")
            try:
                del cmd[j:j+2]
            except Exception:
                del cmd[j]
        for key in ("--geno_minDepth", "--genoDepth", "-geno_minDepth", "--geno-minDepth"):
            while key in cmd:
                j = cmd.index(key)
                if j+1 < len(cmd) and not cmd[j+1].startswith("-"):
                    del cmd[j:j+2]
                else:
                    del cmd[j:j+1]
    ##NGSRELATE_BLOCK##

    # Proper pairs policy (default 1; override to 0 with --allow-improper)
    if getattr(args, "allow_improper", False):
        cmd += ["-only_proper_pairs", "0"]
    else:
        cmd += ["-only_proper_pairs", "1"]


    # Major/minor selection
    if args.polarize:
        if not args.anc:
            raise SystemExit("--polarize requires --anc.")
        cmd += ["-doMajorMinor", "5", "-anc", args.anc]
    else:
        # If both --sites and --fixed are provided, use -doMajorMinor 3; else 1.
        if args.sites and args.fixed:
            cmd += ["-doMajorMinor", "3", "-sites", args.sites]
        else:
            cmd += ["-doMajorMinor", "1"]
        if args.anc:
            cmd += ["-anc", args.anc]


    # Only include minInd if requested
    if minInd_override is not None:
        cmd += ["-minInd", str(minInd_override)]

    # Only include per-individual depth thresholds if requested
    if depth_pair is not None:
        dmin, dmax = depth_pair
        cmd += ["-setMinDepthInd", str(dmin), "-setMaxDepthInd", str(dmax)]

    if getattr(args, "plink", False):
        cmd += ["-doPlink", "2"]
    return cmd
def build_angsd_cmd_perchr(args, chr_name: str, outprefix_chr: Path,
                           minInd_override: Optional[int], depth_pair: Optional[Tuple[int, int]]) -> List[str]:
    return base_angsd_args(args, minInd_override, depth_pair) + ["-out", str(outprefix_chr), "-r", chr_name]


def build_angsd_cmd_wholegenome(args, minInd_override: Optional[int], depth_pair: Optional[Tuple[int, int]]) -> List[str]:
    return base_angsd_args(args, minInd_override, depth_pair) + ["-out", args.out]


async def _stream_subprocess(cmd: List[str], log_file: Path, env: dict) -> int:
    """
    Run a subprocess and stream both stdout and stderr to `log_file` as they arrive.
    Returns the process return code.
    """
    with log_file.open("w") as lf:
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env=env,
        )

        async def _pump(stream, prefix=""):
            while True:
                line = await stream.readline()
                if not line:
                    break
                try:
                    s = line.decode("utf-8", "replace")
                except Exception:
                    s = str(line)
                if prefix:
                    lf.write(prefix + s)
                else:
                    lf.write(s)
                lf.flush()

        # FIX: use gather (not wait with raw coroutines)
        await asyncio.gather(
            _pump(proc.stdout),
            _pump(proc.stderr, prefix=""),
        )
        rc = await proc.wait()
        return rc


async def run_one_chr_async(cmd: List[str], chr_name: str, perchr_logfile: Path,
                            progress_log: Path, threads_per_job: int) -> Tuple[str, int]:
    log_progress(progress_log, f"[run ] {chr_name}: {' '.join(cmd)}")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads_per_job)
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")

    try:
        rc = await _stream_subprocess(cmd, perchr_logfile, env=env)
        if rc != 0:
            log_progress(progress_log, f"[error] {chr_name} failed (code {rc}); see {perchr_logfile}")
        else:
            log_progress(progress_log, f"[ok  ] {chr_name} completed; log: {perchr_logfile}")
        return (chr_name, rc)
    except Exception as e:
        with perchr_logfile.open("a") as lf:
            lf.write(f"\n[exception] {e}\n")
        log_progress(progress_log, f"[error] {chr_name} exception: {e}; see {perchr_logfile}")
        return (chr_name, 1)


async def run_whole_genome_async(cmd: List[str], threads_per_job: int) -> int:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads_per_job)
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        env=env,
    )

    async def _pipe(src, dst):
        while True:
            line = await src.readline()
            if not line:
                break
            dst.buffer.write(line)
            dst.flush()

    # FIX: use gather (not wait with raw coroutines)
    await asyncio.gather(
        _pipe(proc.stdout, sys.stdout),
        _pipe(proc.stderr, sys.stderr),
    )
    return await proc.wait()


def concat_gz_text(inputs: Iterable[Path], output: Path):
    inputs = [p for p in inputs if p.exists() and p.stat().st_size > 0]
    if not inputs:
        return
    with gzip.open(output, "wt") as outfh:
        wrote_header = False
        for path in inputs:
            with gzip.open(path, "rt") as infh:
                for j, line in enumerate(infh):
                    if j == 0:
                        if not wrote_header:
                            outfh.write(line)
                            wrote_header = True
                    else:
                        outfh.write(line)


def maybe_concat(outprefix_dir: Path, outprefix: str, chrs: List[str], progress_log: Path):
    mafs_inputs = [outprefix_dir / f"{outprefix}.{c}.mafs.gz" for c in chrs]
    beagle_inputs = [outprefix_dir / f"{outprefix}.{c}.beagle.gz" for c in chrs]
    mafs_out = outprefix_dir / f"{outprefix}.autosomes.mafs.gz"
    beagle_out = outprefix_dir / f"{outprefix}.autosomes.beagle.gz"

    log_progress(progress_log, "[post] Concatenating .mafs.gz ...")
    concat_gz_text(mafs_inputs, mafs_out)
    log_progress(progress_log, f"[post] Wrote {mafs_out}")

    log_progress(progress_log, "[post] Concatenating .beagle.gz ...")
    concat_gz_text(beagle_inputs, beagle_out)
    log_progress(progress_log, f"[post] Wrote {beagle_out}")


async def main_async():
    args = parse_args()

    if args.polarize and args.fixed_sites:
        raise SystemExit("--polarize and --fixed-sites are mutually exclusive.")

    n_bams = count_bams(args.bamlist)
    minInd_override = compute_minInd(args.missingness, n_bams) if args.missingness is not None else None
    depth_pair = depth_from_mean(args.mean_cov)

    chrs = load_chromosomes(args)
    if chrs:
        perchr_dir = Path(f"{args.out}.perchr")
        logs_dir = perchr_dir / "logs"
        ensure_dir(perchr_dir)
        ensure_dir(logs_dir)
        progress_log = perchr_dir / f"{args.out}.progress.log"

        log_progress(progress_log, f"[info] Mode: per-chromosome from --chr-list ({len(chrs)} entries)")
        log_progress(progress_log, f"[info] Concurrency: {args.procs} ; threads per job (-P): {args.threads_per_job}")

        jobs: List[Tuple[str, List[str], Path]] = []
        for c in chrs:
            outprefix_chr = perchr_dir / f"{args.out}.{c}"
            if (not args.overwrite) and per_chr_outputs_exist(outprefix_chr):
                log_progress(progress_log, f"[skip] {c}: outputs already present.")
                continue
            cmd = build_angsd_cmd_perchr(args, c, outprefix_chr, minInd_override, depth_pair)
            perchr_logfile = logs_dir / f"{args.out}.{c}.log"
            jobs.append((c, cmd, perchr_logfile))

        sem = asyncio.Semaphore(args.procs)
        failures: List[str] = []
        results: List[Tuple[str, int]] = []

        async def worker(c: str, cmd: List[str], logfile: Path):
            async with sem:
                rc = await run_one_chr_async(cmd, c, logfile, progress_log, args.threads_per_job)
                results.append(rc)

        tasks = [asyncio.create_task(worker(c, cmd, lf)) for (c, cmd, lf) in jobs]
        if tasks:
            await asyncio.gather(*tasks)

        for c, rc in results:
            if rc != 0:
                failures.append(c)

        succeeded = [c for (c, _, _) in jobs if c not in failures]

        if args.concat and succeeded:
            maybe_concat(perchr_dir, args.out, succeeded, progress_log)

        if failures:
            log_progress(progress_log, f"[warn] Failed chromosomes: {', '.join(sorted(failures))}")
            log_progress(progress_log, "[done] Completed with failures.")
            raise SystemExit(1)

        log_progress(progress_log, "[done] All requested chromosomes processed.")
    else:
        print(f"[info] Mode: whole-genome (no --chr-list provided)")
        cmd = build_angsd_cmd_wholegenome(args, minInd_override, depth_pair)
        rc = await run_whole_genome_async(cmd, args.threads_per_job)
        if rc != 0:
            raise SystemExit(rc)
        print("[done] Whole-genome job completed.")


def main():
    try:
        asyncio.run(main_async())
    except KeyboardInterrupt:
        print("\n[abort] Interrupted by user.", file=sys.stderr)
        sys.exit(130)


if __name__ == "__main__":
    main()