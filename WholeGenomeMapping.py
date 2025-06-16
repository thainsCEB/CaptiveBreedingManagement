#!/usr/bin/env python3

import argparse
import subprocess
import csv
from pathlib import Path
import shutil

def run_command(cmd, log_file=None):
    print("Running command: {}".format(cmd))
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "w") as lf:
            process = subprocess.run(cmd, shell=True, stdout=lf, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True)

def bwa_index_exists(reference):
    """Check if BWA index files exist for a given reference"""
    index_extensions = [".bwt", ".pac", ".ann", ".amb", ".sa"]
    return all(Path(reference + ext).exists() for ext in index_extensions)

def bam_has_reads(bam_path):
    result = subprocess.run(
        ["samtools", "idxstats", str(bam_path)], capture_output=True, text=True, check=True
    )
    total_mapped = sum(int(line.split('\t')[2]) for line in result.stdout.strip().split('\n') if line)
    return total_mapped > 0

def main():
    parser = argparse.ArgumentParser(description="Whole Genome Mapping Pipeline")
    parser.add_argument("-t", "--threads", type=int, required=True, help="Number of threads to use")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format")
    parser.add_argument("-s", "--samplesheet", required=True, help="Sample Sheet (TSV: ID, read1, read2)")
    parser.add_argument("-o", "--output", required=True, help="Output directory path")
    args = parser.parse_args()

    ref = args.reference
    threads = args.threads
    samplesheet = args.samplesheet
    outdir = Path(args.output)
    tempdir = outdir / "tempbam"
    finalbams = outdir / "final_BAMs"
    stats_dir = outdir / "stats"
    metrics_dir = stats_dir / "metrics"
    logs_dir = outdir / "logs"

    # Create output directories
    for directory in [outdir, tempdir, finalbams, stats_dir, metrics_dir, logs_dir]:
        if not directory.exists():
            directory.mkdir(parents=True)

    # Setup paths for reference support files
    ref_path = Path(ref)
    ref_base = ref_path.stem.replace('.fa', '').replace('.fna', '').replace('.fasta', '')
    ref_dir = ref_path.parent
    dict_path = ref_dir / "{}.dict".format(ref_base)

    # Check for BWA index
    bwa_log = logs_dir / "bwa_index.log"
    if not bwa_index_exists(ref):
        print("BWA index files not found for reference. Indexing with BWA...")
        run_command(f"bwa index {ref}", log_file=bwa_log)
    else:
        print("BWA index files found. Proceeding...")

    # Check for sequence dictionary
    dict_log = logs_dir / "picard_create_dict.log"
    if not dict_path.exists():
        print("Sequence Dictionary not found. Generating with Picard...")
        run_command(f"picard CreateSequenceDictionary R={ref} O={dict_path}", log_file=dict_log)
    else:
        print("Sequence Dictionary found. Proceeding...")

    # Check for FASTA index (.fai)
    fai_log = logs_dir / "samtools_faidx.log"
    fai_path = Path("{}.fai".format(ref))
    if not fai_path.exists():
        print("FASTA index (.fai) not found. Generating with Samtools...")
        run_command(f"samtools faidx {ref}", log_file=fai_log)
    else:
        print("FASTA index (.fai) found. Proceeding...")

    # Process each sample
    with open(samplesheet) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                print("Skipping incomplete row: {}".format(row))
                continue
            prefix, read1, read2 = row[0], row[1], row[2]
            print(f"\nStarting the read alignment for {prefix}")

            bam = tempdir / f"{prefix}.bam"
            markdup_bam = tempdir / f"{prefix}.markdup.bam"
            realn_bam = tempdir / f"{prefix}.realn.bam"
            filt_bam = tempdir / f"{prefix}.filt.bam"
            fix_bam = finalbams / f"{prefix}.fix.bam"
            metrics_file = metrics_dir / f"{prefix}.metrics.txt"

            # Logs per step
            logs = {
                "bwa_mem": logs_dir / f"{prefix}_bwa_mem.log",
                "markdup": logs_dir / f"{prefix}_markdup.log",
                "leftalign": logs_dir / f"{prefix}_leftalign.log",
                "printreads": logs_dir / f"{prefix}_printreads.log",
                "fixmate": logs_dir / f"{prefix}_fixmate.log",
                "index": logs_dir / f"{prefix}_index.log",
                "quickcheck": logs_dir / f"{prefix}_quickcheck.log",
                "idxstats": logs_dir / f"{prefix}_idxstats.log",
            }

            # BWA alignment + sorting
            run_command(
                f"bwa mem -t {threads} -R '@RG\\tID:{prefix}\\tSM:{prefix}\\tLB:ILLUMINA\\tPL:ILLUMINA' {ref} {read1} {read2} "
                f"| samtools view -@ {threads} -bhS "
                f"| samtools sort -@ {threads} -o {bam} -",
                log_file=logs["bwa_mem"]
            )

            # Mark duplicates
            run_command(
                f"picard MarkDuplicates -Xmx20G I={bam} O={markdup_bam} METRICS_FILE={metrics_file} REMOVE_DUPLICATES=TRUE",
                log_file=logs["markdup"]
            )

            # GATK LeftAlignIndels
            run_command(
                f"gatk --java-options '-Xmx40G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' LeftAlignIndels "
                f"-R {ref} -I {markdup_bam} -O {realn_bam}",
                log_file=logs["leftalign"]
            )

            # GATK PrintReads
            run_command(
                f"gatk --java-options '-Xmx40G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' PrintReads "
                f"-R {ref} -I {realn_bam} -O {filt_bam}",
                log_file=logs["printreads"]
            )

            # Fixmate and index final BAM
            run_command(
                f"samtools sort -@ {threads} -n {filt_bam} "
                f"| samtools fixmate -@ {threads} - - "
                f"| samtools sort -@ {threads} -o {fix_bam} -",
                log_file=logs["fixmate"]
            )
            run_command(f"samtools index -@ {threads} {fix_bam}", log_file=logs["index"])

            # BAM integrity check
            bam_ok = False
            try:
                run_command(f"samtools quickcheck {fix_bam}", log_file=logs["quickcheck"])
                # For idxstats, capture output to check reads
                result = subprocess.run(
                    ["samtools", "idxstats", str(fix_bam)],
                    capture_output=True, text=True, check=True
                )
                with open(logs["idxstats"], "w") as lf:
                    lf.write(result.stdout)
                total_mapped = sum(int(line.split('\t')[2]) for line in result.stdout.strip().split('\n') if line)
                if total_mapped > 0:
                    bam_ok = True
                else:
                    print(f"Warning: {fix_bam} is indexed but contains no mapped reads!")
            except subprocess.CalledProcessError:
                print(f"Warning: {fix_bam} may be corrupted or incomplete.")

            if bam_ok:
                print(f"{fix_bam} is indexed and contains mapped reads.")
                # Delete all files in tempbam directory
                for temp_file in tempdir.iterdir():
                    try:
                        if temp_file.is_file() or temp_file.is_symlink():
                            temp_file.unlink()
                        elif temp_file.is_dir():
                            shutil.rmtree(temp_file)
                    except Exception as e:
                        print(f"Warning: Could not delete {temp_file}: {e}")

            print(f"Read Mapping to the Reference with BWA mem is completed for {prefix}")

if __name__ == "__main__":
    main()
