#!/usr/bin/env python3

import argparse
import os
import subprocess
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess Whole Genome Data: Adapter Trimming (fastp) and Decontamination (Kraken2)"
    )
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use")
    parser.add_argument("-s", required=True, help="Sample sheet in TSV format")
    parser.add_argument("-f", help="Path to fastp")
    parser.add_argument("-k", help="Path to Kraken2")
    parser.add_argument("-d", help="Path to Kraken2 database")
    parser.add_argument("-o", required=True, help="Output directory")
    parser.add_argument("--summary-only", action="store_true", help="Only summarize existing Kraken2 and fastp reports")
    parser.add_argument("--summary-output", type=str, help="File name to save contamination summary table")
    return parser.parse_args()

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_fastp(fastp, threads, rawread1, rawread2, sampleID, outdir):
    trimmed_1 = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_1.trimmed.fastq.gz")
    trimmed_2 = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_2.trimmed.fastq.gz")
    json_out = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_fastp.json")
    html_out = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_fastp.html")

    subprocess.run([
        fastp, "-w", str(threads), "-i", rawread1, "-I", rawread2,
        "-o", trimmed_1, "-O", trimmed_2, "-j", json_out, "-h", html_out
    ], check=True)
    return json_out

def extract_fastp_statistics(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            total_reads = data['summary']['before_filtering']['total_reads']
            filtered_reads = data['summary']['after_filtering']['total_reads']
            percent_retained = (filtered_reads / total_reads) * 100 if total_reads > 0 else 0
            return total_reads, filtered_reads, percent_retained
    except Exception as e:
        print(f"[Warning] Could not read fastp JSON: {json_file}. Error: {e}")
        return None, None, None

def run_kraken2(kraken2, kraken2_db, threads, sampleID, outdir):
    trimmed_1 = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_1.trimmed.fastq.gz")
    trimmed_2 = os.path.join(outdir, "trimmed_fastq", f"{sampleID}_2.trimmed.fastq.gz")
    report = os.path.join(outdir, "decon_fastq", f"{sampleID}_kraken2_report.txt")
    unclassified = os.path.join(outdir, "decon_fastq", f"{sampleID}#.decon.fastq")

    subprocess.run([
        kraken2, "--threads", str(threads), "--db", kraken2_db,
        "--report", report, "--use-names", "--unclassified-out", unclassified,
        "--paired", trimmed_1, trimmed_2
    ], check=True)

    for file in os.listdir(os.path.join(outdir, "decon_fastq")):
        if file.startswith(sampleID) and file.endswith(".fastq"):
            subprocess.run(["pigz", "-p", str(threads), os.path.join(outdir, "decon_fastq", file)], check=True)
    return report

def extract_classification_percentages(report_file, filtered_reads):
    try:
        with open(report_file, 'r') as f:
            for line in f:
                if "unclassified" in line:
                    unclassified_pct = float(line.strip().split()[0])
                    contaminated_pct = 100.0 - unclassified_pct
                    uncontaminated_reads = int((unclassified_pct / 100.0) * filtered_reads)
                    contaminated_reads = filtered_reads - uncontaminated_reads
                    return unclassified_pct, contaminated_pct, uncontaminated_reads, contaminated_reads
    except Exception as e:
        print(f"[Warning] Could not read Kraken2 report: {report_file}. Error: {e}")
    return None, None, None, None

def print_summary_table(summary_data, out_file=None):
    header = "{:<15} {:>10} {:>10} {:>10} {:>12} {:>18} {:>12} {:>18}".format(
        "SampleID", "Raw", "Trimmed", "%Trimmed",
        "Uncontam.", "Uncontam. (%)", "Contam.", "Contam. (%)"
    )
    separator = "-" * 120

    output_lines = [header, separator]
    for row in summary_data:
        (sampleID, total, trimmed, retained, good, good_pct, bad, bad_pct) = row
        line = "{:<15} {:>10} {:>10} {:>10.2f} {:>12} {:>18.2f} {:>12} {:>18.2f}".format(
            sampleID,
            total if total is not None else "N/A",
            trimmed if trimmed is not None else "N/A",
            retained if retained is not None else 0,
            good if good is not None else "N/A",
            good_pct if good_pct is not None else 0,
            bad if bad is not None else "N/A",
            bad_pct if bad_pct is not None else 0
        )
        output_lines.append(line)

    # Print to stdout
    print("\nContamination Summary:")
    for line in output_lines:
        print(line)

    # Optionally write to file
    if out_file:
        with open(out_file, "w") as f:
            for line in output_lines:
                f.write(line + "\n")
        print(f"\nSummary saved to: {out_file}")

def main():
    args = parse_args()
    summary_data = []

    if not args.summary_only:
        ensure_directory(args.o)
        ensure_directory(os.path.join(args.o, "trimmed_fastq"))
        ensure_directory(os.path.join(args.o, "decon_fastq"))

    with open(args.s, 'r') as f:
        for line in f:
            sampleID, rawread1, rawread2 = line.strip().split('\t')
            json_file = os.path.join(args.o, "trimmed_fastq", f"{sampleID}_fastp.json")
            report_file = os.path.join(args.o, "decon_fastq", f"{sampleID}_kraken2_report.txt")

            if not args.summary_only:
                json_file = run_fastp(args.f, args.threads, rawread1, rawread2, sampleID, args.o)
                report_file = run_kraken2(args.k, args.d, args.threads, sampleID, args.o)

            total_reads, filtered_reads, retained_pct = extract_fastp_statistics(json_file)
            unclassified_pct, contaminated_pct, unclassified_reads, contaminated_reads = extract_classification_percentages(
                report_file, filtered_reads if filtered_reads else 0)

            summary_data.append((
                sampleID, total_reads, filtered_reads, retained_pct,
                unclassified_reads, unclassified_pct, contaminated_reads, contaminated_pct
            ))

    out_path = os.path.join(args.o, args.summary_output) if args.summary_output else None
    print_summary_table(summary_data, out_file=out_path)

if __name__ == "__main__":
    main()
