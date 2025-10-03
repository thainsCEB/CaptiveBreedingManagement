#!/usr/bin/env python3

import argparse
import tempfile
import subprocess
import os

def parse_roh_line(line):
    """Parse a bcftools roh line with RG prefix and return a BED entry."""
    fields = line.strip().split()
    if fields[0] != "RG" or len(fields) < 5:
        raise ValueError("Invalid ROH line format")
    sample = fields[1]
    chrom = fields[2]
    start = int(fields[3]) - 1  # Convert to 0-based BED
    end = int(fields[4])        # BED end is exclusive
    return (chrom, start, end, sample)

def write_bed(input_file, output_file, include_header, min_length):
    """Parse ROH file and write BED file."""
    with open(input_file) as infile, open(output_file, 'w') as out:
        if include_header:
            out.write("chrom\tstart\tend\tsample\n")
        for line in infile:
            if line.startswith("#") or not line.startswith("RG"):
                continue
            try:
                chrom, start, end, sample = parse_roh_line(line)
                if (end - start) >= min_length:
                    out.write(f"{chrom}\t{start}\t{end}\t{sample}\n")
            except Exception as e:
                print(f"Skipping line: {line.strip()} ({e})")

def merge_bed(input_bed, output_bed):
    """Sort and merge a BED file using bedtools."""
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
        sorted_file = tmp.name
    try:
        # Remove header if present
        with open(input_bed) as f_in, open(sorted_file, "w") as f_out:
            for line in f_in:
                if line.startswith("chrom"):  # skip header
                    continue
                f_out.write(line)
        # Sort
        sort_cmd = f"sort -k1,1 -k2,2n {sorted_file}"
        sorted_bed = sorted_file + ".sorted"
        with open(sorted_bed, "w") as f:
            subprocess.run(sort_cmd, shell=True, stdout=f, check=True)

        # Merge
        merge_cmd = f"bedtools merge -i {sorted_bed}"
        with open(output_bed, "w") as f:
            subprocess.run(merge_cmd, shell=True, stdout=f, check=True)
    finally:
        os.remove(sorted_file)
        if os.path.exists(sorted_file + ".sorted"):
            os.remove(sorted_file + ".sorted")

def main():
    parser = argparse.ArgumentParser(description="Convert bcftools ROH output (RG format) to BED format.")
    parser.add_argument("-i", "--input", required=True, help="Input ROH file from bcftools")
    parser.add_argument("-o", "--output", required=True, help="Output BED file for all samples")
    parser.add_argument("--include-header", action="store_true", help="Include BED header line")
    parser.add_argument("--merge", action="store_true", help="Also produce merged BED file across samples")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum ROH length to retain (default: 0)")
    args = parser.parse_args()

    # Write raw BED
    write_bed(args.input, args.output, args.include_header, args.min_length)

    # Optionally merge
    if args.merge:
        merged_output = os.path.splitext(args.output)[0] + ".merged.bed"
        print(f"Generating merged BED: {merged_output}")
        merge_bed(args.output, merged_output)

if __name__ == "__main__":
    main()
