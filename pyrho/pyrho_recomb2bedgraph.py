#!/usr/bin/env python3
import sys
import argparse
import gzip
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Convert recombination rate to BEDGRAPH format for plotting.")
    parser.add_argument("recomb_file", help="Input recombination file (TSV: chrom start end rate)")
    parser.add_argument("--window", type=int, required=True, help="Window size in bp")
    parser.add_argument("--out", help="Output BEDGRAPH file (TSV format)")
    parser.add_argument("--agp", help="Optional AGP file for liftover (not implemented)")
    parser.add_argument("--chain", help="Optional CHAIN or CHAIN.GZ file for liftover (not implemented)")
    return parser.parse_args()

def open_chain_file(chain_path):
    if chain_path.endswith(".gz"):
        return gzip.open(chain_path, "rt")
    return open(chain_path, "r")

def emit_window(out_stream, chrom, start, end, window_data):
    if not window_data:
        return
    total_len = sum(e - b for b, e, r in window_data)
    if total_len == 0:
        return
    weighted_rate = sum(r * (e - b) for b, e, r in window_data) / total_len
    print(f"{chrom}\t{start}\t{end}\t{weighted_rate:.6f}", file=out_stream)

def main():
    args = parse_args()
    window_size = args.window
    out_stream = open(args.out, "w") if args.out else sys.stdout

    if args.agp:
        print(f"[INFO] AGP file provided: {args.agp}", file=sys.stderr)
        # TODO: Parse AGP for coordinate transformation

    if args.chain:
        print(f"[INFO] CHAIN file provided: {args.chain}", file=sys.stderr)
        try:
            with open_chain_file(args.chain) as chain_file:
                # Placeholder: read lines for later transformation
                chain_lines = chain_file.readlines()
                # TODO: implement coordinate liftover
        except Exception as e:
            print(f"[ERROR] Could not open chain file: {e}", file=sys.stderr)
            sys.exit(1)

    curr_chr = ''
    curr_window = []
    window_start = 0

    with open(args.recomb_file, 'r') as in_file:
        for line in in_file:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            try:
                chrom, begin, end, rate = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
            except ValueError:
                continue  # header or malformed line

            # Placeholder: if liftover active, apply remapping to chrom, begin, end

            if chrom != curr_chr:
                if curr_chr and curr_window:
                    emit_window(out_stream, curr_chr, window_start, curr_window[-1][1], curr_window)
                curr_chr = chrom
                curr_window = []
                window_start = 0

            while end > window_start + window_size:
                curr_window.append((begin, window_start + window_size, rate))
                emit_window(out_stream, curr_chr, window_start, window_start + window_size, curr_window)
                begin = window_start + window_size
                window_start += window_size
                curr_window = []

            curr_window.append((begin, end, rate))

    if curr_window:
        emit_window(out_stream, curr_chr, window_start, curr_window[-1][1], curr_window)

    if args.out:
        out_stream.close()

if __name__ == "__main__":
    main()
