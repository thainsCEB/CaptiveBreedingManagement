#!/usr/bin/env python3
"""
Read a two-column table (ACCESSION, SAMPLE_NAME), peek the first FASTQ header
for each accession with fastq-dump -X 1 -Z, and write a three-column TSV:

  accession  sampleName  headerType

where headerType ∈ {instrument, numeric, unknown}.

- instrument : colon-delimited instrument-style tokens (e.g., HWI:...:... or known prefixes)
- numeric    : headers like @12345 or @12345/1
- unknown    : anything else or errors (network, permissions, etc.)

Usage:
  python classify_sra_headers_3col.py -t runs.tsv -o runs.3col.tsv \
      --threads 8 --sratoolkit /path/to/sratoolkit   # optional

Requires:
  - SRA Toolkit (fastq-dump) available on PATH or via --sratoolkit
"""
import argparse, csv, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed

INSTRUMENT_PREFIXES = (
    "HWI","HWUSI","MISEQ","HISEQ","NEXTSEQ","NOVASEQ","NS","NB",
    "DNBSEQ","BGISEQ","ILLUMINA","A00","MN","MINION","GRIDION",
    "PROMETHION","PROM"
)
SRA_RUN_RE = re.compile(r'^[SED]RR[0-9]+$', re.IGNORECASE)

def parse_args():
    ap = argparse.ArgumentParser(description="Classify SRA header style into 3-column TSV.")
    ap.add_argument("-t","--table", required=True, help="Two-column input: ACCESSION SAMPLE_NAME")
    ap.add_argument("-o","--out-tsv", required=True, help="Output 3-column TSV path")
    ap.add_argument("--sratoolkit", default="", help="Path to SRA toolkit root or fastq-dump")
    ap.add_argument("--threads", type=int, default=4, help="Parallel workers (default: 4)")
    return ap.parse_args()

def resolve_fastq_dump(hint: str) -> str:
    cand = (os.path.join(hint, "bin", "fastq-dump") if hint and os.path.isdir(hint) else hint) if hint else None
    exe = None
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        exe = cand
    else:
        exe = shutil.which("fastq-dump") if not cand else shutil.which(cand)
    if not exe:
        sys.exit("ERROR: fastq-dump not found (use --sratoolkit or add to PATH).")
    return exe

def read_pairs(path):
    with open(path, newline="") as fh:
        text = fh.read()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return []
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    rows = list(csv.reader(lines, delimiter=delim)) if delim else [re.split(r"\s+", ln) for ln in lines]
    pairs = []
    for i,row in enumerate(rows):
        if len(row) < 2: continue
        acc, sample = row[0].strip(), row[1].strip()
        if i == 0 and not SRA_RUN_RE.match(acc):  # skip header row
            continue
        pairs.append((acc, sample))
    return pairs

def peek_header(fastq_dump: str, acc: str):
    cmd = [fastq_dump, "-X","1","-Z","--skip-technical",
           "--defline-seq","@$sn[_$rn]/$ri","--defline-qual","+",
           acc]
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, text=True)
        if p.returncode != 0:
            return "", f"fastq-dump rc={p.returncode}: {p.stderr.strip()[:300]}"
        lines = p.stdout.splitlines()
        return (lines[0].strip() if lines else ""), ""
    except Exception as e:
        return "", f"exception: {e}"

def classify_header_line(h: str) -> str:
    if not h or not h.startswith("@"):
        return "unknown"
    # Instrument-style: allow leading letters OR digits/underscores before colon; require >=2 colon groups
    if re.search(r"^@[^ /:]*(?::[^ /:]*){2,}(?:/[12])?$", h):
        return "instrument"
    first = h[1:].split("/", 1)[0].upper()
    if any(first.startswith(pfx) for pfx in INSTRUMENT_PREFIXES):
        return "instrument"
    # Heuristic: if there are many colon-delimited fields before any whitespace, it's instrument-like
    head = h.split()[0]
    if head.count(":") >= 4:
        return "instrument"
    # Numeric-only
    if re.match(r"^@[0-9]+(?:[/_][12])?$", h):
        return "numeric"
    return "unknown"

def main():
    args = parse_args()
    fastq_dump = resolve_fastq_dump(args.sratoolkit)
    pairs = read_pairs(args.table)
    if not pairs:
        sys.exit("ERROR: No (accession, sampleName) pairs parsed.")

    uniq = sorted({a for a,_ in pairs})
    results = {}
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        fut2acc = {ex.submit(peek_header, fastq_dump, acc): acc for acc in uniq}
        for fut in as_completed(fut2acc):
            acc = fut2acc[fut]
            header, err = fut.result()
            typ = classify_header_line(header) if not err else "unknown"
            results[acc] = typ
            print(f"[classify] {acc}\t{typ}\t{header or err}", file=sys.stderr)

    with open(args.out_tsv, "w", newline="") as out:
        out.write("accession\tsampleName\theaderType\n")
        for acc, sample in pairs:
            out.write(f"{acc}\t{sample}\t{results.get(acc,'unknown')}\n")

if __name__ == "__main__":
    main()
