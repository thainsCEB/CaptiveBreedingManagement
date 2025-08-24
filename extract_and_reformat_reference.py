#!/usr/bin/env python3
import argparse, gzip, io, os, re, sys
from typing import Set, Optional

def smart_open_read(path: str):
    if path == "-" or path is None: return sys.stdin.buffer
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")

def smart_open_write(path: str):
    if path == "-" or path is None: return sys.stdout.buffer
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    return gzip.open(path, "wb") if path.endswith(".gz") else open(path, "wb")

def to_text(stream):
    return stream if isinstance(stream, io.TextIOBase) else io.TextIOWrapper(stream, encoding="utf-8", newline="")

# -------- Cleaning / parsing helpers --------
CHR_PATTS = [
    re.compile(r"\bchr[_\- ]?([0-9]+|[xy])\b", re.IGNORECASE),
    re.compile(r"\bchromosome[_\- ]*([0-9]+|[xy])\b", re.IGNORECASE),
    re.compile(r"\bchrom(?:osome)?([0-9]+|[xy])\b", re.IGNORECASE),  # handles "chrom16"
]
ACC_VER = re.compile(r"^([A-Za-z]{2}_[0-9]+)(?:\.[0-9]+)?$")  # e.g. NC_072414.2 -> NC_072414

def looks_like_mt(header: str) -> bool:
    h = header
    if re.search(r"\bmt\b", h, re.IGNORECASE): return True
    if re.search(r"mitochondr", h, re.IGNORECASE): return True
    if re.search(r"\bmitogenome\b", h, re.IGNORECASE): return True
    return False

def cleaned_name_from_header(header_no_gt: str) -> Optional[str]:
    if looks_like_mt(header_no_gt):
        return "MT"
    for rp in CHR_PATTS:
        m = rp.search(header_no_gt)
        if m:
            tok = m.group(1).upper()
            if tok in {"X","Y"} or tok.isdigit():
                return tok
    return None

def strip_version(token: str) -> str:
    m = ACC_VER.match(token)
    return m.group(1) if m else token

# -------- Targets normalization --------
def normalize_targets(list_path: str) -> Set[str]:
    """
    Produce a LOWERCASED set that includes user inputs plus sensible variants:
    - drop leading '>'
    - add/no 'chr' versions
    - strip accession version (NC_072414.2 -> NC_072414)
    - uppercase short names handled via lowercasing the set
    """
    raw: Set[str] = set()
    with to_text(smart_open_read(list_path)) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"): continue
            raw.add(s)

    norm: Set[str] = set()
    for t in raw:
        t0 = t.strip()
        if t0.startswith(">"): t0 = t0[1:].strip()

        # base
        norm.add(t0.lower())

        # strip version from accession-like token
        t0_nover = strip_version(t0)
        norm.add(t0_nover.lower())

        # add/remove chr prefix for simple chromosomes
        m = re.fullmatch(r"(?:chr)?([0-9]+|[xy]|mt)$", t0, flags=re.IGNORECASE)
        if m:
            short = m.group(1).lower()
            norm.add(short)              # e.g. "16"
            norm.add(("chr"+short).lower())  # "chr16"
    return norm

# -------- Candidate generation for each header --------
def record_candidates(header_no_gt: str) -> Set[str]:
    """
    Return a LOWERCASED set of match candidates derived from the header.
    """
    full = header_no_gt
    first = header_no_gt.split()[0] if header_no_gt else ""

    cands: Set[str] = set()

    # full header and first token (with & without version)
    cands.add(full.lower())
    if first:
        cands.add(first.lower())
        cands.add(strip_version(first).lower())

    # cleaned short name (16/X/Y/MT) if detectable
    cshort = cleaned_name_from_header(header_no_gt)
    if cshort:
        cands.add(cshort.lower())
        cands.add(("chr"+cshort).lower())

    # Also scrape any chr patterns appearing anywhere (adds robustness to odd formatting)
    for rp in CHR_PATTS:
        for m in rp.finditer(header_no_gt):
            tok = m.group(1).upper()
            cands.add(tok.lower())
            cands.add(("chr"+tok).lower())

    # If looks mito, ensure MT is included
    if looks_like_mt(header_no_gt):
        cands.add("mt")

    return cands

# -------- Extraction --------
def run_extract(in_fa: str, out_fa: str, list_path: str, clean: bool, debug: int):
    targets = normalize_targets(list_path)

    fin = smart_open_read(in_fa)
    fout = smart_open_write(out_fa)

    n_seen = 0
    n_written = 0
    debug_printed = 0

    try:
        with to_text(fin) as inp, fout:
            header = None
            seq_lines = []

            def flush():
                nonlocal n_written, debug_printed
                if header is None: return
                cands = record_candidates(header)
                hit = any((c in targets) for c in cands)

                if debug_printed < debug:
                    sys.stderr.write(f"[debug] header: {header}\n")
                    sys.stderr.write(f"[debug] candidates: {sorted(cands)}\n")
                    sys.stderr.write(f"[debug] matched: {hit}\n")
                    if not hit:
                        # show a few of the target samples
                        sample = ", ".join(list(sorted(targets))[:10])
                        sys.stderr.write(f"[debug] targets(sample): {sample}\n")
                    debug_printed += 1

                if hit:
                    out_header = cleaned_name_from_header(header) if clean else header
                    if clean and not out_header:
                        out_header = header
                    fout.write(f">{out_header}\n".encode("utf-8"))
                    for s in seq_lines:
                        fout.write(s.encode("utf-8"))
                    n_written += 1

            for line in inp:
                if line.startswith(">"):
                    flush()
                    n_seen += 1
                    header = line[1:].rstrip("\r\n")
                    seq_lines = []
                else:
                    if header is None: continue
                    seq_lines.append(line.rstrip("\r\n") + "\n")
            flush()
    finally:
        try: fin.close()
        except Exception: pass
        try: fout.close()
        except Exception: pass

    sys.stderr.write(f"[extract_fasta_simple] Records seen: {n_seen}\n")
    sys.stderr.write(f"[extract_fasta_simple] Records written: {n_written}\n")

def main():
    ap = argparse.ArgumentParser(description="Extract FASTA records by name list; optional header cleanup.")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA (.fa/.fasta or .gz). Use '-' for stdin.")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA (.fa/.fasta or .gz). Use '-' for stdout.")
    ap.add_argument("-l", "--list", required=True, help="Text file of targets (one per line).")
    ap.add_argument("--clean", action="store_true", help="Rewrite headers to compact names (e.g., 16, X, MT).")
    ap.add_argument("--debug", type=int, default=0, help="Print match diagnostics for up to N records.")
    args = ap.parse_args()
    run_extract(args.input, args.output, args.list, args.clean, args.debug)

if __name__ == "__main__":
    main()
