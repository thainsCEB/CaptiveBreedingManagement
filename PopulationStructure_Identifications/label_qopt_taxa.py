#!/usr/bin/env python3
import argparse, re, sys, subprocess
from pathlib import Path
from typing import List, Optional, Dict, Tuple

def tokenize_ws(s: str) -> List[str]:
    if s is None:
        return []
    s = s.replace("\u00A0", " ").replace("\u2007", " ").replace("\u202F", " ")
    s = re.sub(r"\s+", " ", s.strip())
    return s.split(" ") if s else []

def looks_like_header(tokens: List[str]) -> bool:
    if not tokens: return False
    nonnum = any(not re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and kstyle

def is_non_numeric_header(tokens: List[str]) -> bool:
    if not tokens: return False
    nonnum = any(not re.match(r"^[0-9.+-eE]+$", t) for t in tokens)
    kstyle = all(re.match(r"^K\d+$", t) for t in tokens)
    return nonnum and not kstyle

def is_k_header_line(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#"): return False
    return re.fullmatch(r'(?:\s*K\d+\s+)+K\d+\s*', s) is not None

def read_taxa_map(path: str) -> Tuple[List[str], Optional[Dict[str, str]], Optional[List[Tuple[str,str]]]]:
    labels: List[str] = []
    mapping: Dict[str, str] = {}
    multicol = False
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"): continue
            parts = tokenize_ws(s)
            if len(parts) == 1 and not multicol:
                labels.append(parts[0])
            elif len(parts) >= 2:
                multicol = True
                mapping[parts[0]] = parts[-1]
                labels.append(parts[-1])
    if multicol and not mapping: raise SystemExit(f"[error] No valid (old,new) pairs found in: {path}")
    if not multicol and not labels: raise SystemExit(f"[error] No taxa labels found in: {path}")
    return labels, (mapping if multicol else None), (list(mapping.items()) if multicol else None)

def get_taxa_labels(args) -> Tuple[List[str], Optional[Dict[str, str]], Optional[List[Tuple[str,str]]]]:
    if args.taxa: return args.taxa, None, None
    if args.taxa_map: return read_taxa_map(args.taxa_map)
    raise SystemExit("[error] Provide taxa via --taxa or --taxa-map")

def load_file_lines(p: Path) -> List[str]:
    return p.read_text(encoding="utf-8").splitlines()

def write_lines(p: Path, lines: List[str]) -> None:
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")

def _sed_escape_pattern(s: str) -> str:
    return s.replace("\\", "\\\\").replace("|", "\\|")

def _sed_escape_repl(s: str) -> str:
    return s.replace("\\", "\\\\").replace("&", "\\&").replace("|", "\\|")

def apply_sed_pairs_to_text(text: str, pairs) -> str:
    for old, new in pairs:
        pat = _sed_escape_pattern(old)
        rep = _sed_escape_repl(new)
        cp = subprocess.run(["sed", f"s|{pat}|{rep}|g"], input=text.encode("utf-8"),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if cp.returncode != 0:
            raise SystemExit(f"[error] sed failed for pair ({old}->{new}): {cp.stderr.decode('utf-8','ignore')}")
        text = cp.stdout.decode("utf-8", "strict")
    return text

def relabel_qopt_file(qopt: Path,
                      taxa_labels: List[str],
                      in_place: bool,
                      suffix: str,
                      header_policy: str = "taxa",
                      sep: str = " ",
                      taxa_map: Optional[Dict[str, str]] = None,
                      taxa_pairs: Optional[List[Tuple[str,str]]] = None) -> Path:
    lines = load_file_lines(qopt)

    # 1) Remove any K-header lines early
    lines = [ln for ln in lines if not looks_like_header(tokenize_ws(ln.strip())) and not is_k_header_line(ln)]

    # 2) Apply sed-style replacements across the entire file (always-on)
    if taxa_pairs:
        full_text = "\n".join(lines) + "\n"
        full_text = apply_sed_pairs_to_text(full_text, taxa_pairs)
        lines = full_text.splitlines()

    # 3) Remove K-headers again in case sed introduced/left one
    lines = [ln for ln in lines if not looks_like_header(tokenize_ws(ln.strip())) and not is_k_header_line(ln)]

    # 4) Detect existing header (non-numeric first non-comment line)
    first_idx: Optional[int] = None
    for idx, ln in enumerate(lines):
        s = ln.strip()
        if not s or s.startswith("#"): continue
        first_idx = idx; break

    existing_header_idx: Optional[int] = None
    existing_header_tokens: Optional[List[str]] = None
    if first_idx is not None:
        toks0 = tokenize_ws(lines[first_idx].strip())
        if is_non_numeric_header(toks0):
            existing_header_idx = first_idx
            existing_header_tokens = toks0

    # 5) Build replacement labels (rep_labels)
    rep_labels = list(taxa_labels)
    if taxa_map is not None:
        if existing_header_tokens:
            rep_labels = [taxa_map.get(tok, tok) for tok in existing_header_tokens]
        else:
            rep_labels = list(taxa_labels)

    # 6) Compose output: single header + body (drop any second header line)
    if header_policy == "none":
        out_lines = [ln for i, ln in enumerate(lines) if i != existing_header_idx]
    elif header_policy == "keep":
        out_lines = lines[:]
    else:
        body = [ln for i, ln in enumerate(lines) if i != existing_header_idx]
        out_lines = [sep.join(rep_labels)] + body
        # Trim any additional header line immediately following
        trimmed = []
        header_seen = False
        for ln in out_lines:
            s = ln.strip()
            if not s or s.startswith('#'):
                trimmed.append(ln); continue
            toks = tokenize_ws(s)
            if not header_seen:
                trimmed.append(ln); header_seen = True
            else:
                if is_k_header_line(ln) or is_non_numeric_header(toks):
                    continue
                trimmed.append(ln)
        out_lines = trimmed

    # 7) Warn if K mismatch
    K_guess: Optional[int] = None
    for ln in out_lines:
        s = ln.strip()
        if not s or s.startswith("#"): continue
        toks = tokenize_ws(s)
        if toks == rep_labels: continue
        nums = [t for t in toks if re.match(r"^[0-9.+-eE]+$", t)]
        if not nums: continue
        K_guess = len(nums) if K_guess is None else max(K_guess, len(nums))
    if K_guess is not None and header_policy == "taxa" and len(rep_labels) != K_guess:
        print(f"[warn] {qopt.name}: data has {K_guess} columns; taxa_map/taxa has {len(rep_labels)} labels",
              file=sys.stderr)

    outp = qopt if in_place else qopt.with_suffix(qopt.suffix + suffix)
    write_lines(outp, out_lines)
    return outp

def main():
    ap = argparse.ArgumentParser(description="Relabel .qopt headers (K1..KN) with taxa/population names.")
    ap.add_argument("--qopt-dir", required=False, help="Directory containing .qopt files")
    ap.add_argument("--file", required=False, help="Operate on a single .qopt file")
    ap.add_argument("--taxa-map", required=False, help="File with taxa labels (1-col or 2-col original->new)")
    ap.add_argument("--taxa", required=False, nargs="+", help="Inline taxa labels in order (space-separated)")
    ap.add_argument("--pattern", default="*.qopt", help="Glob pattern for qopt files (default: *.qopt)")
    ap.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    ap.add_argument("--in-place", action="store_true", help="Overwrite files in-place (default: write copies)")
    ap.add_argument("--suffix", default=".labeled", help="Suffix for output copies when not in-place (default: .labeled)")
    ap.add_argument("--sep", choices=["tab","space"], default="space", help="Separator for header (default: space)")
    ap.add_argument("--header", choices=["taxa","none","keep"], default="taxa",
                    help="Header policy: 'taxa' ensure first line is taxa labels; 'none' remove header; 'keep' preserve existing non-numeric header.")
    args = ap.parse_args()

    taxa_labels, taxa_mapping, taxa_pairs = get_taxa_labels(args)

    files: List[Path] = []
    if args.file:
        pth = Path(args.file)
        if not pth.exists(): raise SystemExit(f"[error] File not found: {pth}")
        files = [pth]
    else:
        if not args.qopt_dir: raise SystemExit("[error] Provide either --file or --qopt-dir")
        qdir = Path(args.qopt_dir)
        if not qdir.is_dir(): raise SystemExit(f"[error] Not a directory: {qdir}")
        files = list(qdir.rglob(args.pattern)) if args.recursive else list(qdir.glob(args.pattern))
        if not files: raise SystemExit(f"[error] No files matched pattern {args.pattern} in {qdir}")

    sep_char = "\t" if args.sep == "tab" else " "
    n_ok = 0
    for f in sorted(files):
        outp = relabel_qopt_file(
            f, taxa_labels,
            in_place=args.in_place,
            suffix=args.suffix,
            header_policy=args.header,
            sep=sep_char,
            taxa_map=taxa_mapping,
            taxa_pairs=taxa_pairs
        )
        action = "overwrote" if args.in_place else "wrote"
        print(f"[labeled] {action} {outp}")
        n_ok += 1
    print(f"[done] labeled {n_ok} file(s).")

if __name__ == "__main__":
    main()
