#!/usr/bin/env python3
import argparse, os, shutil, subprocess, re
from datetime import datetime

# ---------- Step registry ----------
STEPS = [
    "step1",           # gff_cleaner -> step1.clean.gff
    "step2",           # gt gff3 sort/tidy -> step2.gff
    "copy_cleaned",    # copy to final.clean.gff
    "copy_final",      # copy to dated final.gff
    "manage_ids",      # agat_sp_manage_IDs.pl
    "sed_rename",      # sed replacements (AGAT/AnnotationFinalizer_/species M->T)
    "add_pid",         # add protein_id to mRNA/CDS
    "gffread",         # gffread -> pep/cds/tx
    "norm_fasta",      # normalize FASTA headers
    "keep_longest",    # agat_sp_keep_longest_isoform.pl
    "canonical_flags"  # annotate canonical=1 back into final_gff
]
STEP_HELP = "\n".join(f"  - {s}" for s in STEPS)

# ---------- Utilities ----------
def run_cmd(cmd, log_file=None):
    print(f"[CMD] {cmd}")
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as p:
        out, _ = p.communicate()
        if log_file:
            with open(log_file, "a") as lf:
                lf.write(out or "")
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, cmd)

def _state_path(outdir):
    return os.path.join(outdir, ".gemoma_cleanup.state.json")

def _state_load(outdir):
    sp = _state_path(outdir)
    if os.path.exists(sp):
        import json
        with open(sp) as f:
            try:
                return json.load(f)
            except Exception:
                return {}
    return {}

def _state_save(outdir, d):
    sp = _state_path(outdir)
    import json
    with open(sp, "w") as f:
        json.dump(d, f, indent=2)

def _mark_done(state, key, value=True):
    state[key] = bool(value)
    return state

def _is_done(state, key):
    return state.get(key, False)

# ---------- Feature helpers ----------
def _normalize_fasta_headers(fpath, convert_tid_to_pid=False):
    updated = []
    with open(fpath) as fin:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                header = header.replace("AnnotationFinalizer_", "GeMoMa")
                header = re.sub(r"\bAGAT\b", "GeMoMa", header)
                header = header.replace("AGAT_", "GeMoMa_")
                if convert_tid_to_pid:
                    tid = header.split()[0]
                    pid = tid.replace("T", "P") if "T" in tid else f"{tid}_P"
                    header = f"{pid} transcript_id={tid}"
                updated.append(f">{header}\n")
            else:
                updated.append(line)
    with open(fpath, "w") as fout:
        fout.writelines(updated)

def _add_protein_id_to_gff(gff_path):
    def parse_attrs(attr_str):
        attrs = {}
        for part in filter(None, [x.strip() for x in attr_str.strip().split(";") ]):
            if "=" in part:
                k, v = part.split("=", 1)
                attrs[k] = v
            else:
                attrs[part] = ""
        return attrs

    def attrs_to_str(attrs):
        common = ["ID", "Parent", "Name", "gene_id", "transcript_id", "protein_id", "canonical"]
        parts = []
        for k in common:
            if k in attrs:
                parts.append(f"{k}={attrs[k]}")
        for k in sorted(k for k in attrs.keys() if k not in set(common)):
            parts.append(f"{k}={attrs[k]}")
        return ";".join(parts)

    # pass 1: build mRNA/transcript -> protein_id
    mrna_to_pid = {}
    lines = []
    with open(gff_path) as fin:
        for line in fin:
            lines.append(line)
            if line.startswith("#") or "\t" not in line: continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            ftype = cols[2]
            attrs = parse_attrs(cols[8])
            if ftype in ("mRNA", "transcript"):
                tid = attrs.get("ID") or attrs.get("transcript_id")
                if tid:
                    pid = tid.replace("T", "P") if "T" in tid else f"{tid}_P"
                    mrna_to_pid[tid] = pid

    # pass 2: write with protein_id on mRNA and CDS
    out_lines = []
    for line in lines:
        if line.startswith("#") or "\t" not in line:
            out_lines.append(line); continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            out_lines.append(line); continue
        ftype = cols[2]
        attrs = parse_attrs(cols[8])

        if ftype in ("mRNA", "transcript"):
            tid = attrs.get("ID") or attrs.get("transcript_id")
            if tid:
                pid = mrna_to_pid.get(tid) or (tid.replace("T", "P") if "T" in tid else f"{tid}_P")
                attrs["protein_id"] = pid
                if "transcript_id" not in attrs and "ID" in attrs:
                    attrs["transcript_id"] = attrs["ID"]
            cols[8] = attrs_to_str(attrs)
            out_lines.append("\t".join(cols) + "\n")
            continue

        if ftype == "CDS":
            parents = []
            if "Parent" in attrs:
                parents = [p.strip() for p in attrs["Parent"].split(",") if p.strip()]
            pid = None
            for par in parents:
                if par in mrna_to_pid:
                    pid = mrna_to_pid[par]; break
            if not pid:
                tid = attrs.get("transcript_id") or (parents[0] if parents else None)
                if tid:
                    pid = tid.replace("T", "P") if ("T" in tid) else f"{tid}_P"
            if pid:
                attrs["protein_id"] = pid
            cols[8] = attrs_to_str(attrs)
            out_lines.append("\t".join(cols) + "\n")
            continue

        out_lines.append(line)

    with open(gff_path, "w") as fout:
        fout.writelines(out_lines)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(
        description="GFF3 Cleanup Pipeline with resume and named-step support\n\nSteps:\n" + STEP_HELP,
        epilog="Use --start-at to begin at a specific step; --list-steps to print them."
    )
    ap.add_argument("--gff", "-g", required=True)
    ap.add_argument("--fasta", "-f", required=True)
    ap.add_argument("--species", "-s", required=True)
    ap.add_argument("--prefix", "-pfx", required=False, help="Prefix for outputs (defaults to species)")
    ap.add_argument("--outdir", "-o", required=True)
    ap.add_argument("--cleanup", "-c", action="store_true")
    ap.add_argument("--resume", action="store_true", help="Resume by skipping finished steps (tracked in .gemoma_cleanup.state.json)")
    ap.add_argument("--start-at", choices=STEPS, default=STEPS[0], help="Start (or resume) at a named step.")
    ap.add_argument("--list-steps", action="store_true", help="List step names and exit")
    args = ap.parse_args()

    if args.list_steps:
        print("Steps:\n" + STEP_HELP)
        return

    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = os.path.join(args.outdir, "tmpGFF")
    os.makedirs(tmpdir, exist_ok=True)

    state = _state_load(args.outdir) if args.resume else {}

    # helpers for gating
    start_idx = STEPS.index(args.start_at)
    def should_run(step):
        return (STEPS.index(step) >= start_idx) and (not args.resume or not _is_done(state, step))
    def mark_done(step):
        nonlocal state
        state = _mark_done(state, step); _state_save(args.outdir, state)

    agat_log = os.path.join(args.outdir, "agat_id_log.txt")
    step1 = os.path.join(tmpdir, "step1.clean.gff")
    step2 = os.path.join(tmpdir, "step2.gff")
    cleaned_final = os.path.join(tmpdir, "final.clean.gff")
    today = datetime.today().strftime("%m%d%y")
    prefix = args.prefix if args.prefix else args.species
    final_gff = os.path.join(args.outdir, f"{prefix}.final.cleaned.{today}.gff")
    pep = os.path.join(args.outdir, f"{prefix}.proteins.faa")
    cds = os.path.join(args.outdir, f"{prefix}.cds.fa")
    tx  = os.path.join(args.outdir, f"{prefix}.transcripts.fa")
    canonical_gff = os.path.join(args.outdir, f"{prefix}.canonical.gff")
    canonical_pep = os.path.join(args.outdir, f"{prefix}.canonical.proteins.fa")
    canonical_tx  = os.path.join(args.outdir, f"{prefix}.canonical.transcripts.fa")

    # step1: clean
    if should_run("step1"):
        run_cmd(f"gff_cleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates {args.gff} -o {step1}")
        mark_done("step1")
    else:
        print("[resume] Skipping step1")

    # step2: sort/tidy
    if should_run("step2"):
        run_cmd(f"gt gff3 -sort -tidy -retainids -addintrons {step1} > {step2}")
        mark_done("step2")
    else:
        print("[resume] Skipping step2")

    # copy to cleaned_final
    if should_run("copy_cleaned"):
        shutil.copy2(step2, cleaned_final)
        mark_done("copy_cleaned")
    else:
        print("[resume] Skipping copy_cleaned")

    # copy to dated final_gff
    if should_run("copy_final"):
        shutil.copy2(cleaned_final, final_gff)
        mark_done("copy_final")
    else:
        print("[resume] Skipping copy_final")

    # manage IDs
    if should_run("manage_ids"):
        run_cmd(f"agat_sp_manage_IDs.pl -f {final_gff} --ensembl --prefix {args.species} -o {final_gff}.tmp", log_file=agat_log)
        mark_done("manage_ids")
    else:
        print("[resume] Skipping manage_ids")

    # sed rename
    if should_run("sed_rename"):
        run_cmd(f"mv {final_gff}.tmp {final_gff} && sed -i 's/{args.species}M/{args.species}T/g;s/AGAT/GeMoMa/g;s/AnnotationFinalizer_/GeMoMa/g' {final_gff}")
        mark_done("sed_rename")
    else:
        print("[resume] Skipping sed_rename")

    # add protein_id to mRNA/CDS
    if should_run("add_pid"):
        _add_protein_id_to_gff(final_gff)
        mark_done("add_pid")
    else:
        print("[resume] Skipping add_pid")

    # gffread to FASTAs
    if should_run("gffread"):
        run_cmd(f"gffread {final_gff} -g {args.fasta} -y {pep} -x {cds} -w {tx}")
        mark_done("gffread")
    else:
        print("[resume] Skipping gffread")

    # normalize headers
    if should_run("norm_fasta"):
        _normalize_fasta_headers(pep, convert_tid_to_pid=True)
        _normalize_fasta_headers(cds)
        _normalize_fasta_headers(tx)
        mark_done("norm_fasta")
    else:
        print("[resume] Skipping norm_fasta")

    # keep longest isoform (canonical) and annotate back
    if should_run("keep_longest"):
        run_cmd(f"agat_sp_keep_longest_isoform.pl -gff {final_gff} -o {canonical_gff}", log_file=agat_log)
        mark_done("keep_longest")
    else:
        print("[resume] Skipping keep_longest")

    # extract IDs and annotate canonical=1 on transcripts
    if should_run("canonical_flags"):
        canonical_ids = set()
        if os.path.exists(canonical_gff):
            with open(canonical_gff) as fin:
                for line in fin:
                    if line.startswith("#") or "\t" not in line: continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 9 and parts[2] in {"mRNA", "transcript"}:
                        attrs = parts[8]
                        for kv in attrs.split(";"):
                            if kv.startswith("ID="):
                                canonical_ids.add(kv.split("=",1)[1])
                                break
        # annotate
        updated = []
        with open(final_gff) as fin:
            for line in fin:
                if line.startswith("#") or "\t" not in line:
                    updated.append(line); continue
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] in {"mRNA","transcript"}:
                    attrs = parts[8]
                    kv = [a for a in attrs.split(";") if "=" in a]
                    d = dict(a.split("=",1) for a in kv)
                    tid = d.get("ID")
                    if tid in canonical_ids:
                        d["canonical"] = "1"
                        parts[8] = ";".join([f"{k}={v}" for k,v in d.items()])
                        updated.append("\t".join(parts) + "\n")
                        continue
                updated.append(line)
        with open(final_gff, "w") as fo:
            fo.writelines(updated)
        mark_done("canonical_flags")
    else:
        print("[resume] Skipping canonical_flags")

    # optional cleanup
    if args.cleanup:
        shutil.rmtree(tmpdir, ignore_errors=True)
        print("[cleanup] Removed tmp dir.")

if __name__ == '__main__':
    main()
