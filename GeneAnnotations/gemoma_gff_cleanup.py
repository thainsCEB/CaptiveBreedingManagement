#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
import os
import shutil
from collections import defaultdict
from datetime import datetime

def run_cmd(cmd, log_file=None, silent=False):
    print(f"Running: {cmd}")
    stderr = subprocess.DEVNULL if silent else None
    try:
        if log_file:
            with open(log_file, "a") as lf:
                subprocess.run(cmd, shell=True, stdout=lf, stderr=lf, check=True)
        else:
            subprocess.run(cmd, shell=True, stderr=stderr, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Command failed: {cmd}")
        if log_file:
            print(f"Check log file for details: {log_file}")
        else:
            print("No log file specified.")
        raise

def validate_gff(gff_file):
    hierarchy = {'mRNA': 'gene', 'exon': 'mRNA', 'CDS': 'mRNA'}
    features = {}
    parent_to_children = defaultdict(list)
    child_to_parent = {}

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] not in hierarchy and parts[2] != 'gene':
                continue
            attrs = dict([tuple(a.strip().split('=')) for a in parts[8].split(';') if '=' in a])
            feat_id = attrs.get('ID')
            if not feat_id:
                print(f"Skipping entry without ID: {line}")
                continue
            if parts[2] == 'CDS':
                feat_id += f".{parts[0]}:{parts[3]}-{parts[4]}"
            parent_id = attrs.get('Parent')
            features[feat_id] = {
                'type': parts[2], 'lend': int(parts[3]), 'rend': int(parts[4]),
                'orient': parts[6], 'contig': parts[0], 'parent': parent_id
            }
            if parent_id:
                parent_to_children[parent_id].append(feat_id)
                child_to_parent[feat_id] = parent_id

    for parent, children in parent_to_children.items():
        if parent not in features:
            print(f"Missing parent feature: {parent}")
            continue
        pfeat = features[parent]
        for child in children:
            cfeat = features[child]
            if pfeat['contig'] != cfeat['contig']:
                print(f"Different contig: {parent} {child}")
            if pfeat['orient'] != cfeat['orient']:
                print(f"Orientation mismatch: {parent} {child}")
            if not (pfeat['lend'] <= cfeat['lend'] and pfeat['rend'] >= cfeat['rend']):
                print(f"Not nested: {parent} {child}")
            if hierarchy.get(cfeat['type']) != pfeat['type']:
                print(f"Bad hierarchy: {parent} {child}")

    for fid, feat in features.items():
        if feat['type'] == 'gene':
            continue
        parent_id = feat['parent']
        if not parent_id:
            print(f"No parent for {fid}")
        elif parent_id in features:
            if hierarchy.get(feat['type']) != features[parent_id]['type']:
                print(f"Bad parent type for {fid}")

    for fid, feat in features.items():
        if feat['type'] != 'mRNA':
            continue
        children = parent_to_children.get(fid, [])
        exons = [features[cid] for cid in children if features[cid]['type'] == 'exon']
        cds = [features[cid] for cid in children if features[cid]['type'] == 'CDS']
        if not cds:
            print(f"mRNA {fid} has no CDS")
        used = set()
        for cd in cds:
            found = False
            for ex in exons:
                if cd['lend'] >= ex['lend'] and cd['rend'] <= ex['rend']:
                    if (ex['lend'], ex['rend']) in used:
                        print(f"CDS {cd} maps to multiple exons")
                    used.add((ex['lend'], ex['rend']))
                    found = True
                    break
            if not found:
                print(f"CDS {cd} not encapsulated")

def main():
    parser = argparse.ArgumentParser(description="GFF3 Cleanup Pipeline with Canonical Isoform Annotation")
    parser.add_argument("--gff", "-g", required=True, help="Input GFF3 file")
    parser.add_argument("--fasta", "-f", required=True, help="Reference genome FASTA")
    parser.add_argument("--species", "-s", required=True, help="Species prefix for IDs")
    parser.add_argument("--prefix", "-pfx", required=False, help="Prefix for all output files (overrides species in filenames)")
    parser.add_argument("--outdir", "-o", required=True, help="Output directory")
    parser.add_argument("--cleanup", "-c", action="store_true", help="Delete temp directory after run")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = os.path.join(args.outdir, "tmpGFF")
    os.makedirs(tmpdir, exist_ok=True)

    agat_log = os.path.join(args.outdir, "agat_id_log.txt")
    step1 = os.path.join(tmpdir, "step1.clean.gff")
    run_cmd(f"gff_cleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= {args.gff} -o {step1}")

    step2 = os.path.join(tmpdir, "step2.gff")
    run_cmd(f"gt gff3 -sort -tidy -retainids -addintrons {step1} > {step2}")

    cleaned_final = os.path.join(tmpdir, "final.clean.gff")
    run_cmd(f"cp {step2} {cleaned_final}")

    today = datetime.today().strftime('%m%d%y')
    prefix = args.prefix if args.prefix else args.species
    final_gff = os.path.join(args.outdir, f"{prefix}.final.cleaned.{today}.gff")
    run_cmd(f"cp {cleaned_final} {final_gff}")

    print("Running GFF3 structure validator...")
    validate_gff(final_gff)

    run_cmd(f"agat_sp_manage_IDs.pl -f {final_gff} --ensembl --prefix {args.species} -o {final_gff}.tmp", log_file=agat_log)
    run_cmd(f"mv {final_gff}.tmp {final_gff} && sed -i 's/{args.species}M/{args.species}T/g;s/AGAT/GeMoMa/g' {final_gff}")

    pep = os.path.join(args.outdir, f"{prefix}.proteins.faa")
    cds = os.path.join(args.outdir, f"{prefix}.cds.fa")
    tx = os.path.join(args.outdir, f"{prefix}.transcripts.fa")
    run_cmd(f"gffread {final_gff} -g {args.fasta} -y {pep} -x {cds} -w {tx}")

    # Rewrite protein FASTA headers to include both protein ID and transcript ID
    print("Rewriting protein FASTA headers with protein and transcript IDs...")
    updated_proteins = []
    with open(pep) as fin:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                tid = header.split()[0]
                pid = tid.replace("T", "P") if "T" in tid else f"{tid}_P"
                updated_proteins.append(f">{pid} transcript_id={tid}")
            else:
                updated_proteins.append(line)
    with open(pep, "w") as fout:
        fout.writelines(updated_proteins)

    canonical_gff = os.path.join(args.outdir, f"{prefix}.canonical.gff")
    run_cmd(f"agat_sp_keep_longest_isoform.pl -gff {final_gff} -o {canonical_gff}", log_file=agat_log)

    canonical_pep = os.path.join(args.outdir, f"{prefix}.canonical.proteins.fa")
    canonical_tx = os.path.join(args.outdir, f"{prefix}.canonical.transcripts.fa")
    run_cmd(f"gffread {canonical_gff} -g {args.fasta} -y {canonical_pep} -w {canonical_tx}")

    print("Annotating canonical isoforms in GFF...")
    canonical_ids = set()
    with open(canonical_gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] not in {"mRNA", "transcript"}:
                continue
            attrs = dict([tuple(a.strip().split('=')) for a in parts[8].split(';') if '=' in a])
            tid = attrs.get("ID")
            if tid:
                canonical_ids.add(tid)

    updated_lines = []
    with open(final_gff) as fin:
        for line in fin:
            if line.startswith("#") or "\t" not in line:
                updated_lines.append(line)
                continue
            parts = line.strip().split("\t")
            if parts[2] in {"mRNA", "transcript"}:
                attrs = parts[8]
                kv_pairs = [a.strip() for a in attrs.split(';') if '=' in a]
                attr_dict = dict([a.split('=') for a in kv_pairs])
                tid = attr_dict.get("ID")
                if tid and tid in canonical_ids:
                    attr_dict["canonical"] = "1"
                    parts[8] = ";".join([f"{k}={v}" for k, v in attr_dict.items()])
            updated_lines.append("\t".join(parts) + "\n")

    with open(final_gff, "w") as fout:
        fout.writelines(updated_lines)

    print(f"Canonical isoform annotations written to: {final_gff}")

    if args.cleanup:
        shutil.rmtree(tmpdir)
        print("Temporary files cleaned up.")

if __name__ == "__main__":
    main()
