#!/usr/bin/env python3

import os
import subprocess
import argparse
import shutil
import urllib.request

def check_and_download_lua(script="filter_protein_match.lua"):
    if not os.path.exists(script):
        print(f"[INFO] {script} not found. Downloading...")
        url = "https://raw.githubusercontent.com/satta/filterama/master/filter_protein_match.lua"
        urllib.request.urlretrieve(url, script)
        print(f"[INFO] Downloaded {script}")

def run(cmd, outfile=None):
    if outfile and os.path.exists(outfile):
        if os.path.getsize(outfile) > 0:
            print(f"[✔] Skipping: {outfile} already exists and is non-empty.")
            return
        else:
            print(f"[!] Re-running: {outfile} exists but is empty.")
    else:
        print(f"[→] Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Repeat annotation pipeline with checkpointing and cleanup")
    parser.add_argument("-a", "--assembly", required=True, help="Input genome FASTA")
    parser.add_argument("-r", "--repeatmodeler", required=True, help="RepeatModeler consensi.fa.classified")
    parser.add_argument("-e", "--edta", required=True, help="EDTA TE FASTA file")
    parser.add_argument("-s", "--sines", required=True, help="SINEs.bnk file")
    parser.add_argument("-m", "--hmmdir", required=True, help="Directory of LTRdigest HMMs")
    parser.add_argument("-u", "--uniprot", required=True, help="Uniprot filtered reviewed FASTA")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for final output library")
    parser.add_argument("-t", "--threads", type=int, default=10, help="Number of threads")
    parser.add_argument("-c", "--cleanup", action="store_true", help="Delete intermediate files after completion")
    args = parser.parse_args()

    check_and_download_lua()

    assembly = os.path.abspath(args.assembly)
    base = os.path.splitext(os.path.basename(assembly))[0]
    prefix = args.prefix
    final_output = f"{prefix}.Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa"
    intermediates = []

    # TransposonPSI
    tpsi_base = f"{assembly}.TPSI.allHits.chains.bestPerLocus"
    tpsi_fa = f"{tpsi_base}.fa"
    if os.path.exists(tpsi_fa) and os.path.getsize(tpsi_fa) > 0:
        print(f"[✔] Skipping TransposonPSI: {tpsi_fa} already exists and is non-empty.")
    else:
        run(f"transposonPSI.pl {assembly} nuc")
        run(f"transposonPSI_2fasta.pl {assembly} {tpsi_base}", outfile=tpsi_fa)
    intermediates.append(tpsi_fa)

    # LTRharvest/LTRdigest
    ltr_fasta = f"{base}.ltrh.sorted.ltrd.filtered.sequences.fasta"
    if os.path.exists(ltr_fasta) and os.path.getsize(ltr_fasta) > 0:
        print(f"[✔] Skipping LTRharvest and LTRdigest: {ltr_fasta} already exists and is non-empty.")
    else:
        ltrh_gff3 = f"{base}.ltrh.gff3"
        run(f"gt suffixerator -dna -db {assembly} -lcp -ssp -suf -tis -des -lossless")
        run(f"gt ltrharvest -index {assembly} -tabout no -seqids -md5 > {ltrh_gff3}", outfile=ltrh_gff3)
        run(f"gt stat {ltrh_gff3} > stat.{base}.gff3.out")
        sorted_gff3 = f"{base}.ltrh.sorted.gff3"
        run(f"gt gff3 -sort {ltrh_gff3} > {sorted_gff3}", outfile=sorted_gff3)
        ltrd_gff3 = f"{base}.ltrh.sorted.ltrd.gff3"
        run(f"gt -j {args.threads} ltrdigest -hmms {args.hmmdir}/*hmm -encseq {assembly} -matchdescstart < {sorted_gff3} > {ltrd_gff3}", outfile=ltrd_gff3)
        filtered_gff3 = f"{base}.ltrh.sorted.ltrd.filtered.gff3"
        run(f"gt select -rule_files filter_protein_match.lua < {ltrd_gff3} > {filtered_gff3}", outfile=filtered_gff3)
        run(f"gt stat {filtered_gff3} >> stat.{base}.gff3.out")
        run(f"gt extractfeat -type LTR_retrotransposon -encseq {assembly} -matchdescstart < {filtered_gff3} > {ltr_fasta}", outfile=ltr_fasta)
    intermediates.append(ltr_fasta)

    # Combine + seqtk + usearch
    combined_fa = "Combined.rep.library.fasta"
    minlen_fa = "Combined.rep.library.minlen50.fasta"
    clustered_fa = "Combined.rep.library.minlen50.cdhit.fasta"

    if os.path.exists(clustered_fa) and os.path.getsize(clustered_fa) > 0:
        print(f"[✔] Skipping sequence filtering and clustering: {clustered_fa} already exists and is non-empty.")
    else:
        if os.path.exists(combined_fa) and os.path.getsize(combined_fa) > 0 and \
           os.path.exists(tpsi_fa) and os.path.getsize(tpsi_fa) > 0 and \
           os.path.exists(ltr_fasta) and os.path.getsize(ltr_fasta) > 0:
            print(f"[✔] Skipping combination: {combined_fa} already exists and all inputs are complete.")
        else:
            run(f"cat {args.repeatmodeler} {args.edta} {ltr_fasta} {tpsi_fa} {args.sines} > {combined_fa}", outfile=combined_fa)
        intermediates.append(combined_fa)

        run(f"seqtk seq -L 50 {combined_fa} > {minlen_fa}", outfile=minlen_fa)
        intermediates.append(minlen_fa)

        run(f"usearch -cluster_fast {minlen_fa} -id 0.8 -consout {clustered_fa}", outfile=clustered_fa)
    intermediates.append(clustered_fa)

    run(f"RepeatClassifier -consensi {clustered_fa}")
    classified_fa = f"{clustered_fa}.classified"
    intermediates.append(classified_fa)

    # Filter unknowns
    unknown_ids = "UnknownIDs.txt"
    run(f"grep 'Unknown' {classified_fa} > {unknown_ids}", outfile=unknown_ids)
    run(f"sed -i 's/>//g' {unknown_ids}")
    run(f"seqtk subseq {classified_fa} {unknown_ids} > Unknown_repeats.fasta", outfile="Unknown_repeats.fasta")
    intermediates += [unknown_ids, "Unknown_repeats.fasta"]

    run(f"makeblastdb -in {args.uniprot} -dbtype prot")
    run(f"blastx -query Unknown_repeats.fasta -db {args.uniprot} -evalue 1e-10 -num_threads {args.threads} -max_target_seqs 1 "
        f"-outfmt '6 qseqid sseqid evalue bitscore sgi sacc stitle' -out Blast_out.txt", outfile="Blast_out.txt")

    run("awk -F '\t' '{print $1,$7}' Blast_out.txt | sort | uniq "
        "| grep -i -v -E 'transposon|Copia protein|mobile element|transposable|transposase' "
        "| awk '{print $1}' > Unknowns_with_Port_hit.txt", outfile="Unknowns_with_Port_hit.txt")

    run(f"faSomeRecords -exclude {classified_fa} Unknowns_with_Port_hit.txt {final_output}", outfile=final_output)

    print(f"\n✅ Completed! Final repeat library: {final_output}")

    if args.cleanup:
        print("\n🧹 Cleaning up intermediate files...")
        for f in intermediates + ["Blast_out.txt", "Unknowns_with_Port_hit.txt"]:
            try:
                os.remove(f)
                print(f"Removed: {f}")
            except FileNotFoundError:
                continue

if __name__ == "__main__":
    main()
