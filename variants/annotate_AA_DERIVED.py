#!/usr/bin/env python3
import pysam
import sys
import os

if len(sys.argv) != 3:
    sys.stderr.write(f"Usage: {os.path.basename(sys.argv[0])} <input.vcf.gz> <ancestral.fa.gz>\n")
    sys.exit(1)

vcf_path = sys.argv[1]
fasta_path = sys.argv[2]

# Open files
vcf_in = pysam.VariantFile(vcf_path)
fasta = pysam.FastaFile(fasta_path)

# Add INFO fields if not present
if "AA" not in vcf_in.header.info:
    vcf_in.header.info.add("AA", 1, "String", "Ancestral allele from FASTA")
if "DERIVED" not in vcf_in.header.info:
    vcf_in.header.info.add("DERIVED", ".", "String", "Derived allele(s) relative to ancestral")

# Output VCF
vcf_out = pysam.VariantFile("-", 'w', header=vcf_in.header)

# Count annotated sites
count = 0
for rec in vcf_in:
    try:
        base = fasta.fetch(rec.chrom, rec.pos - 1, rec.pos).upper()
        rec.info["AA"] = base

        # Determine derived alleles: any REF or ALT not matching AA
        derived = []
        if rec.ref != base:
            derived.append(rec.ref)
        for alt in rec.alts:
            if alt != base:
                derived.append(alt)
        if derived:
            rec.info["DERIVED"] = ",".join(derived)
        else:
            rec.info["DERIVED"] = "."

        count += 1
    except KeyError:
        # Chromosome not found in FASTA
        rec.info["AA"] = "."
        rec.info["DERIVED"] = "."

    vcf_out.write(rec)

sys.stderr.write(f"[INFO] Annotated {count} sites with AA and DERIVED\n")
