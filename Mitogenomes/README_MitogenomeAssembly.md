# 🧬 Mitogenome Assembly and Evaluation Workflow

**Audience:** Wet-lab & bioinformatics collaborators working on mtDNA assemblies  
**Maintainer:** Taylor Hains  
**Last updated:** 2025-10-23

This README is a step-by-step playbook for assembling and validating mitochondrial genomes (mitogenomes) from WGS / capture datasets. It covers:

1. Estimating **mitochondrial read proportion (MDR)**  
2. Running **MitoZ** via the batch wrapper **`mitoForger.py`** (recommended when MDR is healthy)  
3. Reconstructing mitogenomes via **ANGSD** using **`angsd_mt_pipeline.py`** (recommended for low-coverage / low-MDR)  
4. QC, annotation, and curation steps

---

## Project layout (recommended)

```
project/
├── 00_raw/                   # raw reads, or links to reads
├── 01_mdr/                   # MDR estimation outputs
├── 02_mitoz/                 # MitoZ assemblies (mitoForger.py)
├── 03_angsd/                 # ANGSD consensus assemblies (angsd_mt_pipeline.py)
├── 04_qc/                    # coverage plots, summary tables
├── references/
│   ├── mt_reference.fasta    # mitochondrial reference (close relative)
│   └── mt_reference.fasta.fai
├── bed/                      # BED file for MT genes (for ANGSD pipeline)
│   └── mt_genes.bed
├── scripts/                  # mitoForger.py, assess_mdr.py, angsd_mt_pipeline.py
└── README_MitogenomeAssembly.md
```

---

## 1) Assessing Mitochondrial Read Proportion (MDR)

**Why:** MDR is the single best predictor of success for de novo mitogenome assembly.

| Outcome | MDR (fraction of reads) | Interpretation |
|--------:|:-------------------------|----------------|
| Excellent | ≥ 0.005 (0.5%) | Full circular assembly very likely |
| Adequate | 0.001–0.005 (0.1–0.5%) | Often assembles; minor gaps possible |
| Marginal | 0.0005–0.001 (0.05–0.1%) | Partial assembly common (gaps, missing tRNAs) |
| Poor | < 0.0005 (<0.05%) | Use ANGSD or enrichment-first strategies |

### Script: `assess_mdr.py`

- **BAM mode:** counts reads on mt contigs via `samtools idxstats`  
- **FASTQ mode:** aligns a subset of reads to a mt reference (minimap2 or bwa-mem2) and counts mapped reads

**Examples**
```bash
# BAM mode
python scripts/assess_mdr.py -b sample1.bam sample2.bam -o 01_mdr/MDR_summary.tsv

# FASTQ mode
python scripts/assess_mdr.py   -s sample_sheet.tsv   -m references/mt_reference.fasta   -n 2000000   -o 01_mdr/MDR_summary.tsv
```
`sample_sheet.tsv` (tab-delimited):
```
SampleID	R1.fq.gz	R2.fq.gz
```

**Low-MDR rule of thumb:** If MDR < 0.1%, prefer ANGSD-based consensus; if ≥ 0.1%, try MitoZ first.

---

## 2) Assembling mitogenomes with MitoZ via `mitoForger.py`

**Use when:** MDR is ≥ 0.1% and/or you want a full de novo assembly and automatic annotation.

**Sample sheet format (TSV):**
```
sampleID	read1	read2	speciesName
```
- `speciesName` must be **scientific name** (e.g., `Ara glaucogularis`), passed verbatim to MitoZ `--species`.

**Example run**
```bash
python scripts/mitoForger.py   -s sample_sheet.tsv   -o 02_mitoz   -t 16   -g 2   -p 1.0   -S /path/to/MitoZ_v3.6.sif
```
- `-p 1.0` = no subsampling. If you *must* subsample, use a fraction (e.g., `-p 0.5`) not a percent.
- Containerized runs: provide `-S` with a Singularity `.sif`; the script auto-binds inputs/outputs.

**Outputs (per sample)** — placed in `02_mitoz/mitogenomes/`:
- `<sample>.mitogenome.fasta` (header normalized to `>MT mitochondrion`)
- `<sample>.mito.genes.pcg.fasta` (protein-coding genes; simplified headers)
- `<sample>.mito.genes.all.fasta` (all annotated genes; simplified headers)
- `<sample>.gbf` (GenBank format)
- Logs: `02_mitoz/batch_mitoz.log`, and a global `mito.log` collecting “Potential missing genes:”

**Coverage & MDR notes (embedded in script):**
- MitoZ performs best when **MDR ≥ 0.1%**, ideally **0.3–0.5%**.
- If **MDR < 0.05%** or mt coverage < 30–50×, assemblies may be fragmented/missing tRNAs.
- Avoid overly aggressive subsampling when uncertain.

---

## 3) When to use ANGSD instead of MitoZ

Prefer **ANGSD** when:
- MDR < 0.1% or read depth is low/uneven
- Samples are **ancient/museum** with short, damaged fragments
- You already have BAMs from population pipelines and want a robust consensus quickly

ANGSD creates a **posterior consensus FASTA** from mapped reads (no de novo step), which is more stable at low coverage.

---

## 4) Assembling mitogenomes with ANGSD via `angsd_mt_pipeline.py`

**Correct flags (from `-h/--help`):**
```
-h, --help                      show this help message and exit
-s SAMPLE_SHEET, --sample-sheet SAMPLE_SHEET
                                Two-column TSV: sample_name<TAB>bam_path
-r REFERENCE, --reference REFERENCE
                                Reference genome FASTA
-b BED, --bed BED               BED file of MT genes
-m MT_REGION, --mt-region MT_REGION
                                MT region name/range (e.g., MT or chrM or 'MT:1-16569')
-t THREADS, --threads THREADS   Threads for ANGSD [4]
-o OUTDIR, --outdir OUTDIR      Output directory (will contain 'mitogenome/' and temp dirs)
```

### Inputs
- `--sample-sheet`: two-column TSV with **sample_name** and **bam_path** (absolute or relative)
- `--reference`: mitochondrial **or** whole-reference FASTA (ensure MT contig name matches `--mt-region`)
- `--bed`: BED file listing MT genes (used for per-gene stats/plots in the pipeline if implemented)
- `--mt-region`: region/contig name or explicit span (e.g., `MT`, `chrM`, or `MT:1-16569`)
- `--threads`: CPU threads for ANGSD steps
- `--outdir`: destination directory; pipeline creates `mitogenome/` and any temp/work subfolders

### Example
```bash
python scripts/angsd_mt_pipeline.py   -s 03_angsd/sample_bams.tsv   -r references/mt_reference.fasta   -b bed/mt_genes.bed   -m MT:1-16569   -t 12   -o 03_angsd
```

**Notes**
- Ensure the reference FASTA and `--mt-region` use the **same contig naming** (e.g., `MT` vs `chrM`).  
- If using a whole-genome reference, include only one mitochondrial contig in `--mt-region` (or give the span).  
- The pipeline writes per-sample consensus FASTAs under `OUTDIR/mitogenome/` (and logs/work files under temp dirs).  
- For annotation, you can optionally run `mitoz annotate` on each consensus FASTA afterward (see Section 6).

---

## 5) Low-coverage / Low-MDR strategies

If MitoZ fails or mtDNA is incomplete:

**A. Enrich mt reads before assembly**
```bash
minimap2 -ax sr references/mt_reference.fasta R1.fq.gz R2.fq.gz |   samtools view -b -F 4 - | samtools sort -o 04_qc/mito_enriched.bam
samtools fastq 04_qc/mito_enriched.bam -1 mito_R1.fq.gz -2 mito_R2.fq.gz
# Re-run mitoForger.py on mito_R1/mito_R2
```

**B. Reference-guided tools**
- **NOVOPlasty**
  ```bash
  perl NOVOPlasty.pl -c config.txt
  ```
- **GetOrganelle**
  ```bash
  get_organelle_from_reads.py -1 mito_R1.fq.gz -2 mito_R2.fq.gz -R 10 -k 21,45,65,85,105 -F animal_mt -o getorganelle_mt
  ```

**C. Iterative map/extend** (very low depth, ancient DNA)
- Repeatedly map → extract → remap to extend covered regions until plateau

---

## 6) QC, curation, and re-annotation

| Task | Commands / Tools | What to expect |
|-----:|-------------------|----------------|
| Coverage depth | `samtools depth sample.bam | awk '{sum+=$3} END {print sum/NR}'` | Mean depth ~30×+ preferred |
| Circularization | `circlator` or visual inspection | One contiguous circular contig |
| Re-annotation | `mitoz annotate` or **MITOS2** | ~37 genes: 13 PCGs, 22 tRNAs, 2 rRNAs |
| Concordance | BLAST vs close references | Similar length/order; no large gaps |
| Graphs | `Bandage`, `IGV` | Spot misassemblies/NUMTs |

---

## 7) Decision tree (who runs what)

1. **Run MDR**: `assess_mdr.py` → If **MDR ≥ 0.1%**, go to 2; else go to 3.  
2. **MitoZ path**: `mitoForger.py` → Review results; if fragmented, try enrichment or ANGSD.  
3. **ANGSD path**: `angsd_mt_pipeline.py` → Re-annotate FASTAs; if short/incomplete, consider NOVOPlasty/GetOrganelle.

---

## 8) Dependencies

- `samtools ≥ 1.15`
- `minimap2 ≥ 2.26` (or `bwa-mem2 ≥ 2.2`)
- `MitoZ ≥ 3.6` (+ Singularity `.sif` if containerized)
- `ANGSD ≥ 0.940`
- Optional: `NOVOPlasty`, `GetOrganelle`, `circlator`, `MITOS2`

---

## 9) Quick copy-paste block for new users

```bash
# 0) Make folders
mkdir -p 01_mdr 02_mitoz 03_angsd 04_qc references bed scripts

# 1) Assess MDR
python scripts/assess_mdr.py -s sample_sheet.tsv -m references/mt_reference.fasta -n 2000000 -o 01_mdr/MDR_summary.tsv

# 2) If MDR ≥ 0.1%, try MitoZ first
python scripts/mitoForger.py -s sample_sheet.tsv -o 02_mitoz -t 16 -g 2 -p 1.0 -S /path/MitoZ_v3.6.sif

# 3) If MDR < 0.1% (or MitoZ incomplete), use ANGSD pipeline (correct flags)
python scripts/angsd_mt_pipeline.py -s 03_angsd/sample_bams.tsv -r references/mt_reference.fasta -b bed/mt_genes.bed -m MT:1-16569 -t 12 -o 03_angsd

# 4) Re-annotate any consensus FASTAs as needed (optional)
mitoz annotate --genetic_code 2 --clade Chordata --genome 03_angsd/mitogenome/SampleA.mt.fasta --outprefix 03_angsd/annotation/SampleA.mt --thread_number 8
```

---

## 10) Notes

- **Species names** in `sample_sheet.tsv` must be scientific names (e.g., *Genus species*), and are passed to MitoZ via `--species` exactly as written.
- Avoid aggressive subsampling unless sure about high MDR; small fractions can kill de novo assemblies.
- NUMTs can confuse assemblies—validate with coverage plots and BLAST.

---

**Questions?** Ping Taylor. This README is source-of-truth for the team’s mtDNA pipeline.
