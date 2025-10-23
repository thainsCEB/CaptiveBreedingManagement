
# DeepVariant + GLnexus Batch Runner — Setup & Usage

**Script:** `deepvariant_glnexus_batch.py`  
**Author:** Taylor Hains  
**Last updated:** 2025-10-18

This guide explains how to set up dependencies, prepare inputs, run the batch pipeline, and interpret outputs for the DeepVariant + GLnexus cohort calling workflow.

---

## 1) What this does

- Runs **DeepVariant** per sample to produce:
  - Genotyped VCF (`*.vcf.gz`) and
  - gVCF (`*.g.vcf.gz`)
  - Optional **visual report** (`*.visual_report.html`)
- Optionally **cleans up** DeepVariant intermediate files to save disk space.
- Merges **all gVCFs** with **GLnexus** into a **cohort BCF**, then converts to **bgzipped VCF** and **indexes** with tabix.
- Supports **Singularity** or **Docker** for DeepVariant. GLnexus can run either in a container or **on the host** (`glnexus_cli` in your `PATH`).

---

## 2) Requirements / Installation

### A. Container runtime
Choose **one** for DeepVariant:
- **Singularity/Apptainer**: `singularity --version` or `apptainer --version`
- **Docker**: `docker --version` and user has permission to run `docker`

### B. DeepVariant image
- **Singularity**: Path to a `.sif` image (e.g., `/opt/deepvariant/deepvariant-<ver>.sif`).
- **Docker**: An available tag (e.g., `google/deepvariant:1.6.0`).

### C. GLnexus
You have two options:
1. **Host install** (recommended on HPC with modules): make sure `glnexus_cli` is in `PATH`.  
   - Example: `module load glnexus/1.4.1`
2. **Container**: provide a GLnexus container image (e.g., `quay.io/mlin/glnexus:v1.4.1`).

### D. bcftools + tabix (host)
Used to convert the GLnexus BCF to VCF.GZ and index with tabix.
- `bcftools --version`
- `tabix --version`

### E. Reference indexing
Reference FASTA must have:
- `genome.fa.fai` (via `samtools faidx genome.fa`)
- BWA index is **not** required by DeepVariant.
- If your DeepVariant image requires a sequence dictionary: `samtools dict` or `picard CreateSequenceDictionary` (rarely needed; DeepVariant typically only needs `.fai`).

### F. Input BAMs
Each BAM should be:
- **Coordinate-sorted** and accompanied by an index (`.bam.bai`).
- Aligned to the **same reference** you pass with `-r/--reference`.

---

## 3) Inputs

### A. Samples list (`samples.tsv` or CSV)
Two columns: `sample_id` and absolute path to BAM. Header optional.

**TSV example:**
```
NA12878    /data/bams/NA12878.sorted.bam
SAMPLE_B   /data/bams/SAMPLE_B.sorted.bam
```

**CSV example:**
```
NA12878,/data/bams/NA12878.sorted.bam
SAMPLE_B,/data/bams/SAMPLE_B.sorted.bam
```

### B. Optional region and BEDs
- `-g/--region "chr:start-end"` for a single interval (tiny tests).
- `-R/--regions-bed regions.bed` for a BED of intervals.
- `-p/--par-bed par_regions.bed` to pass PAR regions to DeepVariant (sex chromosomes).

---

## 4) Basic usage

### A. Singularity for DeepVariant, GLnexus on host
```bash
./deepvariant_glnexus_batch.py   -s samples.tsv   -r /data/ref/genome.fa   -o dv_out   -E singularity -I /opt/deepvariant/deepvariant.sif   -H host   -M WGS   -t 32
```

### B. Docker for DeepVariant, GLnexus in container
```bash
./deepvariant_glnexus_batch.py   -s samples.tsv   -r /data/ref/genome.fa   -o dv_out   -E docker -I google/deepvariant:1.6.0   -H container -G quay.io/mlin/glnexus:v1.4.1   -M WGS   -t 48
```

### C. Small-region test with a visual report & intermediates to /scratch
```bash
./deepvariant_glnexus_batch.py   -s samples.tsv   -r /input/ucsc.hg19.chr20.unittest.fasta   -o /output   -E docker -I google/deepvariant:${BIN_VERSION}   -H host   -M WGS   -g "chr20:10,000,000-10,010,000"   -S   -i /scratch/dv_intermediate   -O output   -t 1
```

---

## 5) Key options

- `-E/--engine`: `singularity` or `docker` for **DeepVariant**.
- `-I/--dv-image`: DeepVariant image (`.sif` or Docker tag).
- `-H/--glx-engine`: `host` to use local `glnexus_cli`; `container` to run in a container.
- `-G/--glx-image`: GLnexus container image (required if `-H container`).
- `-M/--model-type`: `WGS` (default), `WES`, `PACBIO`, or `ONT`.
- `-n/--shards`: DeepVariant shards (`--num_shards`). Defaults to CPU count.
- `-t/--threads`: Threads for **GLnexus** and **bcftools/tabix**. Defaults to CPU count.
- `-g/--region`: One interval (e.g., `chr1:1,000,000-2,000,000`).
- `-R/--regions-bed`: BED of intervals.
- `-p/--par-bed`: PAR BED for sex-chr handling.
- `-S/--vcf-stats-report`: Enable DeepVariant’s `--vcf_stats_report=true`.
- `-i/--intermediate-dir`: Root for DeepVariant intermediates (per-sample subdir created).
- `-K/--keep-intermediates`: Do **not** delete intermediates (default is to delete).
- `-O/--output-prefix`: Controls naming (see below).
- `-x/--dv-extra` and `-X/--glnexus-extra`: Pass-through to DeepVariant and GLnexus.
- `-f/--force`: Re-run even if outputs exist.
- `-N/--dry-run`: Print commands without executing.

---

## 6) Outputs & naming

### Per-sample (written to `OUTDIR/<sample>/`)

If you set `-O PREFIX`:
```
PREFIX.<sampleID>.dv.vcf.gz
PREFIX.<sampleID>.dv.vcf.gz.tbi
PREFIX.<sampleID>.dv.g.vcf.gz
PREFIX.<sampleID>.dv.g.vcf.gz.tbi
PREFIX.<sampleID>.dv.visual_report.html            # when -S is set and report is produced
```
If you omit `-O`, defaults are `<sampleID>.deepvariant.*`.

**Intermediates**
- Default location: `OUTDIR/<sample>/intermediate_results/`
- Or `-i /path/to/root` → `/path/to/root/<sample>/`
- **Cleaned up** after DeepVariant unless `-K` is provided.

### Cohort (written to `OUTDIR/`)

If you set `-O PREFIX`:
```
PREFIX.glnexus.bcf
PREFIX.glnexus.vcf.gz
PREFIX.glnexus.vcf.gz.tbi
```
If you omit `-O`, defaults to `--cohort-name` (default: `cohort` → `cohort.glnexus.*`).

**Logs**
- `OUTDIR/logs/deepvariant/` — per-sample DeepVariant logs, tabix/index logs, preserve-report logs.
- `OUTDIR/logs/glnexus/` — GLnexus logs.
- `OUTDIR/logs/bcftools/` — bcftools conversion logs.

**Manifest**
- `OUTDIR/manifest.deepvariant_glnexus.json` describes the run (inputs, images, outputs).

---

## 7) Resume behavior

- The script **skips** re-running a sample’s DeepVariant if its gVCF already exists (unless `-f/--force`).
- GLnexus is skipped if the cohort BCF already exists (unless `-f/--force`).

---

## 8) Tips & troubleshooting

- **Permissions / mounts**: Make sure all input and output paths are absolute and readable/writable. The script mounts each path into the container at the **same absolute path**.
- **GPU**: Add `-u/--gpu` if your DeepVariant image supports GPU and you have NVIDIA drivers set up.
- **Reference mismatch**: All BAMs must be aligned to the same `-r/--reference` you pass to the script.
- **Regions vs whole genome**: For full WGS, omit `-g` and `-R`. For quick smoke tests, use a small `-g` region.
- **GLnexus memory**: Adjust with `-X/--glnexus-extra`, e.g. `--mem-gbytes 64`.
- **bcftools errors**: Ensure `bcftools` and `tabix` are installed on host; the script runs them outside containers.

---

## 9) Minimal end-to-end example

```bash
# 0) Prepare sample list
printf "NA12878	/data/bams/NA12878.sorted.bam
" > samples.tsv

# 1) Run DeepVariant (Singularity) and GLnexus (host), 32 threads, WGS
./deepvariant_glnexus_batch.py   -s samples.tsv   -r /data/ref/hg38.fa   -o dv_out   -E singularity -I /opt/deepvariant/deepvariant.sif   -H host   -M WGS   -t 32   -O cohort1

# 2) Check outputs
ls -1 dv_out/NA12878/
# -> cohort1.NA12878.dv.vcf.gz, cohort1.NA12878.dv.g.vcf.gz, indexes, optional visual_report.html

ls -1 dv_out/
# -> cohort1.glnexus.bcf, cohort1.glnexus.vcf.gz, cohort1.glnexus.vcf.gz.tbi
```

---

## 10) Reproducibility & citation notes

- Record the exact **DeepVariant** and **GLnexus** versions you used (image tags or module versions).
- The script writes a **JSON manifest** with key settings and paths for future reference.
- Consider pushing your `samples.tsv`, command line, and manifest into version control with your analysis notes.

---

Happy calling! If you need a Snakemake/Nextflow wrapper around this, we can scaffold one to parallelize per-sample DeepVariant and handle cohort runs automatically.
