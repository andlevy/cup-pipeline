# CUP RNA Pipeline -- Environment Setup Log
## Installation Record and Reproducibility Guide

**Machine:** cssd028012  
**OS:** Ubuntu 20.04.x LTS (Focal Fossa)  
**User:** al871  
**Large drive:** `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/` (7.3TB, ~6TB free)  
**Date completed:** May 2026  

See also: `InitialSetupLog.md` for base environment setup (Miniconda, nfcore conda env, Java 17, Apptainer, Nextflow).

---

## Prerequisites

The following from `InitialSetupLog.md` are required before proceeding:

- Section 1: Miniconda installed at `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/`
- Section 2: `nfcore` conda environment with nf-core 4.0.2
- Section 3: Java 17 (in `nfcore` env)
- Section 4: Apptainer (Singularity) 1.3.6 system-wide
- Section 5: Nextflow 25.10.2

---

## 1. nf-core/rnaseq Pipeline Download

Downloaded via nf-core tools with Singularity containers. A convenience script is provided:

```bash
bash /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/CUP_RNA_Pipeline/scripts/rnaseq/rnaseq-local-download.sh
```

Or run manually:

```bash
conda activate nfcore

export NXF_SINGULARITY_CACHEDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/Singularity-cache/nf-rnaseq/rnaseq-3.23.0

mkdir -p ${NXF_SINGULARITY_CACHEDIR}

nf-core pipelines download rnaseq \
    -r 3.23.0 \
    --container-system singularity \
    --compress none \
    --outdir /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/CUP_RNA_Pipeline/rnaseq-3.23.0
```

**Pipeline version:** 3.23.0  
**Pipeline location:** `CUP_Pipeline/CUP_RNA_Pipeline/rnaseq-3.23.0/3_23_0/main.nf`  
**Singularity images:** `CUP_Pipeline/Singularity-cache/nf-rnaseq/rnaseq-3.23.0/`  
**Approximate container size:** ~5–10GB

---

## 2. Reference Files

No download needed. The Ensembl 115 FASTA and GTF were already present on the machine (see `InitialSetupLog.md` section 14).

**FASTA:** `CUP_Pipeline/References/ensembl_115/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`  
**GTF:** `CUP_Pipeline/References/ensembl_115/Homo_sapiens.GRCh38.115.gtf.gz`

**Note on iGenomes:** The `References/igenomes/` directory contains the GATK bundle only (BWAmem2 indices, gnomAD VCFs, etc.) — it does not contain the UCSC/Ensembl iGenomes STAR index or Salmon index used by nf-core/rnaseq. For this reason, the run script uses explicit `--fasta`/`--gtf` paths rather than `--genome GRCh38`, which would attempt to fetch pre-built indices from AWS S3 at runtime (incompatible with offline mode).

**STAR genome index:** Built automatically on the first pipeline run via the `STAR_GENOMEGENERATE` process (~45 minutes, ~30GB output). Cached in the Nextflow work directory and reused on all subsequent `-resume` runs.

---

## 3. Local Configuration

**File:** `CUP_Pipeline/CUP_RNA_Pipeline/configs/local.config`

Key settings:
- Executor: `local`, 52 CPUs, 100GB (capped from 64 CPUs / 125GB to leave headroom)
- Singularity: enabled with `--bind /mnt` and `autoMounts = true`
- `workDir` on big drive to avoid filling root partition
- Process overrides: `STAR_ALIGN` at 16 CPUs / 40GB (allows 2 concurrent alignments within the 100GB cap)

---

## 4. Samplesheet Format

RNA-seq samplesheets use a different format than WES/Sarek. Required columns:

```csv
sample,fastq_1,fastq_2,strandedness
SAMPLE_01,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,auto
```

`strandedness: auto` instructs the pipeline to detect strand orientation automatically using a subset of reads. Place the samplesheet at:

```
CUP_Pipeline/Data/ICR-Test/Pool_XXXX/samplesheet.csv
```

---

## 5. Environment Variables

Set before each pipeline run. Baked into the run scripts rather than `.bashrc` to keep them self-contained.

```bash
export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_TEMP=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export JAVA_OPTS="-Djava.io.tmpdir=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp"
export NXF_SINGULARITY_CACHEDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/Singularity-cache/nf-rnaseq/rnaseq-3.23.0
export NXF_OFFLINE=true
```

---

## 6. Run Scripts

| Script | Purpose |
|--------|---------|
| `scripts/rnaseq/rnaseq-local-download.sh` | One-time download of pipeline + containers (requires internet) |
| `scripts/rnaseq/rnaseq-local-test.sh` | Test run with synthetic data to verify containers and wiring (requires internet) |
| `scripts/rnaseq/rnaseq-local-run-0.1.0.sh` | Production run (fully offline after download) |

---

## Directory Structure

```
CUP_Pipeline/
├── CUP_RNA_Pipeline/
│   ├── configs/
│   │   └── local.config                    # Nextflow local executor + Singularity config
│   ├── rnaseq-3.23.0/                       # Pipeline code (after download)
│   │   └── 3_23_0/
│   │       └── main.nf
│   └── scripts/
│       └── rnaseq/
│           ├── rnaseq-local-download.sh     # One-time download script
│           ├── rnaseq-local-run-0.1.0.sh    # Production run script
│           ├── rnaseq-local-test.sh         # Test run script
│           └── legacy/
│               ├── rnaseq-run.sh            # Original Docker-based script (kept for reference)
│               └── get_genome.sh            # Ensembl reference download (not needed, refs present)
├── References/
│   └── ensembl_115/                         # GRCh38 FASTA + GTF (Ensembl v115)
├── Singularity-cache/
│   └── nf-rnaseq/
│       └── rnaseq-3.23.0/                   # Singularity container images
├── Data/
│   └── ICR-Test/
│       └── Pool_XXXX/
│           └── samplesheet.csv              # RNA-seq samplesheet (per run)
└── Results/
    └── Pool_XXXX/                           # Pipeline outputs (per run)
```

---

## 7. System File Descriptor Limit

The default Linux `fs.file-max` of 65,535 is too low for Nextflow + Singularity on this machine. With 225 container images in the Singularity cache plus concurrent processes, the system exhausts all available file handles and the pipeline fails with `Too many open files in system`.

Increased permanently:

```bash
sudo sysctl -w fs.file-max=2000000
echo "fs.file-max = 2000000" | sudo tee -a /etc/sysctl.conf
```

**Verify current usage:**
```bash
cat /proc/sys/fs/file-nr
# [in_use]  [unused]  [max]
```

**Note:** This must be set on any new machine before running the pipeline.

---

## 8. Run History

| Date | Run | Outcome |
|------|-----|---------|
| 27 May 2026 | `rnaseq-local-test.sh` (test profile) | 219/219 tasks completed |
| 29 May 2026 | `rnaseq-local-run-0.1.0.sh`, Pool_3493 | 49 tasks completed; `RSEQC_READDUPLICATION` killed by out-of-memory (exit 137) |
| 18 Jun 2026 | Resumed run, Pool_3493 | Completed after adding `RSEQC_READDUPLICATION` memory override (20GB × attempt) to `local.config` |

**Note:** `nf-core pipelines download` also places a copy of the container images in `rnaseq-3.23.0/singularity-images/` (~12GB). The run scripts use `Singularity-cache/nf-rnaseq/rnaseq-3.23.0/`; the in-pipeline copy is excluded from version control.

---

## Outstanding Actions

- [x] Run `rnaseq-local-download.sh` to download pipeline and containers
- [x] Run `rnaseq-local-test.sh` to verify environment before first real run
- [x] Create samplesheet for first RNA pool (Pool_3493) and complete first production run
