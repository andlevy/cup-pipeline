# CUP Pipeline

Genomic analysis pipelines for **Cancer of Unknown Primary (CUP)** tumour samples. The project combines:

- **Whole-exome sequencing (WES):** somatic SNV/indel and copy-number calling from FFPE tumour-only samples.
- **RNA-seq:** gene and transcript expression quantification.
- **Downstream analysis:** from mutations and expression to gene dependencies, drugs and clinical trials.

The pipelines run on a single local Linux workstation using Nextflow with Singularity/Apptainer containers. After a one-time download, all production runs are fully **offline**.

---

## Repository Layout

```
CUP_Pipeline/
├── CUP_WES_Pipeline/          # WES: nf-core/sarek 3.8.1 + PureCN + post-pipeline filtering
│   ├── sarek-3.8.1/3_8_1/     #   Pipeline source (nf-core download, unmodified)
│   ├── configs/local.config   #   Local executor resources + Singularity settings
│   └── scripts/dnaseq/        #   Run, PureCN and filter scripts (legacy/ = superseded)
├── CUP_RNA_Pipeline/          # RNA-seq: nf-core/rnaseq 3.23.0
│   ├── rnaseq-3.23.0/3_23_0/  #   Pipeline source (nf-core download, unmodified)
│   ├── configs/local.config
│   └── scripts/rnaseq/        #   Download, test and production run scripts
├── CUP_DependANT_Pipeline/    # Downstream: LoF/GoF calling, essential-gene prediction,
│                              #   gene→drug and clinical-trial matching
├── CUP_GeneTargetDatabase/    # Open Targets Platform dataset download
├── envs/                      # Conda environment definitions
└── Documents/                 # Setup logs / reproducibility guides
```

Not version-controlled (see `.gitignore`): `Data/` (patient FASTQs), `References/` (genomes, VEP cache, COSMIC), `Results/`, `Singularity-cache/`, Nextflow `work/` directories and logs.

---

## WES Pipeline (CUP_WES_Pipeline)

**Input:** FFPE tumour-only WES (KAPA HyperExome V2 capture, NovaSeq), GRCh38.

| Step | Tool | Script |
|------|------|--------|
| Trimming, alignment, duplicate marking | fastp, BWA-MEM2, GATK MarkDuplicates | `dnaseq-local-run-0.1.0.sh` |
| Somatic SNV/indel calling | GATK Mutect2 (tumour-only, 1000G PoN, gnomAD germline resource) | 〃 |
| Copy-number calling | CNVkit | 〃 |
| Annotation | Ensembl VEP v110 | 〃 |
| Purity / ploidy | PureCN | `dnaseq-purecn-setup.sh` (once), `dnaseq-purecn-run.sh` |
| Variant filtering | `dnaseq-filtervariants.py` | `dnaseq-filter-run.sh` or `dnaseq-filtervariants-run.sh` |
| CNV filtering | `dnaseq-filtercnv.py` | `dnaseq-filter-run.sh` or `dnaseq-filtercnv-run.sh` |

**Variant filter** (`dnaseq-filtervariants.py`): filters the VEP-annotated Mutect2 VCF on VAF (default ≥5%), depth (default ≥20×), Mutect2 FILTER flags, VEP consequence (protein-altering by default) and, optionally, COSMIC Cancer Gene Census membership or tier. Outputs a filtered VCF and TSV.

**CNV filter** (`dnaseq-filtercnv.py`): filters CNVkit segments on probe count, |log2| deviation and segment weight, with optional COSMIC filtering. Outputs one row per gene per segment.

Both filters take `--report` to produce a QC PDF/PNG report (a filtering funnel plus diagnostic plots) and a stats TSV, so every threshold choice can be audited.

Key design decisions:
- `--skip_tools baserecalibrator`: BQSR is skipped for FFPE tumour-only data.
- `--dont_use_soft_clipped_bases`: reduces FFPE artefact calls.
- Resource overrides in `configs/local.config` replace Sarek's WGS-sized defaults with WES-appropriate values.

## RNA-seq Pipeline (CUP_RNA_Pipeline)

nf-core/rnaseq 3.23.0 runs with Trim Galore, STAR alignment and Salmon quantification (`--aligner star_salmon`) against local Ensembl GRCh38 v115 FASTA/GTF. Explicit `--fasta`/`--gtf` are used instead of `--genome` so runs stay offline; the STAR index is built on the first run and reused via `-resume`.

## Downstream Analysis (CUP_DependANT_Pipeline)

Python application that takes patient mutation and expression data and:
1. Identifies loss- and gain-of-function mutations (VEP, AlphaMissense).
2. Predicts essential (dependency) genes with a trained classifier.
3. Maps essential genes to drugs (ChEMBL).
4. Finds matching clinical trials (ClinicalTrials.gov API).

See [CUP_DependANT_Pipeline/dep-scripts/README.md](CUP_DependANT_Pipeline/dep-scripts/README.md).

## Gene Target Database (CUP_GeneTargetDatabase)

Downloads a curated set of Open Targets Platform datasets (Parquet) for target prioritisation. See [CUP_GeneTargetDatabase/README.md](CUP_GeneTargetDatabase/README.md).

---

## Getting Started

Full step-by-step setup, with versions and the problems hit along the way:

1. [Documents/InitialSetupLog.md](Documents/InitialSetupLog.md): base environment (Conda, Java 17, Apptainer, Nextflow), Sarek download, references, VEP cache, COSMIC, PureCN.
2. [Documents/RNAseqSetupLog.md](Documents/RNAseqSetupLog.md): RNA-seq pipeline download, config, samplesheet format, run history.

Conda environments:

```bash
conda env create -f envs/nfcore.yml    # nf-core tools + Java 17 (Nextflow runs)
conda env create -f envs/purecn.yml    # PureCN (R 4.3)
conda env create -f envs/cupwes.yml    # WES post-pipeline filtering scripts
conda env create -f envs/OTP_GeneSelection.yml   # Open Targets work
```

Typical WES run:

```bash
# 1. Edit POOL in the script, then:
bash CUP_WES_Pipeline/scripts/dnaseq/dnaseq-local-run-0.1.0.sh
# 2. Purity/ploidy
bash CUP_WES_Pipeline/scripts/dnaseq/dnaseq-purecn-run.sh
# 3. Filter variants and CNVs, with QC reports
bash CUP_WES_Pipeline/scripts/dnaseq/dnaseq-filter-run.sh
```

Typical RNA-seq run:

```bash
bash CUP_RNA_Pipeline/scripts/rnaseq/rnaseq-local-run-0.1.0.sh
```

## Software Versions

| Component | Version |
|-----------|---------|
| Nextflow | 25.10.2 |
| nf-core tools | 4.0.2 |
| nf-core/sarek | 3.8.1 |
| nf-core/rnaseq | 3.23.0 |
| Apptainer | 1.3.6 |
| Java (OpenJDK) | 17 |
| Ensembl VEP cache | 110 (GRCh38) |
| Reference (WES) | GATK GRCh38 iGenomes bundle |
| Reference (RNA) | Ensembl GRCh38 release 115 |

## Validation Status

| Run | Status |
|-----|--------|
| Sarek test profile | Completed |
| WES Pool_3827 (Sarek → PureCN → filters) | Completed |
| RNA-seq test profile | Completed |
| RNA-seq Pool_3493 | Completed |

## Notes

- Paths in run scripts are absolute to the development workstation. Change `BASE` at the top of each script to run elsewhere.
- Pipeline source under `sarek-3.8.1/` and `rnaseq-3.23.0/` is the unmodified nf-core release. All customisation lives in `configs/` and `scripts/`.
