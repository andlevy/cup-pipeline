# CUP Pipeline -- Environment Setup Log
## Master Installation Record and Reproducibility Guide

**Machine:** cssd028012  
**OS:** Ubuntu 20.04.x LTS (Focal Fossa)  
**User:** al871  
**Large drive:** `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/` (7.3TB, ~6TB free)  
**Date completed:** May 2026  

---

## Notes Before Starting

The root partition (`/dev/nvme0n1p2`, 1.9TB) is full due to `/home/jeb49` consuming 1.2TB. All large installations must go on the big drive. Temp directories must also be redirected.

The network home directory is at `/its/home/al871` -- this is separate from the root partition and is fine for small config files and scripts.

**Action required:** Contact sysadmin to clear `/home/jeb49`.

---

## 1. Miniconda

Installed to the big drive to avoid root partition space issues.

```bash
# Download installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install to big drive
bash Miniconda3-latest-Linux-x86_64.sh -b -p \
    /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3

# Initialise
/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/bin/conda init bash
source ~/.bashrc
```

**Install location:** `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/`

---

## 2. nfcore Conda Environment

Python 3.11 required -- Ubuntu 20.04 system Python (3.8/3.9) is insufficient for nf-core dependencies.

```bash
conda create -n nfcore python=3.11 -y
conda activate nfcore
pip install nf-core
```

**nf-core version installed:** 4.0.2  
**Environment location:** `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/envs/nfcore/`

---

## 3. Java 17

Java 11 is the system default -- Nextflow 25.x requires Java 17 minimum. Installed into the nfcore conda environment.

```bash
conda activate nfcore
conda install -c conda-forge openjdk=17 -y
```

Verify:
```bash
java -version
# Should report: openjdk version "17.x.x"
```

---

## 4. Apptainer (Singularity)

Installed system-wide via .deb package. Apptainer is the open-source fork of Singularity and defines `singularity` as an alias, so Nextflow's `-profile singularity` works correctly.

```bash
# Download .deb from https://github.com/apptainer/apptainer/releases
# Version installed: 1.3.6
sudo dpkg -i apptainer_1.3.6_amd64.deb
```

**Version:** 1.3.6  
Note: Minor AppArmor warning on Ubuntu 20.04 -- non-fatal and does not affect functionality.

Verify:
```bash
apptainer --version
# apptainer version 1.3.6
```

---

## 5. Nextflow

Already installed on the machine prior to this setup.

**Version:** 25.10.2  
**Location:** `/its/home/al871/.local/bin/nextflow`

Verify:
```bash
nextflow -version
```

---

## 6. nf-core/sarek Pipeline Download

Downloaded via nf-core tools with Singularity containers.

```bash
conda activate nfcore

export NXF_SINGULARITY_CACHEDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/Singularity-cache/nf-sarek/sarek-3.8.1

nf-core pipelines download sarek \
    -r 3.8.1 \
    --container-system singularity \
    --compress none \
    --outdir /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/CUP_WES_Pipeline/sarek-3.8.1
```

**Pipeline version:** 3.8.1  
**Pipeline location:** `CUP_Pipeline/CUP_WES_Pipeline/sarek-3.8.1/3_8_1/main.nf`  
**Singularity images:** `CUP_Pipeline/Singularity-cache/nf-sarek/sarek-3.8.1/`

**Note:** 50/53 images downloaded successfully. Three failures (non-critical):
- `gatk:4.5.0.0` -- replaced by `gatk4-spark:4.6.1.0` which is present
- `clara-parabricks:4.6.0-1` -- requires NVIDIA NGC account, not needed
- `deepvariant:1.9.0` -- Google container access restricted, not needed

---

## 7. Reference Genomes -- iGenomes GATK.GRCh38 Bundle

Downloaded from AWS S3.

```bash
aws s3 sync \
    s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ \
    /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/igenomes/Homo_sapiens/GATK/GRCh38/ \
    --no-sign-request
```

**Location:** `CUP_Pipeline/References/igenomes/`

Key files confirmed present:
- FASTA + .fai + .dict: `Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta`
- BWAmem2 indices: `Sequence/BWAmem2Index/`
- gnomAD germline resource: `Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz`
- dbSNP 146: `Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz`
- Mills indels: `Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
- 1000G PoN: `Annotation/GATKBundle/1000g_pon.hg38.vcf.gz`
- MSIsensor2 models, ASCAT files, chromosome lengths, intervals

---

## 8. VEP Cache

Ensembl VEP cache version 110 for GRCh38. Sarek 3.8.1 expects version 110.

```bash
cd /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/vep_cache

wget https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

tar xzf homo_sapiens_vep_110_GRCh38.tar.gz

rm homo_sapiens_vep_110_GRCh38.tar.gz
```

**Location:** `CUP_Pipeline/References/vep_cache/homo_sapiens/110_GRCh38/`  
**Size:** ~20-30GB

---

## 9. COSMIC Cancer Gene Census

Downloaded from https://cancer.sanger.ac.uk/census (requires free COSMIC account).

**File:** `Census_all.csv`  
**Location:** `CUP_Pipeline/References/cosmic/Census_all.csv`  
**Genes:** 763 total (Tier 1 + Tier 2)

---

## 10. PureCN Conda Environment

R-based tool for tumour purity and ploidy estimation from tumour-only WES data.

```bash
conda create -n purecn -c conda-forge -c bioconda r-base=4.3 bioconductor-purecn -y
```

Additional R package required:
```bash
conda activate purecn
Rscript -e "install.packages('optparse', repos='https://cloud.r-project.org')"
```

**R version:** 4.3.3  
**Environment location:** `/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/envs/purecn/`

Verify:
```bash
conda activate purecn
R -e "library(PureCN); cat('PureCN loaded successfully\n')"
```

---

## 11. Python Packages (base environment)

Required for the post-VEP variant filtering script.

```bash
pip install cyvcf2 pandas --break-system-packages
```

---

## 12. PureCN Interval File (one-time setup)

Generates a GC-corrected interval file from the KAPA HyperExome V2 capture targets BED file. Only needs to be run once per capture kit.

```bash
conda activate purecn

INTERVAL_SCRIPT=$(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")

mkdir -p /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/purecn

Rscript ${INTERVAL_SCRIPT} \
    --infile  /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/KAPA-HyperExomeV2/HyperExomeV2_hg38/HyperExomeV2_capture_targets.bed \
    --fasta   /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
    --outfile /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/References/purecn/HyperExomeV2_hg38_intervals.txt \
    --genome  hg38
```

**Output:** `CUP_Pipeline/References/purecn/HyperExomeV2_hg38_intervals.txt`

---

## 13. Environment Variables (add to run scripts)

These must be set before each pipeline run. They are baked into the run scripts rather than `.bashrc` to keep them self-contained.

```bash
# Redirect temp directories to big drive (avoids filling root partition)
mkdir -p /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_TEMP=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export JAVA_OPTS="-Djava.io.tmpdir=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp"

# Singularity image cache
export NXF_SINGULARITY_CACHEDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline/Singularity-cache/nf-sarek/sarek-3.8.1

# Prevent Nextflow checking for updates online
export NXF_OFFLINE=true
```

---

## 14. Pre-existing References

The following were already present on the machine before this setup:

- KAPA HyperExome V2 BED files: `References/KAPA-HyperExomeV2/HyperExomeV2_hg38/`
- Ensembl 115 FASTA + GTF: `References/ensembl_115/`
- Additional gnomAD + ExAC VCFs: `References/gatk-bundle/`

---

## Directory Structure

```
/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/
├── miniconda3/                          # Conda installation
│   └── envs/
│       ├── nfcore/                      # nf-core + Java 17
│       └── purecn/                      # R 4.3 + PureCN
├── tmp/                                 # Nextflow and Java temp files
└── CUP_Pipeline/
    ├── Data/                            # Patient sample data
    │   └── ICR-Test/
    │       └── Pool_3827/               # Test FASTQ files
    ├── Results/                         # Pipeline outputs per pool
    │   └── Pool_3827/
    ├── References/                      # All reference files (shared)
    │   ├── igenomes/                    # GATK GRCh38 iGenomes bundle
    │   ├── vep_cache/                   # Ensembl VEP cache v110
    │   ├── cosmic/                      # COSMIC Cancer Gene Census
    │   ├── purecn/                      # PureCN interval file
    │   ├── KAPA-HyperExomeV2/           # Capture kit BED files
    │   ├── ensembl_115/                 # Ensembl 115 FASTA + GTF
    │   └── gatk-bundle/                 # Additional GATK resources
    ├── Singularity-cache/
    │   └── nf-sarek/
    │       └── sarek-3.8.1/             # Singularity container images
    ├── MasterScripts/                   # Shared utility scripts
    ├── CUP_WES_Pipeline/
    │   ├── sarek-3.8.1/                 # Pipeline code
    │   ├── configs/
    │   │   └── local.config             # Nextflow local executor config
    │   └── scripts/
    │       ├── dnaseq/                  # Run scripts
    │       └── filter_variants.py       # Post-VEP filtering script
    ├── CUP_RNA_Pipeline/
    │   └── scripts/
    └── CUP_Downstream_Pipeline/
        └── scripts/
```

---

## Outstanding Actions

- [ ] Contact sysadmin to clear `/home/jeb49` (1.2TB on root partition)
- [ ] Set up Seqera Tower for remote pipeline monitoring
- [ ] Complete PureCN interval file generation
- [ ] Write PureCN run script
- [ ] Download nf-core/rnaseq pipeline when needed