#!/bin/bash
# =============================================================================
# CUP WES Pipeline -- Production Run Script
# nf-core/sarek 3.8.1
# =============================================================================
# FFPE tumour-only WES somatic variant and CNV calling
#
# Sample type:   FFPE tumour-only
# Assay:         WES (KAPA HyperExome V2)
# Sequencer:     NovaSeq
# Tools:         Mutect2, CNVKit, VEP
# Reference:     GRCh38 (GATK iGenomes bundle)
#
# Usage: bash run_cup_wes.sh
# Version: v0.1.0
# =============================================================================

set -e  # Exit immediately on any error

# --- Base paths --------------------------------------------------------------

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
REFS=${BASE}/References

# --- Environment setup -------------------------------------------------------

# Activate conda environment (provides Java 17 and nf-core tools)
source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate nfcore

# Redirect temp directories to big drive to avoid filling root partition
mkdir -p /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export NXF_TEMP=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp
export JAVA_OPTS="-Djava.io.tmpdir=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp"

# Point Nextflow to the local Singularity image cache
export NXF_SINGULARITY_CACHEDIR=${BASE}/Singularity-cache/nf-sarek/sarek-3.8.1

# Prevent Nextflow from checking for pipeline updates online
export NXF_OFFLINE=true

# --- Sample configuration ----------------------------------------------------

# Pool to process -- update this for each run
POOL=Pool_3827

SAMPLE_DIR=${BASE}/Data/ICR-Test/${POOL}
SAMPLESHEET=${SAMPLE_DIR}/samplesheet.csv
OUTDIR=${BASE}/Results/${POOL}

# --- Pipeline and references -------------------------------------------------

PIPELINE=${BASE}/CUP_WES_Pipeline/sarek-3.8.1/3_8_1/main.nf

IGENOMES_BASE=${REFS}/igenomes
INTERVALS=${REFS}/KAPA-HyperExomeV2/HyperExomeV2_hg38/HyperExomeV2_capture_targets.bed
PON=${REFS}/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
GERMLINE_RESOURCE=${REFS}/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
VEP_CACHE=${REFS}/vep_cache

# --- Run ---------------------------------------------------------------------

echo "=============================================="
echo " CUP WES Pipeline"
echo " Pool:    ${POOL}"
echo " Outdir:  ${OUTDIR}"
echo " Started: $(date)"
echo "=============================================="
echo ""

nextflow run ${PIPELINE} \
    -profile singularity \
    -c ${BASE}/CUP_WES_Pipeline/configs/local.config \
    -resume \
    --input ${SAMPLESHEET} \
    --outdir ${OUTDIR} \
    --genome GATK.GRCh38 \
    --igenomes_base ${IGENOMES_BASE} \
    --intervals ${INTERVALS} \
    --pon ${PON} \
    --germline_resource ${GERMLINE_RESOURCE} \
    --vep_cache ${VEP_CACHE} \
    --vep_genome GRCh38 \
    --vep_species homo_sapiens \
    --vep_cache_version 110 \
    --aligner bwa-mem2 \
    --wes \
    --trim_fastq \
    --trim_nextseq \
    --tools mutect2,cnvkit,vep \
    --skip_tools baserecalibrator \
    --dont_use_soft_clipped_bases \
    -offline

echo ""
echo "=============================================="
echo " Pipeline complete: $(date)"
echo "=============================================="