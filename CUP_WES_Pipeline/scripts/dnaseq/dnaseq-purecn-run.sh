#!/bin/bash
# =============================================================================
# CUP WES Pipeline -- PureCN Run Script
# =============================================================================
# Estimates tumour purity and ploidy from tumour-only WES data using:
#   - CNVKit log2 ratio file (.cnr)
#   - Mutect2 unfiltered VCF
#   - PureCN GC-corrected interval file
#
# Usage: bash dnaseq-purecn-run.sh
# Version: v0.1.0
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
POOL=Pool_3827
SAMPLE=S1

# --- Activate PureCN conda environment ---------------------------------------

source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate purecn

# --- Redirect temp to big drive ----------------------------------------------

mkdir -p /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp/Rtmp
export TMPDIR=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/tmp/Rtmp

# --- Paths -------------------------------------------------------------------

CNR=${BASE}/Results/${POOL}/variant_calling/cnvkit/${SAMPLE}/${SAMPLE}.cnr
VCF=${BASE}/Results/${POOL}/variant_calling/mutect2/${SAMPLE}/${SAMPLE}.mutect2.vcf.gz
INTERVALS=${BASE}/References/purecn/HyperExomeV2_hg38_intervals.txt
OUTDIR=${BASE}/Results/${POOL}/purecn/${SAMPLE}

# Get PureCN main script location
PURECN_SCRIPT=$(Rscript -e "cat(system.file('extdata', 'PureCN.R', package='PureCN'))")

# --- Create output directory -------------------------------------------------

mkdir -p ${OUTDIR}

# --- Run ---------------------------------------------------------------------

echo "=============================================="
echo " PureCN -- Purity and Ploidy Estimation"
echo " Pool:    ${POOL}"
echo " Sample:  ${SAMPLE}"
echo " Outdir:  ${OUTDIR}"
echo " Started: $(date)"
echo "=============================================="
echo ""

Rscript ${PURECN_SCRIPT} \
    --out             ${OUTDIR}/${SAMPLE} \
    --sampleid        ${SAMPLE} \
    --tumor           ${CNR} \
    --vcf             ${VCF} \
    --intervals       ${INTERVALS} \
    --genome          hg38 \
    --fun-segmentation Hclust \
    --min-af          0.03 \
    --error           0.001 \
    --post-optimize \
    --seed            123

echo ""
echo "=============================================="
echo " PureCN complete: $(date)"
echo " Results: ${OUTDIR}"
echo "=============================================="
