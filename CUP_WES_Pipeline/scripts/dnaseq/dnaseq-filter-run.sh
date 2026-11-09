#!/bin/bash
# =============================================================================
# CUP WES Pipeline -- Somatic Variant & CNV Filter Runner
# =============================================================================
# Version: v0.1.0 Andrew Levy
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
POOL=Pool_3827
SAMPLE=S1

# Project conda env providing cyvcf2, pandas, matplotlib, numpy
source /mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/miniconda3/etc/profile.d/conda.sh
conda activate cupwes

COSMIC=${BASE}/References/cosmic/Census_all.csv

VARIANT_SCRIPT=${BASE}/CUP_WES_Pipeline/scripts/dnaseq/dnaseq-filtervariants.py
CNV_SCRIPT=${BASE}/CUP_WES_Pipeline/scripts/dnaseq/dnaseq-filtercnv.py

VCF=${BASE}/Results/${POOL}/annotation/mutect2/${SAMPLE}/${SAMPLE}.mutect2.filtered_VEP.ann.vcf.gz
CNS=${BASE}/Results/${POOL}/variant_calling/cnvkit/${SAMPLE}/${SAMPLE}.tumor_only.call.cns

VARIANT_OUTDIR=${BASE}/Results/${POOL}/annotation/filtered
CNV_OUTDIR=${BASE}/Results/${POOL}/cnv_filtered

EXCLUDE_FILTERS="germline panel_of_normals strand_bias base_qual weak_evidence clustered_events haplotype slippage multiallelic"

echo "=============================================="
echo " CUP WES -- Somatic Variant & CNV Filter"
echo " Pool:   ${POOL}"
echo " Sample: ${SAMPLE}"
echo " Started: $(date)"
echo "=============================================="
echo ""

# =============================================================================
# Somatic variant filter
# =============================================================================

echo "--- Somatic Variant Filter ---"
python3 ${VARIANT_SCRIPT} \
    --vcf             ${VCF} \
    --outdir          ${VARIANT_OUTDIR} \
    --min_vaf         0.05 \
    --min_depth       20 \
    --exclude_filters ${EXCLUDE_FILTERS} \
    --no_consequence_filter \
    --report

# =============================================================================
# CNV filter -- Run 1: all COSMIC tiers
# =============================================================================

echo "--- CNV Filter: All COSMIC tiers ---"
python3 ${CNV_SCRIPT} \
    --cns        ${CNS} \
    --outdir     ${CNV_OUTDIR} \
    --min_probes 10 \
    --min_log2   0.3 \
    --min_weight 0.0 \
    --report


echo ""
echo "Completed: $(date)"
