#!/bin/bash
# =============================================================================
# CUP WES Pipeline -- CNV Filter Runner
# =============================================================================
# Version: v0.1.0
# =============================================================================

set -e

BASE=/mnt/913bddcf-fa31-4a34-a763-fbc307b90cde/CUP_Pipeline
POOL=Pool_3827
SAMPLE=S1

CNS=${BASE}/Results/${POOL}/variant_calling/cnvkit/${SAMPLE}/${SAMPLE}.tumor_only.call.cns
COSMIC=${BASE}/References/cosmic/Census_all.csv
OUTDIR=${BASE}/Results/${POOL}/cnv_filtered
SCRIPT=${BASE}/CUP_WES_Pipeline/scripts/dnaseq/dnaseq-filtercnv.py

echo "=============================================="
echo " CUP WES -- CNV Filter"
echo " Pool:    ${POOL}"
echo " Started: $(date)"
echo "=============================================="
echo ""

# --- Run 1: All COSMIC tiers ---
echo "Run 1: All COSMIC tiers..."
python3 ${SCRIPT} \
    --cns        ${CNS} \
    --cosmic     ${COSMIC} \
    --outdir     ${OUTDIR} \
    --min_probes 10 \
    --min_log2   0.3 \
    --min_weight 0.0 \
    --report

# --- Run 2: Tier 1 only ---
echo "Run 2: Tier 1 only..."
python3 ${SCRIPT} \
    --cns        ${CNS} \
    --cosmic     ${COSMIC} \
    --outdir     ${OUTDIR} \
    --min_probes 10 \
    --min_log2   0.3 \
    --min_weight 50.0 \
    --tier       1 \
    --report

echo ""
echo "Completed: $(date)"
