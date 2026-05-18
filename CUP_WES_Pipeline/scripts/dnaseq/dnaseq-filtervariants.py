#!/usr/bin/env python3
"""
CUP WES Pipeline -- Post-VEP Somatic Variant Filter
=====================================================
Filters VEP-annotated Mutect2 VCF by:
  - Minimum VAF (default 5%)
  - Minimum depth (default 20x)
  - Removes synonymous and non-coding consequences
  - Filters to COSMIC Cancer Gene Census genes
  - Excludes variants with specified Mutect2 filter flags
  - Outputs filtered VCF and TSV table

Usage:
    python3 filter_variants.py \
        --vcf <vep_annotated.vcf.gz> \
        --cosmic <Census_all.csv> \
        --outdir <output_dir> \
        [--min_vaf 0.05] \
        [--min_depth 20] \
        [--tier 1] \
        [--exclude_filters germline panel_of_normals strand_bias base_qual weak_evidence]

Requirements:
    pip install cyvcf2 pandas --break-system-packages
"""

import argparse
import os
import sys
import pandas as pd
import cyvcf2

# =============================================================================
# Consequences to keep (protein-altering, high and moderate impact)
# =============================================================================

KEEP_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "splice_region_variant",
}

# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Post-VEP somatic variant filter")
    parser.add_argument("--vcf",            required=True,  help="VEP annotated VCF (gz or uncompressed)")
    parser.add_argument("--cosmic",         required=True,  help="COSMIC Cancer Gene Census CSV (Census_all.csv)")
    parser.add_argument("--outdir",         required=True,  help="Output directory")
    parser.add_argument("--min_vaf",        type=float, default=0.05, help="Minimum VAF (default: 0.05)")
    parser.add_argument("--min_depth",      type=int,   default=20,   help="Minimum depth (default: 20)")
    parser.add_argument("--tier",           type=int,   default=None, help="Filter to COSMIC tier (1 or 2). Default: both")
    parser.add_argument("--exclude_filters",nargs="+",  default=[],
                        help="Mutect2 filter flags to exclude. E.g: germline panel_of_normals strand_bias base_qual weak_evidence clustered_events haplotype slippage multiallelic")
    return parser.parse_args()

# =============================================================================
# Load COSMIC Cancer Gene Census
# =============================================================================

def load_cosmic_genes(cosmic_csv, tier=None):
    df = pd.read_csv(cosmic_csv)
    df.columns = df.columns.str.strip()
    if tier:
        df = df[df["Tier"] == tier]
        print(f"  Loaded {len(df)} COSMIC Tier {tier} genes")
    else:
        print(f"  Loaded {len(df)} COSMIC genes (Tier 1 + Tier 2)")
    return set(df["Gene Symbol"].dropna().str.strip())

# =============================================================================
# Parse VEP CSQ field
# =============================================================================

def parse_csq(variant, csq_fields):
    try:
        csq_raw = variant.INFO.get("CSQ")
        if not csq_raw:
            return []
        annotations = []
        for entry in csq_raw.split(","):
            values = entry.split("|")
            ann = dict(zip(csq_fields, values))
            annotations.append(ann)
        return annotations
    except Exception:
        return []

# =============================================================================
# Main filtering logic
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("\n=== CUP WES -- Somatic Variant Filter ===\n")
    print(f"  Input VCF:        {args.vcf}")
    print(f"  COSMIC list:      {args.cosmic}")
    print(f"  Min VAF:          {args.min_vaf}")
    print(f"  Min depth:        {args.min_depth}")
    print(f"  COSMIC tier:      {args.tier if args.tier else 'All'}")
    print(f"  Exclude filters:  {args.exclude_filters if args.exclude_filters else 'None'}\n")

    # Load COSMIC genes
    print("Loading COSMIC Cancer Gene Census...")
    cancer_genes = load_cosmic_genes(args.cosmic, tier=args.tier)

    # Open VCF
    vcf = cyvcf2.VCF(args.vcf)

    # Parse VEP CSQ field names from header
    csq_fields = []
    for header_line in vcf.header_iter():
        info = header_line.info(extra=True)
        if info.get("ID") == "CSQ":
            desc = info.get("Description", "")
            if "Format:" in desc:
                csq_fields = desc.split("Format:")[1].strip().strip('"').split("|")
            break

    if not csq_fields:
        print("ERROR: Could not parse VEP CSQ field from VCF header.")
        sys.exit(1)

    print(f"  Found {len(csq_fields)} VEP annotation fields\n")

    # Output paths
    sample_name  = os.path.basename(args.vcf).replace(".vcf.gz", "").replace(".vcf", "")
    tier_suffix  = f".tier{args.tier}" if args.tier else ""
    out_vcf_path = os.path.join(args.outdir, f"{sample_name}{tier_suffix}.filtered.vcf")
    out_tsv_path = os.path.join(args.outdir, f"{sample_name}{tier_suffix}.filtered.tsv")

    # Write filtered VCF
    out_vcf = cyvcf2.Writer(out_vcf_path, vcf)

    # Counters
    total              = 0
    pass_filter_flags  = 0
    pass_vaf           = 0
    pass_depth         = 0
    pass_csq           = 0
    pass_gene          = 0
    kept               = 0

    rows = []

    print("Filtering variants...\n")

    for variant in vcf:
        total += 1

        # --- Mutect2 filter flag exclusion ---
        if args.exclude_filters:
            variant_filters = variant.FILTER if variant.FILTER else ""
            if isinstance(variant_filters, str):
                variant_filter_set = set(variant_filters.split(";"))
            else:
                variant_filter_set = set(variant_filters)
            if variant_filter_set.intersection(set(args.exclude_filters)):
                continue
        pass_filter_flags += 1

        # --- VAF filter ---
        try:
            af = variant.format("AF")[0][0]
        except Exception:
            continue
        if af < args.min_vaf:
            continue
        pass_vaf += 1

        # --- Depth filter ---
        try:
            dp = variant.format("DP")[0][0]
        except Exception:
            continue
        if dp < args.min_depth:
            continue
        pass_depth += 1

        # --- Parse VEP annotations ---
        annotations = parse_csq(variant, csq_fields)
        if not annotations:
            continue

        # --- Consequence and gene filter ---
        kept_annotations = []
        for ann in annotations:
            consequence  = ann.get("Consequence", "")
            gene         = ann.get("SYMBOL", "")
            consequences = set(consequence.split("&"))

            if not consequences.intersection(KEEP_CONSEQUENCES):
                continue
            pass_csq += 1

            if gene not in cancer_genes:
                continue
            pass_gene += 1

            kept_annotations.append(ann)

        if not kept_annotations:
            continue

        kept += 1
        out_vcf.write_record(variant)

        # Build TSV rows
        for ann in kept_annotations:
            rows.append({
                "CHROM":            variant.CHROM,
                "POS":              variant.POS,
                "REF":              variant.REF,
                "ALT":              ",".join(variant.ALT),
                "VAF":              round(af, 4),
                "DEPTH":            dp,
                "FILTER":           variant.FILTER,
                "GENE":             ann.get("SYMBOL", ""),
                "CONSEQUENCE":      ann.get("Consequence", ""),
                "HGVSc":            ann.get("HGVSc", ""),
                "HGVSp":            ann.get("HGVSp", ""),
                "IMPACT":           ann.get("IMPACT", ""),
                "EXISTING_VARIANT": ann.get("Existing_variation", ""),
                "AF_gnomAD":        ann.get("gnomADg_AF", ""),
                "CLIN_SIG":         ann.get("CLIN_SIG", ""),
                "SIFT":             ann.get("SIFT", ""),
                "PolyPhen":         ann.get("PolyPhen", ""),
            })

    out_vcf.close()
    vcf.close()

    # Write TSV
    df = pd.DataFrame(rows)
    df.to_csv(out_tsv_path, sep="\t", index=False)

    # ==========================================================================
    # Summary
    # ==========================================================================

    print("=== Filtering Summary ===\n")
    print(f"  Total variants input:            {total:>8,}")
    print(f"  Pass filter flag exclusion:      {pass_filter_flags:>8,}")
    print(f"  Pass VAF >= {args.min_vaf}:              {pass_vaf:>8,}")
    print(f"  Pass depth >= {args.min_depth}:             {pass_depth:>8,}")
    print(f"  Pass consequence filter:         {pass_csq:>8,}")
    print(f"  Pass COSMIC gene filter:         {pass_gene:>8,}")
    print(f"  Final variants kept:             {kept:>8,}")
    print(f"  TSV rows (one per annotation):   {len(rows):>8,}")
    print(f"\n  Output VCF: {out_vcf_path}")
    print(f"  Output TSV: {out_tsv_path}\n")

if __name__ == "__main__":
    main()