#!/usr/bin/env python3
"""
CUP WES Pipeline -- CNV Segment Filter
=======================================
Filters CNVKit tumour-only call segments by:
  - Minimum number of probes per segment
  - Minimum log2 ratio deviation from diploid
  - Minimum segment weight
  - Filters to COSMIC Cancer Gene Census genes
  - Outputs filtered TSV with one row per gene per segment

Usage:
    python3 filter_cnv.py \
        --cns <S1.tumor_only.call.cns> \
        --cosmic <Census_all.csv> \
        --outdir <output_dir> \
        [--min_probes 10] \
        [--min_log2 0.3] \
        [--min_weight 0.0] \
        [--tier 1]

Requirements:
    pip install pandas --break-system-packages
"""

import argparse
import os
import sys
import re
import pandas as pd

# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="CNV segment filter")
    parser.add_argument("--cns",        required=True, help="CNVKit tumour-only call file (S1.tumor_only.call.cns)")
    parser.add_argument("--cosmic",     required=True, help="COSMIC Cancer Gene Census CSV (Census_all.csv)")
    parser.add_argument("--outdir",     required=True, help="Output directory")
    parser.add_argument("--min_probes", type=int,   default=10,  help="Minimum probes per segment (default: 10)")
    parser.add_argument("--min_log2",   type=float, default=0.3, help="Minimum absolute log2 ratio deviation (default: 0.3)")
    parser.add_argument("--min_weight", type=float, default=0.0, help="Minimum segment weight (default: 0.0)")
    parser.add_argument("--tier",       type=int,   default=None, help="Filter to COSMIC tier (1 or 2). Default: both")
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
# Parse gene symbols from CNVKit gene column
# =============================================================================

def parse_gene_symbols(gene_field):
    """Extract unique gene symbols from CNVKit's semicolon-delimited gene field."""
    if pd.isna(gene_field) or gene_field == ".":
        return []
    symbols = set()
    for entry in gene_field.split(","):
        match = re.search(r'gene_symbol=([^;,]+)', entry)
        if match:
            symbol = match.group(1).strip()
            if symbol:
                symbols.add(symbol)
    return sorted(symbols)

# =============================================================================
# Main filtering logic
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("\n=== CUP WES -- CNV Segment Filter ===\n")
    print(f"  Input CNS:    {args.cns}")
    print(f"  COSMIC list:  {args.cosmic}")
    print(f"  Min probes:   {args.min_probes}")
    print(f"  Min |log2|:   {args.min_log2}")
    print(f"  Min weight:   {args.min_weight}")
    print(f"  COSMIC tier:  {args.tier if args.tier else 'All'}\n")

    # Load COSMIC genes
    print("Loading COSMIC Cancer Gene Census...")
    cancer_genes = load_cosmic_genes(args.cosmic, tier=args.tier)

    # Load CNS file
    print(f"Loading CNS file...")
    df = pd.read_csv(args.cns, sep="\t")
    print(f"  Loaded {len(df)} segments\n")

    # Counters
    total       = len(df)
    pass_probes = 0
    pass_log2   = 0
    pass_weight = 0
    kept        = 0
    rows        = []

    print("Filtering segments...\n")

    for _, row in df.iterrows():

        # --- Probe filter ---
        if row["probes"] < args.min_probes:
            continue
        pass_probes += 1

        # --- Log2 ratio filter ---
        if abs(row["log2"]) < args.min_log2:
            continue
        pass_log2 += 1

        # --- Weight filter ---
        if row["weight"] < args.min_weight:
            continue
        pass_weight += 1

        # --- Parse gene symbols ---
        gene_symbols = parse_gene_symbols(row["gene"])
        if not gene_symbols:
            continue

        # --- Filter to COSMIC cancer genes ---
        cosmic_genes_in_segment = [g for g in gene_symbols if g in cancer_genes]
        if not cosmic_genes_in_segment:
            continue

        kept += 1

        # Build one row per COSMIC gene in this segment
        call_type = "GAIN" if row["log2"] > 0 else "LOSS"
        for gene in cosmic_genes_in_segment:
            rows.append({
                "CHROM":      row["chromosome"],
                "START":      row["start"],
                "END":        row["end"],
                "GENE":       gene,
                "LOG2":       round(row["log2"], 4),
                "CN":         int(row["cn"]) if not pd.isna(row["cn"]) else "",
                "CALL":       call_type,
                "DEPTH":      round(row["depth"], 1),
                "PROBES":     int(row["probes"]),
                "WEIGHT":     round(row["weight"], 4),
                "CI_HI":      round(row["ci_hi"], 4),
                "CI_LO":      round(row["ci_lo"], 4),
            })

    # Write output
    sample_name = os.path.basename(args.cns).replace(".tumor_only.call.cns", "")
    tier_suffix = f".tier{args.tier}" if args.tier else ""
    out_tsv = os.path.join(args.outdir, f"{sample_name}.cnv{tier_suffix}.filtered.tsv")

    out_df = pd.DataFrame(rows)
    if not out_df.empty:
        out_df = out_df.sort_values(["CHROM", "START"])
    out_df.to_csv(out_tsv, sep="\t", index=False)

    # Summary
    print("=== Filtering Summary ===\n")
    print(f"  Total segments input:            {total:>8,}")
    print(f"  Pass probes >= {args.min_probes}:               {pass_probes:>8,}")
    print(f"  Pass |log2| >= {args.min_log2}:              {pass_log2:>8,}")
    print(f"  Pass weight >= {args.min_weight}:              {pass_weight:>8,}")
    print(f"  Segments with COSMIC gene hits:  {kept:>8,}")
    print(f"  Output rows (one per gene):      {len(rows):>8,}")
    print(f"\n  Output TSV: {out_tsv}\n")

if __name__ == "__main__":
    main()