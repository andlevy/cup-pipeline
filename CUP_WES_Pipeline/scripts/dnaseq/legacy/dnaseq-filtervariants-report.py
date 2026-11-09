#!/usr/bin/env python3
"""
CUP WES Pipeline -- Somatic Variant Filter QC Report
=====================================================
Generates diagnostic figures for each filter step applied by
dnaseq-filtervariants.py. Accepts identical arguments so it can be
run alongside the filter script with the same command-line call.

Outputs (under {outdir}/reports/dependantfilter/):
  {sample}.filter_qc.pdf          -- all figures, one per page
  {sample}.filter_qc_stats.tsv    -- per-step variant counts
  png/{01..07}_*.png              -- individual figures
  svg/{01..07}_*.svg              -- individual figures

Usage:
    python3 dnaseq-filtervariants-report.py \
        --vcf <vep_annotated.vcf.gz> \
        --cosmic <Census_all.csv> \
        --outdir <output_dir> \
        [--min_vaf 0.05] \
        [--min_depth 20] \
        [--tier 1] \
        [--exclude_filters germline panel_of_normals strand_bias base_qual weak_evidence]

Requirements:
    pip install cyvcf2 pandas matplotlib numpy --break-system-packages
"""

import argparse
import collections
import os
import sys

import cyvcf2
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


# =============================================================================
# Consequences kept by dnaseq-filtervariants.py (must stay in sync)
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

PALETTE = {
    "pass":    "#2ecc71",
    "fail":    "#e74c3c",
    "neutral": "#95a5a6",
    "blue":    "#3498db",
    "orange":  "#e67e22",
    "tier1":   "#2ecc71",
    "tier2":   "#3498db",
}


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Somatic variant filter QC report")
    p.add_argument("--vcf",             required=True,  help="VEP-annotated VCF (gz or uncompressed)")
    p.add_argument("--cosmic",          required=True,  help="COSMIC Cancer Gene Census CSV (Census_all.csv)")
    p.add_argument("--outdir",          required=True,  help="Output directory (figures go in reports/dependantfilter/)")
    p.add_argument("--min_vaf",         type=float, default=0.05, help="Minimum VAF threshold (default: 0.05)")
    p.add_argument("--min_depth",       type=int,   default=20,   help="Minimum depth threshold (default: 20)")
    p.add_argument("--tier",            type=int,   default=None, help="COSMIC tier filter (1 or 2; default: both)")
    p.add_argument("--exclude_filters", nargs="+",  default=[],
                   help="Mutect2 filter flags to exclude (e.g. germline panel_of_normals strand_bias)")
    return p.parse_args()


# =============================================================================
# COSMIC loader -- returns gene -> tier dict (not a set, to enable tier plots)
# =============================================================================

def load_cosmic_genes(cosmic_csv, tier=None):
    df = pd.read_csv(cosmic_csv)
    df.columns = df.columns.str.strip()
    if tier:
        df = df[df["Tier"] == tier]
    gene_tier = (
        df.dropna(subset=["Gene Symbol"])
          .assign(**{"Gene Symbol": lambda d: d["Gene Symbol"].str.strip()})
          .set_index("Gene Symbol")["Tier"]
          .to_dict()
    )
    tier_label = f"Tier {tier}" if tier else "Tier 1 + 2"
    print(f"  Loaded {len(gene_tier)} COSMIC genes ({tier_label})")
    return gene_tier


# =============================================================================
# VEP CSQ parser
# =============================================================================

def parse_csq(variant, csq_fields):
    try:
        csq_raw = variant.INFO.get("CSQ")
        if not csq_raw:
            return []
        return [
            dict(zip(csq_fields, entry.split("|")))
            for entry in csq_raw.split(",")
        ]
    except Exception:
        return []


# =============================================================================
# Single-pass metric collection
# =============================================================================

def collect_metrics(args, csq_fields, gene_tier):
    """
    One pass through the VCF. Applies filter steps in order, recording
    per-step counts and collecting raw values for plotting.
    """
    exclude_set = set(args.exclude_filters)

    # Raw values for continuous distributions (collected after flag filter)
    all_vafs   = []
    all_depths = []

    # (vaf, depth, pass_vaf, pass_depth) tuples for the scatter plot
    scatter = []

    # Flag counts across ALL variants
    flag_counts = collections.Counter()

    # Consequence counts among variants passing flag + VAF + depth
    csq_counts = collections.Counter()

    # Gene counts (all genes seen) among variants passing flag + VAF + depth + consequence
    gene_counts = collections.Counter()

    # Funnel step counters
    n_total      = 0
    n_pass_flags = 0
    n_pass_vaf   = 0
    n_pass_depth = 0
    n_pass_csq   = 0
    n_pass_gene  = 0

    vcf = cyvcf2.VCF(args.vcf)
    for variant in vcf:
        n_total += 1

        # --- Collect Mutect2 filter flags for every variant ---
        raw_filter = variant.FILTER
        if not raw_filter or raw_filter == "PASS":
            flag_counts["PASS"] += 1
            flags = set()
        else:
            if isinstance(raw_filter, str):
                flags = set(raw_filter.split(";")) - {"", "PASS"}
            else:
                flags = set(raw_filter) - {None, "", "PASS"}
            for f in flags:
                flag_counts[f] += 1

        # --- Filter: Mutect2 flag exclusion ---
        if exclude_set and flags.intersection(exclude_set):
            continue
        n_pass_flags += 1

        # --- Read VAF and depth (needed for both plots and filters) ---
        try:
            af = variant.format("AF")[0][0]
        except Exception:
            continue
        try:
            dp = int(variant.format("DP")[0][0])
        except Exception:
            continue

        all_vafs.append(af)
        all_depths.append(dp)

        pv = af >= args.min_vaf
        pd_ = dp >= args.min_depth
        scatter.append((af, dp, pv, pd_))

        # --- Filter: VAF ---
        if not pv:
            continue
        n_pass_vaf += 1

        # --- Filter: depth ---
        if not pd_:
            continue
        n_pass_depth += 1

        # --- Parse VEP annotations ---
        annotations = parse_csq(variant, csq_fields)
        if not annotations:
            continue

        # Consequence counts (before consequence filter)
        for ann in annotations:
            for term in ann.get("Consequence", "").split("&"):
                t = term.strip()
                if t:
                    csq_counts[t] += 1

        # --- Filter: consequence ---
        kept_anns = [
            ann for ann in annotations
            if set(ann.get("Consequence", "").split("&")).intersection(KEEP_CONSEQUENCES)
        ]
        if not kept_anns:
            continue
        n_pass_csq += 1

        # Gene counts (before COSMIC filter)
        genes_seen = set()
        for ann in kept_anns:
            g = ann.get("SYMBOL", "").strip()
            if g and g not in genes_seen:
                gene_counts[g] += 1
                genes_seen.add(g)

        # --- Filter: COSMIC gene ---
        if not any(ann.get("SYMBOL", "").strip() in gene_tier for ann in kept_anns):
            continue
        n_pass_gene += 1

    vcf.close()

    funnel = {
        "Total input":                       n_total,
        "Pass flag exclusion":               n_pass_flags,
        f"Pass VAF ≥ {args.min_vaf}":   n_pass_vaf,
        f"Pass depth ≥ {args.min_depth}x": n_pass_depth,
        "Pass consequence filter":           n_pass_csq,
        "Pass COSMIC gene filter":           n_pass_gene,
    }

    return {
        "funnel":      funnel,
        "all_vafs":    np.array(all_vafs, dtype=float),
        "all_depths":  np.array(all_depths, dtype=float),
        "scatter":     scatter,
        "flag_counts": flag_counts,
        "csq_counts":  csq_counts,
        "gene_counts": gene_counts,
        "gene_tier":   gene_tier,
    }


# =============================================================================
# Figure 1 -- Filtering funnel
# =============================================================================

def fig_funnel(metrics, args):
    funnel = metrics["funnel"]
    labels = list(funnel.keys())
    counts = list(funnel.values())

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [PALETTE["blue"]] + [PALETTE["pass"]] * (len(labels) - 2) + [PALETTE["tier1"]]

    bars = ax.barh(range(len(labels)), counts, color=colors, edgecolor="white", height=0.6)

    for i, (bar, count) in enumerate(zip(bars, counts)):
        prev = counts[i - 1] if i > 0 else count
        pct  = (1 - count / prev) * 100 if (i > 0 and prev > 0) else 0
        label = f"  {count:,}"
        if i > 0:
            label += f"   (removed {pct:.1f}%)"
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                label, va="center", ha="left", fontsize=9)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Variant count")
    ax.set_title("Filtering Funnel — cumulative effect of each step", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, max(counts) * 1.25)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 2 -- Mutect2 filter flag distribution
# =============================================================================

def fig_filter_flags(metrics, args):
    flag_counts = metrics["flag_counts"]
    exclude_set = set(args.exclude_filters)

    sorted_flags = sorted(flag_counts.items(), key=lambda x: x[1], reverse=True)
    labels  = [f for f, _ in sorted_flags]
    counts  = [c for _, c in sorted_flags]
    colors  = [PALETTE["fail"] if f in exclude_set else PALETTE["neutral"] for f in labels]

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.45 + 1)))
    bars = ax.barh(range(len(labels)), counts, color=colors, edgecolor="white", height=0.6)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  {count:,}", va="center", ha="left", fontsize=9)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Variant count")
    ax.set_title("Mutect2 Filter Flag Distribution", fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.2 if counts else 1)

    legend_handles = [
        mpatches.Patch(color=PALETTE["fail"],    label="Excluded (in --exclude_filters)"),
        mpatches.Patch(color=PALETTE["neutral"],  label="Not excluded"),
    ]
    if exclude_set:
        ax.legend(handles=legend_handles, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 3 -- VAF distribution
# =============================================================================

def fig_vaf(metrics, args):
    vafs      = metrics["all_vafs"]
    threshold = args.min_vaf

    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 1, 51)

    ax.hist(vafs[vafs < threshold],  bins=bins, color=PALETTE["fail"], alpha=0.75,
            label=f"Excluded  (VAF < {threshold})")
    ax.hist(vafs[vafs >= threshold], bins=bins, color=PALETTE["pass"], alpha=0.75,
            label=f"Kept  (VAF ≥ {threshold})")

    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}")

    n_excl = int((vafs < threshold).sum())
    n_kept = int((vafs >= threshold).sum())
    ax.text(0.98, 0.97,
            f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Variant Allele Frequency (VAF)")
    ax.set_ylabel("Variant count")
    ax.set_title("VAF Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 4 -- Depth distribution
# =============================================================================

def fig_depth(metrics, args):
    depths    = metrics["all_depths"]
    threshold = args.min_depth

    fig, ax = plt.subplots(figsize=(9, 5))

    cap = float(np.percentile(depths, 99)) if len(depths) > 0 else 500.0
    cap = max(cap, threshold * 2)
    bins = np.linspace(0, cap, 60)

    ax.hist(depths[depths < threshold],  bins=bins, color=PALETTE["fail"], alpha=0.75,
            label=f"Excluded  (depth < {threshold}x)")
    ax.hist(depths[depths >= threshold], bins=bins, color=PALETTE["pass"], alpha=0.75,
            label=f"Kept  (depth ≥ {threshold}x)")

    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}x")

    n_excl = int((depths < threshold).sum())
    n_kept = int((depths >= threshold).sum())
    ax.text(0.98, 0.97,
            f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Read Depth")
    ax.set_ylabel("Variant count")
    ax.set_title("Read Depth Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 5 -- VAF vs depth scatter (joint threshold view)
# =============================================================================

def fig_vaf_depth_scatter(metrics, args):
    scatter = metrics["scatter"]
    if not scatter:
        return None

    vafs    = np.array([x[0] for x in scatter])
    depths  = np.array([x[1] for x in scatter], dtype=float)
    pv_arr  = np.array([x[2] for x in scatter])
    pd_arr  = np.array([x[3] for x in scatter])

    color_map = [
        PALETTE["pass"]   if pv and pd else
        PALETTE["orange"] if (not pv and pd) else
        PALETTE["blue"]   if (pv and not pd) else
        PALETTE["fail"]
        for pv, pd in zip(pv_arr, pd_arr)
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(vafs, depths, c=color_map, alpha=0.35, s=6, linewidths=0)

    ax.axvline(args.min_vaf,   color="black", linestyle="--", linewidth=1.2,
               label=f"VAF = {args.min_vaf}")
    ax.axhline(args.min_depth, color="black", linestyle="-.", linewidth=1.2,
               label=f"Depth = {args.min_depth}x")

    q = {
        "both":       int((pv_arr  & pd_arr).sum()),
        "fail_vaf":   int((~pv_arr & pd_arr).sum()),
        "fail_depth": int((pv_arr  & ~pd_arr).sum()),
        "both_fail":  int((~pv_arr & ~pd_arr).sum()),
    }
    kw = dict(transform=ax.transAxes, fontsize=8, fontweight="bold")
    ax.text(0.98, 0.98, f"Both pass\n{q['both']:,}",      ha="right", va="top",    color=PALETTE["pass"],   **kw)
    ax.text(0.02, 0.98, f"Fail VAF\n{q['fail_vaf']:,}",   ha="left",  va="top",    color=PALETTE["orange"], **kw)
    ax.text(0.98, 0.02, f"Fail depth\n{q['fail_depth']:,}",ha="right", va="bottom", color=PALETTE["blue"],   **kw)
    ax.text(0.02, 0.02, f"Both fail\n{q['both_fail']:,}", ha="left",  va="bottom", color=PALETTE["fail"],   **kw)

    ax.set_xlabel("VAF")
    ax.set_ylabel("Read Depth")
    ax.set_title("VAF vs Read Depth — joint threshold view", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 6 -- VEP consequence distribution
# =============================================================================

def fig_consequences(metrics, args):
    csq_counts = metrics["csq_counts"]
    if not csq_counts:
        return None

    sorted_csq = sorted(csq_counts.items(), key=lambda x: x[1], reverse=True)
    labels  = [c for c, _ in sorted_csq]
    counts  = [n for _, n in sorted_csq]
    colors  = [PALETTE["pass"] if c in KEEP_CONSEQUENCES else PALETTE["fail"] for c in labels]

    fig, ax = plt.subplots(figsize=(10, max(5, len(labels) * 0.38 + 1)))
    bars = ax.barh(range(len(labels)), counts, color=colors, edgecolor="white", height=0.65)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  {count:,}", va="center", ha="left", fontsize=8)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Annotation count")
    ax.set_title(
        "VEP Consequence Distribution\n"
        "(variants passing flag / VAF / depth filters)",
        fontweight="bold"
    )
    ax.set_xlim(0, max(counts) * 1.2)

    legend_handles = [
        mpatches.Patch(color=PALETTE["pass"], label="Kept (protein-altering)"),
        mpatches.Patch(color=PALETTE["fail"], label="Excluded"),
    ]
    ax.legend(handles=legend_handles)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 7 -- COSMIC gene coverage
# =============================================================================

def fig_cosmic(metrics, args):
    gene_counts = metrics["gene_counts"]
    gene_tier   = metrics["gene_tier"]

    if not gene_counts:
        return None

    n_tier1 = sum(c for g, c in gene_counts.items() if gene_tier.get(g) == 1)
    n_tier2 = sum(c for g, c in gene_counts.items() if gene_tier.get(g) == 2)
    n_non   = sum(c for g, c in gene_counts.items() if g not in gene_tier)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # --- Panel A: donut chart ---
    sizes  = [n_tier1, n_tier2, n_non]
    clrs   = [PALETTE["tier1"], PALETTE["tier2"], PALETTE["neutral"]]
    donut_labels = [
        f"COSMIC Tier 1\n({n_tier1:,})",
        f"COSMIC Tier 2\n({n_tier2:,})",
        f"Non-COSMIC\n({n_non:,})",
    ]
    non_zero = [(s, l, c) for s, l, c in zip(sizes, donut_labels, clrs) if s > 0]
    if non_zero:
        s_, l_, c_ = zip(*non_zero)
        ax1.pie(s_, labels=l_, colors=c_, autopct="%1.1f%%",
                startangle=90, wedgeprops=dict(width=0.5))
    ax1.set_title(
        "COSMIC Gene Coverage\n(consequence-passing variants)",
        fontweight="bold"
    )

    # --- Panel B: top 20 genes bar chart ---
    top_n = 20
    top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    gnames  = [g for g, _ in top_genes]
    gcounts = [c for _, c in top_genes]
    gcolors = [
        PALETTE["tier1"] if gene_tier.get(g) == 1 else
        PALETTE["tier2"] if gene_tier.get(g) == 2 else
        PALETTE["neutral"]
        for g in gnames
    ]

    bars = ax2.barh(range(len(gnames)), gcounts, color=gcolors, edgecolor="white", height=0.65)
    for bar, count in zip(bars, gcounts):
        ax2.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                 f"  {count:,}", va="center", ha="left", fontsize=8)

    ax2.set_yticks(range(len(gnames)))
    ax2.set_yticklabels(gnames, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("Variant count")
    ax2.set_title(f"Top {top_n} Genes (consequence-passing)", fontweight="bold")
    ax2.set_xlim(0, max(gcounts) * 1.2 if gcounts else 1)

    legend_handles = [
        mpatches.Patch(color=PALETTE["tier1"],   label="COSMIC Tier 1"),
        mpatches.Patch(color=PALETTE["tier2"],   label="COSMIC Tier 2"),
        mpatches.Patch(color=PALETTE["neutral"], label="Non-COSMIC"),
    ]
    ax2.legend(handles=legend_handles, fontsize=8)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    return fig


# =============================================================================
# Output helpers
# =============================================================================

def save_figure(fig, report_dir, name):
    for fmt, subdir in [("png", "png"), ("svg", "svg")]:
        d = os.path.join(report_dir, subdir)
        os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, f"{name}.{fmt}"),
                    dpi=150 if fmt == "png" else None,
                    bbox_inches="tight")


def write_stats(funnel, report_dir, sample_name):
    rows = []
    prev = None
    for step, count in funnel.items():
        pct = round((1 - count / prev) * 100, 2) if prev else 0.0
        rows.append({"filter_step": step, "variants": count, "pct_removed_from_prev": pct})
        prev = count
    tsv_path = os.path.join(report_dir, f"{sample_name}.filter_qc_stats.tsv")
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    print("\n=== CUP WES -- Somatic Variant Filter QC Report ===\n")
    print(f"  Input VCF:        {args.vcf}")
    print(f"  COSMIC list:      {args.cosmic}")
    print(f"  Min VAF:          {args.min_vaf}")
    print(f"  Min depth:        {args.min_depth}x")
    print(f"  COSMIC tier:      {args.tier if args.tier else 'All'}")
    print(f"  Exclude filters:  {args.exclude_filters if args.exclude_filters else 'None'}\n")

    sample_name = os.path.basename(args.vcf).replace(".vcf.gz", "").replace(".vcf", "")
    tier_suffix = f".tier{args.tier}" if args.tier else ""
    report_dir  = os.path.join(args.outdir, "reports", "dependantfilter")
    os.makedirs(report_dir, exist_ok=True)

    # Load COSMIC
    print("Loading COSMIC Cancer Gene Census...")
    gene_tier = load_cosmic_genes(args.cosmic, args.tier)

    # Parse VEP CSQ field names
    vcf = cyvcf2.VCF(args.vcf)
    csq_fields = []
    for header_line in vcf.header_iter():
        info = header_line.info(extra=True)
        if info.get("ID") == "CSQ":
            desc = info.get("Description", "")
            if "Format:" in desc:
                csq_fields = desc.split("Format:")[1].strip().strip('"').split("|")
            break
    vcf.close()

    if not csq_fields:
        print("ERROR: Could not parse VEP CSQ field from VCF header.")
        sys.exit(1)
    print(f"  Found {len(csq_fields)} VEP annotation fields\n")

    # Single-pass metric collection
    print("Collecting metrics (single pass through VCF)...")
    metrics = collect_metrics(args, csq_fields, gene_tier)
    print("  Done.\n")

    print("=== Filtering Summary ===\n")
    for step, count in metrics["funnel"].items():
        print(f"  {step:<45} {count:>8,}")
    print()

    # Generate figures
    figure_defs = [
        ("01_funnel",               fig_funnel(metrics, args)),
        ("02_filter_flags",         fig_filter_flags(metrics, args)),
        ("03_vaf_distribution",     fig_vaf(metrics, args)),
        ("04_depth_distribution",   fig_depth(metrics, args)),
        ("05_vaf_vs_depth",         fig_vaf_depth_scatter(metrics, args)),
        ("06_consequence_types",    fig_consequences(metrics, args)),
        ("07_cosmic_genes",         fig_cosmic(metrics, args)),
    ]
    figures = [(name, fig) for name, fig in figure_defs if fig is not None]

    # Save individual PNG and SVG
    print("Saving figures...")
    for name, fig in figures:
        save_figure(fig, report_dir, name)
        print(f"  {name}")

    # Save combined PDF
    pdf_path = os.path.join(report_dir, f"{sample_name}{tier_suffix}.filter_qc.pdf")
    with PdfPages(pdf_path) as pdf:
        for _, fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    print(f"\n  PDF:       {pdf_path}")

    # Save stats TSV
    tsv_path = write_stats(metrics["funnel"], report_dir, f"{sample_name}{tier_suffix}")
    print(f"  Stats TSV: {tsv_path}\n")

    for _, fig in figures:
        plt.close(fig)

    print("Done.\n")


if __name__ == "__main__":
    main()
