#!/usr/bin/env python3
"""
CUP WES Pipeline -- Post-VEP Somatic Variant Filter
=====================================================
Filters VEP-annotated Mutect2 VCF by:
  - Minimum VAF (default 5%)
  - Minimum depth (default 20x)
  - Consequence filter: keeps protein-altering consequences by default (configurable or disableable)
  - Optionally filters to COSMIC Cancer Gene Census genes (skipped if --cosmic not provided)
  - Excludes variants with specified Mutect2 filter flags
  - Outputs filtered VCF and TSV table
  - Optionally generates a QC report with filtering funnel and diagnostic figures (--report)

Usage:
    python3 dnaseq-filtervariants.py \
        --vcf <vep_annotated.vcf.gz> \
        --outdir <output_dir> \
        [--cosmic <Census_all.csv>] \
        [--min_vaf 0.05] \
        [--min_depth 20] \
        [--tier 1] \
        [--exclude_filters germline panel_of_normals strand_bias base_qual weak_evidence] \
        [--consequences missense_variant stop_gained frameshift_variant ...] \
        [--no_consequence_filter] \
        [--report]

Requirements:
    pip install cyvcf2 pandas --break-system-packages
    pip install matplotlib numpy --break-system-packages  # only needed for --report
"""

import argparse
import collections
import os
import sys

import cyvcf2
import pandas as pd

# =============================================================================
# Default consequences to keep (protein-altering, high and moderate impact)
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
# Argument parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Post-VEP somatic variant filter")
    parser.add_argument("--vcf",                   required=True, help="VEP annotated VCF (gz or uncompressed)")
    parser.add_argument("--cosmic",                default=None,  help="COSMIC Cancer Gene Census CSV (Census_all.csv). If omitted, gene filter is skipped.")
    parser.add_argument("--outdir",                required=True, help="Output directory")
    parser.add_argument("--min_vaf",               type=float, default=0.05, help="Minimum VAF (default: 0.05)")
    parser.add_argument("--min_depth",             type=int,   default=20,   help="Minimum depth (default: 20)")
    parser.add_argument("--tier",                  type=int,   default=None, help="Filter to COSMIC tier (1 or 2). Default: both")
    parser.add_argument("--exclude_filters",       nargs="+",  default=[],
                        help="Mutect2 filter flags to exclude. E.g: germline panel_of_normals strand_bias base_qual weak_evidence clustered_events haplotype slippage multiallelic")
    parser.add_argument("--consequences",          nargs="+",  default=None,
                        help="Custom list of VEP consequences to keep. Overrides the default set. E.g: missense_variant stop_gained frameshift_variant")
    parser.add_argument("--no_consequence_filter", action="store_true",
                        help="Disable consequence filtering entirely — keep all consequences")
    parser.add_argument("--report",                action="store_true",
                        help="Generate QC report with filtering funnel and diagnostic figures (requires matplotlib)")
    return parser.parse_args()

# =============================================================================
# Load COSMIC Cancer Gene Census — returns gene -> tier dict
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
# Parse VEP CSQ field
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
# Report figures
# =============================================================================

def fig_funnel(funnel, args):
    import matplotlib.pyplot as plt
    labels = list(funnel.keys())
    counts = list(funnel.values())

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [PALETTE["blue"]] + [PALETTE["pass"]] * (len(labels) - 2) + [PALETTE["tier1"]]
    bars = ax.barh(range(len(labels)), counts, color=colors, edgecolor="white", height=0.6)

    for i, (bar, count) in enumerate(zip(bars, counts)):
        prev  = counts[i - 1] if i > 0 else count
        pct   = (1 - count / prev) * 100 if (i > 0 and prev > 0) else 0
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


def fig_filter_flags(flag_counts, args):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    exclude_set   = set(args.exclude_filters)
    sorted_flags  = sorted(flag_counts.items(), key=lambda x: x[1], reverse=True)
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
        mpatches.Patch(color=PALETTE["neutral"], label="Not excluded"),
    ]
    if exclude_set:
        ax.legend(handles=legend_handles, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_vaf(all_vafs, args):
    import matplotlib.pyplot as plt
    import numpy as np
    vafs      = np.array(all_vafs, dtype=float)
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
    ax.text(0.98, 0.97, f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Variant Allele Frequency (VAF)")
    ax.set_ylabel("Variant count")
    ax.set_title("VAF Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_depth(all_depths, args):
    import matplotlib.pyplot as plt
    import numpy as np
    depths    = np.array(all_depths, dtype=float)
    threshold = args.min_depth

    fig, ax = plt.subplots(figsize=(9, 5))
    cap  = float(np.percentile(depths, 99)) if len(depths) > 0 else 500.0
    cap  = max(cap, threshold * 2)
    bins = np.linspace(0, cap, 60)
    ax.hist(depths[depths < threshold],  bins=bins, color=PALETTE["fail"], alpha=0.75,
            label=f"Excluded  (depth < {threshold}x)")
    ax.hist(depths[depths >= threshold], bins=bins, color=PALETTE["pass"], alpha=0.75,
            label=f"Kept  (depth ≥ {threshold}x)")
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}x")

    n_excl = int((depths < threshold).sum())
    n_kept = int((depths >= threshold).sum())
    ax.text(0.98, 0.97, f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Read Depth")
    ax.set_ylabel("Variant count")
    ax.set_title("Read Depth Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_vaf_depth_scatter(scatter, args):
    import matplotlib.pyplot as plt
    import numpy as np
    if not scatter:
        return None

    vafs   = np.array([x[0] for x in scatter])
    depths = np.array([x[1] for x in scatter], dtype=float)
    pv_arr = np.array([x[2] for x in scatter])
    pd_arr = np.array([x[3] for x in scatter])

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
        "both":       int(( pv_arr &  pd_arr).sum()),
        "fail_vaf":   int((~pv_arr &  pd_arr).sum()),
        "fail_depth": int(( pv_arr & ~pd_arr).sum()),
        "both_fail":  int((~pv_arr & ~pd_arr).sum()),
    }
    kw = dict(transform=ax.transAxes, fontsize=8, fontweight="bold")
    ax.text(0.98, 0.98, f"Both pass\n{q['both']:,}",        ha="right", va="top",    color=PALETTE["pass"],   **kw)
    ax.text(0.02, 0.98, f"Fail VAF\n{q['fail_vaf']:,}",     ha="left",  va="top",    color=PALETTE["orange"], **kw)
    ax.text(0.98, 0.02, f"Fail depth\n{q['fail_depth']:,}", ha="right", va="bottom", color=PALETTE["blue"],   **kw)
    ax.text(0.02, 0.02, f"Both fail\n{q['both_fail']:,}",   ha="left",  va="bottom", color=PALETTE["fail"],   **kw)

    ax.set_xlabel("VAF")
    ax.set_ylabel("Read Depth")
    ax.set_title("VAF vs Read Depth — joint threshold view", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_consequences(csq_counts, keep_consequences):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if not csq_counts:
        return None

    sorted_csq = sorted(csq_counts.items(), key=lambda x: x[1], reverse=True)
    labels = [c for c, _ in sorted_csq]
    counts = [n for _, n in sorted_csq]

    if keep_consequences is None:
        colors     = [PALETTE["pass"]] * len(labels)
        kept_label = "Kept (all consequences)"
    else:
        colors     = [PALETTE["pass"] if c in keep_consequences else PALETTE["fail"] for c in labels]
        kept_label = "Kept"

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
        "VEP Consequence Distribution\n(variants passing flag / VAF / depth filters)",
        fontweight="bold"
    )
    ax.set_xlim(0, max(counts) * 1.2)
    legend_handles = [mpatches.Patch(color=PALETTE["pass"], label=kept_label)]
    if keep_consequences is not None:
        legend_handles.append(mpatches.Patch(color=PALETTE["fail"], label="Excluded"))
    ax.legend(handles=legend_handles)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_cosmic(gene_counts, gene_tier):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if not gene_counts:
        return None

    n_tier1 = sum(c for g, c in gene_counts.items() if gene_tier.get(g) == 1)
    n_tier2 = sum(c for g, c in gene_counts.items() if gene_tier.get(g) == 2)
    n_non   = sum(c for g, c in gene_counts.items() if g not in gene_tier)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    sizes        = [n_tier1, n_tier2, n_non]
    clrs         = [PALETTE["tier1"], PALETTE["tier2"], PALETTE["neutral"]]
    donut_labels = [f"COSMIC Tier 1\n({n_tier1:,})", f"COSMIC Tier 2\n({n_tier2:,})", f"Non-COSMIC\n({n_non:,})"]
    non_zero = [(s, lbl, c) for s, lbl, c in zip(sizes, donut_labels, clrs) if s > 0]
    if non_zero:
        s_, lbl_, c_ = zip(*non_zero)
        ax1.pie(s_, labels=lbl_, colors=c_, autopct="%1.1f%%",
                startangle=90, wedgeprops=dict(width=0.5))
    ax1.set_title("COSMIC Gene Coverage\n(consequence-passing variants)", fontweight="bold")

    top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:20]
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
    ax2.set_title("Top 20 Genes (consequence-passing)", fontweight="bold")
    ax2.set_xlim(0, max(gcounts) * 1.2 if gcounts else 1)
    ax2.legend(handles=[
        mpatches.Patch(color=PALETTE["tier1"],   label="COSMIC Tier 1"),
        mpatches.Patch(color=PALETTE["tier2"],   label="COSMIC Tier 2"),
        mpatches.Patch(color=PALETTE["neutral"], label="Non-COSMIC"),
    ], fontsize=8)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


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


def generate_report(metrics, args, keep_consequences, report_dir, sample_name, tier_suffix):
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not installed — skipping report.")
        print("         Run: pip install matplotlib numpy --break-system-packages")
        return

    print("Generating QC report...")

    figure_defs = [
        ("01_funnel",             fig_funnel(metrics["funnel"], args)),
        ("02_filter_flags",       fig_filter_flags(metrics["flag_counts"], args)),
        ("03_vaf_distribution",   fig_vaf(metrics["all_vafs"], args)),
        ("04_depth_distribution", fig_depth(metrics["all_depths"], args)),
        ("05_vaf_vs_depth",       fig_vaf_depth_scatter(metrics["scatter"], args)),
        ("06_consequence_types",  fig_consequences(metrics["csq_counts"], keep_consequences)),
    ]
    if metrics["gene_tier"]:
        figure_defs.append(("07_cosmic_genes", fig_cosmic(metrics["gene_counts"], metrics["gene_tier"])))

    figures = [(name, fig) for name, fig in figure_defs if fig is not None]

    for name, fig in figures:
        save_figure(fig, report_dir, name)
        print(f"  {name}")

    pdf_path = os.path.join(report_dir, f"{sample_name}{tier_suffix}.filter_qc.pdf")
    with PdfPages(pdf_path) as pdf:
        for _, fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    print(f"\n  PDF:       {pdf_path}")

    tsv_path = write_stats(metrics["funnel"], report_dir, f"{sample_name}{tier_suffix}")
    print(f"  Stats TSV: {tsv_path}\n")

    for _, fig in figures:
        plt.close(fig)

# =============================================================================
# Main filtering logic
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Resolve consequence filter
    if args.no_consequence_filter:
        keep_consequences = None
        csq_label = "None (all consequences kept)"
    elif args.consequences:
        keep_consequences = set(args.consequences)
        csq_label = "custom"
    else:
        keep_consequences = KEEP_CONSEQUENCES
        csq_label = "default"

    print("\n=== CUP WES -- Somatic Variant Filter ===\n")
    print(f"  Input VCF:          {args.vcf}")
    print(f"  COSMIC list:        {args.cosmic if args.cosmic else 'None (gene filter disabled)'}")
    print(f"  Min VAF:            {args.min_vaf}")
    print(f"  Min depth:          {args.min_depth}")
    print(f"  COSMIC tier:        {args.tier if args.tier else 'All'}")
    print(f"  Exclude filters:    {args.exclude_filters if args.exclude_filters else 'None'}")
    print(f"  Consequence filter: {csq_label}")
    if keep_consequences:
        for c in sorted(keep_consequences):
            print(f"    - {c}")
    print()

    # Load COSMIC genes (optional)
    if args.cosmic:
        print("Loading COSMIC Cancer Gene Census...")
        gene_tier = load_cosmic_genes(args.cosmic, tier=args.tier)
    else:
        gene_tier = None
        print("No COSMIC file provided — gene filter will be skipped.\n")

    # Open VCF and parse VEP CSQ field names from header
    vcf = cyvcf2.VCF(args.vcf)
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

    out_vcf = cyvcf2.Writer(out_vcf_path, vcf)

    # Funnel counters
    total             = 0
    pass_filter_flags = 0
    pass_vaf          = 0
    pass_depth        = 0
    pass_csq          = 0
    pass_gene         = 0
    kept              = 0

    # Report metric collectors
    all_vafs    = []
    all_depths  = []
    scatter     = []
    flag_counts = collections.Counter()
    csq_counts  = collections.Counter()
    gene_counts = collections.Counter()

    rows = []

    print("Filtering variants...\n")

    for variant in vcf:
        total += 1

        # --- Collect Mutect2 filter flags ---
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

        # --- Mutect2 filter flag exclusion ---
        if args.exclude_filters and flags.intersection(set(args.exclude_filters)):
            continue
        pass_filter_flags += 1

        # --- Read VAF and depth ---
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
        pv  = af >= args.min_vaf
        pd_ = dp >= args.min_depth
        scatter.append((af, dp, pv, pd_))

        if not pv:
            continue
        pass_vaf += 1

        if not pd_:
            continue
        pass_depth += 1

        # --- Parse VEP annotations ---
        annotations = parse_csq(variant, csq_fields)
        if not annotations:
            continue

        # Collect consequence counts before filtering
        for ann in annotations:
            for term in ann.get("Consequence", "").split("&"):
                t = term.strip()
                if t:
                    csq_counts[t] += 1

        # --- Consequence and gene filter ---
        csq_passed       = False
        gene_passed      = False
        kept_annotations = []
        for ann in annotations:
            consequences = set(ann.get("Consequence", "").split("&"))
            gene         = ann.get("SYMBOL", "").strip()

            if keep_consequences is not None and not consequences.intersection(keep_consequences):
                continue
            csq_passed = True

            if gene_tier is not None and gene not in gene_tier:
                continue
            gene_passed = True

            kept_annotations.append(ann)

        if csq_passed:
            pass_csq += 1
        if gene_passed:
            pass_gene += 1

        if not kept_annotations:
            continue

        kept += 1
        out_vcf.write_record(variant)

        # Collect gene counts for report
        genes_seen = set()
        for ann in kept_annotations:
            g = ann.get("SYMBOL", "").strip()
            if g and g not in genes_seen:
                gene_counts[g] += 1
                genes_seen.add(g)

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
    pd.DataFrame(rows).to_csv(out_tsv_path, sep="\t", index=False)

    # ==========================================================================
    # Summary
    # ==========================================================================

    funnel = {
        "Total input":                     total,
        "Pass filter flag exclusion":      pass_filter_flags,
        f"Pass VAF >= {args.min_vaf}":     pass_vaf,
        f"Pass depth >= {args.min_depth}": pass_depth,
        "Pass consequence filter":         pass_csq,
        "Pass COSMIC gene filter":         pass_gene,
        "Final variants kept":             kept,
    }

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

    # ==========================================================================
    # QC Report (optional)
    # ==========================================================================

    if args.report:
        report_dir = os.path.join(args.outdir, "reports", "filter_qc")
        os.makedirs(report_dir, exist_ok=True)
        metrics = {
            "funnel":      funnel,
            "all_vafs":    all_vafs,
            "all_depths":  all_depths,
            "scatter":     scatter,
            "flag_counts": flag_counts,
            "csq_counts":  csq_counts,
            "gene_counts": gene_counts,
            "gene_tier":   gene_tier or {},
        }
        generate_report(metrics, args, keep_consequences, report_dir, sample_name, tier_suffix)


if __name__ == "__main__":
    main()
