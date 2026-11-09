#!/usr/bin/env python3
"""
CUP WES Pipeline -- CNV Segment Filter
=======================================
Filters CNVKit tumour-only call segments by:
  - Minimum number of probes per segment
  - Minimum log2 ratio deviation from diploid
  - Minimum segment weight
  - Optionally filters to COSMIC Cancer Gene Census genes (skipped if --cosmic not provided)
  - Outputs filtered TSV with one row per gene per segment
  - Optionally generates a QC report with filtering funnel and diagnostic figures (--report)

Usage:
    python3 dnaseq-filtercnv.py \
        --cns <S1.tumor_only.call.cns> \
        --outdir <output_dir> \
        [--cosmic <Census_all.csv>] \
        [--min_probes 10] \
        [--min_log2 0.3] \
        [--min_weight 0.0] \
        [--tier 1] \
        [--report]

Requirements:
    pip install pandas --break-system-packages
    pip install matplotlib numpy --break-system-packages  # only needed for --report
"""

import argparse
import os
import re
import sys

import pandas as pd

# =============================================================================
# Constants
# =============================================================================

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]

PALETTE = {
    "pass":    "#2ecc71",
    "fail":    "#e74c3c",
    "neutral": "#95a5a6",
    "blue":    "#3498db",
    "orange":  "#e67e22",
    "gain":    "#2ecc71",
    "loss":    "#e74c3c",
    "tier1":   "#2ecc71",
    "tier2":   "#3498db",
}

# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="CNV segment filter")
    parser.add_argument("--cns",        required=True, help="CNVKit tumour-only call file (S1.tumor_only.call.cns)")
    parser.add_argument("--cosmic",     default=None,  help="COSMIC Cancer Gene Census CSV (Census_all.csv). If omitted, gene filter is skipped.")
    parser.add_argument("--outdir",     required=True, help="Output directory")
    parser.add_argument("--min_probes", type=int,   default=10,  help="Minimum probes per segment (default: 10)")
    parser.add_argument("--min_log2",   type=float, default=0.3, help="Minimum absolute log2 ratio deviation (default: 0.3)")
    parser.add_argument("--min_weight", type=float, default=0.0, help="Minimum segment weight (default: 0.0)")
    parser.add_argument("--tier",       type=int,   default=None, help="Filter to COSMIC tier (1 or 2). Default: both")
    parser.add_argument("--report",     action="store_true",
                        help="Generate QC report with filtering funnel and diagnostic figures (requires matplotlib)")
    return parser.parse_args()

# =============================================================================
# Load COSMIC Cancer Gene Census -- returns gene -> tier dict
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
# Parse gene symbols from CNVKit gene column
# =============================================================================

def parse_gene_symbols(gene_field):
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
# Report figures
# =============================================================================

def fig_funnel(funnel, args):
    import matplotlib.pyplot as plt
    labels = list(funnel.keys())
    counts = list(funnel.values())

    fig, ax = plt.subplots(figsize=(10, 4))
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
    ax.set_xlabel("Segment count")
    ax.set_title("CNV Filtering Funnel — cumulative effect of each step", fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.3 if counts else 1)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_probes(segs, args):
    import matplotlib.pyplot as plt
    import numpy as np
    probes    = segs["probes"].values.astype(float)
    threshold = args.min_probes

    fig, ax = plt.subplots(figsize=(9, 5))
    cap  = max(float(np.percentile(probes, 99)) if len(probes) else threshold * 3, threshold * 3)
    bins = np.linspace(0, cap, 50)

    ax.hist(probes[probes < threshold],  bins=bins, color=PALETTE["fail"], alpha=0.75,
            label=f"Excluded  (probes < {threshold})")
    ax.hist(probes[probes >= threshold], bins=bins, color=PALETTE["pass"], alpha=0.75,
            label=f"Kept  (probes ≥ {threshold})")
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}")

    n_excl = int((probes < threshold).sum())
    n_kept = int((probes >= threshold).sum())
    ax.text(0.98, 0.97, f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Probes per segment")
    ax.set_ylabel("Segment count")
    ax.set_title("Probe Count Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_log2(segs, args):
    import matplotlib.pyplot as plt
    import numpy as np
    log2      = segs["log2"].values.astype(float)
    threshold = args.min_log2

    fig, ax = plt.subplots(figsize=(9, 5))
    bound = max(float(np.abs(log2).max()) * 1.05, threshold * 3) if len(log2) else threshold * 3
    bins  = np.linspace(-bound, bound, 61)

    near_diploid = np.abs(log2) < threshold
    gain_mask    = log2 >= threshold
    loss_mask    = log2 <= -threshold

    ax.hist(log2[near_diploid], bins=bins, color=PALETTE["neutral"], alpha=0.6,
            label="Near-diploid (excluded)")
    ax.hist(log2[gain_mask],    bins=bins, color=PALETTE["gain"],    alpha=0.75, label="Gain (kept)")
    ax.hist(log2[loss_mask],    bins=bins, color=PALETTE["loss"],    alpha=0.75, label="Loss (kept)")
    ax.axvline( threshold, color="black", linestyle="--", linewidth=1.5)
    ax.axvline(-threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold ±{threshold}")
    ax.axhspan(-threshold, threshold, alpha=0.07, color="grey")

    n_gain = int(gain_mask.sum())
    n_loss = int(loss_mask.sum())
    n_excl = int(near_diploid.sum())
    ax.text(0.98, 0.97, f"Gains:    {n_gain:,}\nLosses:   {n_loss:,}\nExcluded: {n_excl:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("log2 ratio")
    ax.set_ylabel("Segment count")
    ax.set_title("Log2 Ratio Distribution — exclusion band shaded", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_weight(segs, args):
    import matplotlib.pyplot as plt
    import numpy as np
    weights   = segs["weight"].values.astype(float)
    threshold = args.min_weight

    fig, ax = plt.subplots(figsize=(9, 5))
    upper = max(float(weights.max()), 1.0) if len(weights) else 1.0
    bins  = np.linspace(0, upper, 50)

    ax.hist(weights[weights < threshold],  bins=bins, color=PALETTE["fail"], alpha=0.75,
            label=f"Excluded  (weight < {threshold})")
    ax.hist(weights[weights >= threshold], bins=bins, color=PALETTE["pass"], alpha=0.75,
            label=f"Kept  (weight ≥ {threshold})")
    ax.axvline(threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}")

    n_excl = int((weights < threshold).sum())
    n_kept = int((weights >= threshold).sum())
    ax.text(0.98, 0.97, f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Segment weight")
    ax.set_ylabel("Segment count")
    ax.set_title("Segment Weight Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_probes_log2_scatter(segs, args):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if segs.empty:
        return None

    probes = segs["probes"].values.astype(float)
    log2   = segs["log2"].values.astype(float)
    stages = segs["stage"].values

    color_map = [
        PALETTE["pass"]    if s == "pass"        else
        PALETTE["fail"]    if s == "fail_probes" else
        PALETTE["neutral"] if s == "fail_log2"   else
        PALETTE["orange"]
        for s in stages
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(probes, log2, c=color_map, alpha=0.4, s=8, linewidths=0)
    ax.axvline(args.min_probes, color="black", linestyle="--", linewidth=1.2,
               label=f"min probes = {args.min_probes}")
    ax.axhline( args.min_log2, color="black", linestyle="-.", linewidth=1.2,
               label=f"|log2| = {args.min_log2}")
    ax.axhline(-args.min_log2, color="black", linestyle="-.", linewidth=1.2)
    ax.axhline(0, color="grey", linewidth=0.7, linestyle=":")

    ax.legend(handles=[
        mpatches.Patch(color=PALETTE["pass"],    label="Pass all filters"),
        mpatches.Patch(color=PALETTE["fail"],    label="Fail probes"),
        mpatches.Patch(color=PALETTE["neutral"], label="Fail log2"),
        mpatches.Patch(color=PALETTE["orange"],  label="Fail weight / COSMIC"),
    ], fontsize=8)
    ax.set_xlabel("Probes per segment")
    ax.set_ylabel("log2 ratio")
    ax.set_title("Probes vs Log2 Ratio — joint threshold view", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_cosmic(gene_counts, gene_tier):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if not gene_counts:
        return None

    n_tier1 = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if gene_tier.get(g) == 1)
    n_tier2 = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if gene_tier.get(g) == 2)
    n_non   = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if g not in gene_tier)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    sizes        = [n_tier1, n_tier2, n_non]
    clrs         = [PALETTE["tier1"], PALETTE["tier2"], PALETTE["neutral"]]
    donut_labels = [f"COSMIC Tier 1\n({n_tier1:,})", f"COSMIC Tier 2\n({n_tier2:,})", f"Non-COSMIC\n({n_non:,})"]
    non_zero = [(s, lbl, c) for s, lbl, c in zip(sizes, donut_labels, clrs) if s > 0]
    if non_zero:
        s_, lbl_, c_ = zip(*non_zero)
        ax1.pie(s_, labels=lbl_, colors=c_, autopct="%1.1f%%",
                startangle=90, wedgeprops=dict(width=0.5))
    ax1.set_title("COSMIC Gene Coverage\n(weight-passing segments)", fontweight="bold")

    top_genes = sorted({g: v["gain"] + v["loss"] for g, v in gene_counts.items()}.items(),
                       key=lambda x: x[1], reverse=True)[:20]
    gnames = [g for g, _ in top_genes]
    gains  = [gene_counts[g]["gain"] for g in gnames]
    losses = [gene_counts[g]["loss"] for g in gnames]

    y = range(len(gnames))
    ax2.barh(list(y), gains,  color=PALETTE["gain"], edgecolor="white", height=0.65, label="Gain")
    ax2.barh(list(y), losses, left=gains, color=PALETTE["loss"], edgecolor="white", height=0.65, label="Loss")

    for i, (g_count, l_count, gname) in enumerate(zip(gains, losses, gnames)):
        total = g_count + l_count
        tier  = gene_tier.get(gname)
        marker = " ★" if tier == 1 else " ◆" if tier == 2 else ""
        ax2.text(total, i, f"  {total:,}{marker}", va="center", ha="left", fontsize=8)

    ax2.set_yticks(list(y))
    ax2.set_yticklabels(gnames, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("Segment count")
    ax2.set_title("Top 20 Genes — Gain/Loss breakdown\n(★ Tier 1  ◆ Tier 2)", fontweight="bold")
    ax2.set_xlim(0, max(g + l for g, l in zip(gains, losses)) * 1.25 if gains else 1)
    ax2.legend(fontsize=8)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def fig_genome(segs, args):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    if segs.empty:
        return None

    df = segs.copy()
    chrom_order     = {c: i for i, c in enumerate(CHR_ORDER)}
    df["chrom_clean"] = df["chrom"].str.replace("chr", "", regex=False)
    df["chrom_idx"]   = df["chrom_clean"].map(chrom_order)
    df = df.dropna(subset=["chrom_idx"])
    df["chrom_idx"] = df["chrom_idx"].astype(int)

    chrom_max = df.groupby("chrom_idx")["midpt"].transform("max").clip(lower=1)
    df["x"] = df["chrom_idx"] + df["midpt"] / chrom_max

    fig, ax = plt.subplots(figsize=(14, 5))

    excl = df[df["stage"] != "pass"]
    if not excl.empty:
        ax.scatter(excl["x"], excl["log2"], c=PALETTE["neutral"], alpha=0.2, s=5, linewidths=0)

    passed = df[df["stage"] == "pass"]
    if not passed.empty:
        gains  = passed[passed["log2"] > 0]
        losses = passed[passed["log2"] <= 0]
        if not gains.empty:
            ax.scatter(gains["x"],  gains["log2"],  c=PALETTE["gain"], alpha=0.7, s=12, linewidths=0, label="Gain")
        if not losses.empty:
            ax.scatter(losses["x"], losses["log2"], c=PALETTE["loss"], alpha=0.7, s=12, linewidths=0, label="Loss")

        cosmic_hits = passed[passed["cosmic_genes"].apply(len) > 0]
        if not cosmic_hits.empty:
            ax.scatter(cosmic_hits["x"], cosmic_hits["log2"],
                       facecolors="none", edgecolors="black",
                       s=45, linewidths=0.8, zorder=5, label="COSMIC gene")

    ax.axhline( args.min_log2, color="black", linestyle="--", linewidth=0.8, alpha=0.45)
    ax.axhline(-args.min_log2, color="black", linestyle="--", linewidth=0.8, alpha=0.45)
    ax.axhline(0, color="grey", linewidth=0.6, linestyle=":")

    chrs_present = sorted(df["chrom_idx"].unique())
    for idx in chrs_present:
        ax.axvline(idx, color="lightgrey", linewidth=0.5, zorder=0)

    ax.set_xticks([i + 0.5 for i in chrs_present])
    ax.set_xticklabels([CHR_ORDER[i] for i in chrs_present], fontsize=7)

    handles, lbls = ax.get_legend_handles_labels()
    handles.append(mpatches.Patch(color=PALETTE["neutral"], alpha=0.5, label="Excluded"))
    ax.legend(handles=handles, fontsize=8, loc="upper right")

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("log2 ratio")
    ax.set_title("Genome-wide CNV Profile", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
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
        rows.append({"filter_step": step, "segments": count, "pct_removed_from_prev": pct})
        prev = count
    tsv_path = os.path.join(report_dir, f"{sample_name}.cnv_filter_qc_stats.tsv")
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


def generate_report(metrics, args, report_dir, sample_name, tier_suffix):
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
    segs = metrics["segs"]

    figure_defs = [
        ("01_funnel",            fig_funnel(metrics["funnel"], args)),
        ("02_probe_distribution",fig_probes(segs, args)),
        ("03_log2_distribution", fig_log2(segs, args)),
        ("04_weight_distribution",fig_weight(segs, args)),
        ("05_probes_vs_log2",    fig_probes_log2_scatter(segs, args)),
        ("07_genome_profile",    fig_genome(segs, args)),
    ]
    if metrics["gene_tier"]:
        figure_defs.insert(5, ("06_cosmic_genes", fig_cosmic(metrics["gene_counts"], metrics["gene_tier"])))

    figures = [(name, fig) for name, fig in figure_defs if fig is not None]

    for name, fig in figures:
        save_figure(fig, report_dir, name)
        print(f"  {name}")

    pdf_path = os.path.join(report_dir, f"{sample_name}{tier_suffix}.cnv_filter_qc.pdf")
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

    print("\n=== CUP WES -- CNV Segment Filter ===\n")
    print(f"  Input CNS:    {args.cns}")
    print(f"  COSMIC list:  {args.cosmic if args.cosmic else 'None (gene filter disabled)'}")
    print(f"  Min probes:   {args.min_probes}")
    print(f"  Min |log2|:   {args.min_log2}")
    print(f"  Min weight:   {args.min_weight}")
    print(f"  COSMIC tier:  {args.tier if args.tier else 'All'}\n")

    # Load COSMIC genes (optional)
    if args.cosmic:
        print("Loading COSMIC Cancer Gene Census...")
        gene_tier = load_cosmic_genes(args.cosmic, tier=args.tier)
    else:
        gene_tier = None
        print("No COSMIC file provided — gene filter will be skipped.\n")

    # Load CNS file
    print("Loading CNS file...")
    df = pd.read_csv(args.cns, sep="\t")
    print(f"  Loaded {len(df)} segments\n")

    # Counters
    total       = len(df)
    pass_probes = 0
    pass_log2   = 0
    pass_weight = 0
    pass_gene   = 0
    kept        = 0
    rows        = []

    # Report metric collectors
    seg_records = []
    gene_counts = {}  # gene -> {"gain": n, "loss": n}

    print("Filtering segments...\n")

    for _, row in df.iterrows():
        probes = int(row["probes"])
        log2   = float(row["log2"])
        weight = float(row["weight"])
        chrom  = str(row["chromosome"]).replace("chr", "")
        midpt  = (float(row["start"]) + float(row["end"])) / 2.0

        # Evaluate each filter
        pv_probes = probes >= args.min_probes
        pv_log2   = abs(log2) >= args.min_log2
        pv_weight = weight >= args.min_weight

        gene_symbols = parse_gene_symbols(row["gene"])
        if gene_tier is not None:
            cosmic_genes_in_seg = [g for g in gene_symbols if g in gene_tier]
            pv_gene = len(cosmic_genes_in_seg) > 0
        else:
            cosmic_genes_in_seg = []
            pv_gene = True

        # Determine stage for report
        if not pv_probes:
            stage = "fail_probes"
        elif not pv_log2:
            stage = "fail_log2"
        elif not pv_weight:
            stage = "fail_weight"
        elif not pv_gene:
            stage = "fail_gene"
        else:
            stage = "pass"

        seg_records.append({
            "chrom":        chrom,
            "midpt":        midpt,
            "log2":         log2,
            "probes":       probes,
            "weight":       weight,
            "stage":        stage,
            "cosmic_genes": cosmic_genes_in_seg,
            "call":         "GAIN" if log2 > 0 else "LOSS",
        })

        # Apply filter cascade
        if not pv_probes:
            continue
        pass_probes += 1

        if not pv_log2:
            continue
        pass_log2 += 1

        if not pv_weight:
            continue
        pass_weight += 1

        # Collect gene counts for report (all genes, not just COSMIC)
        call_key = "gain" if log2 > 0 else "loss"
        for g in gene_symbols:
            gene_counts.setdefault(g, {"gain": 0, "loss": 0})
            gene_counts[g][call_key] += 1

        if not pv_gene:
            continue
        pass_gene += 1

        kept += 1

        # Determine which genes to write — COSMIC-filtered or all
        output_genes = cosmic_genes_in_seg if gene_tier is not None else gene_symbols
        if not output_genes:
            continue

        call_type = "GAIN" if log2 > 0 else "LOSS"
        for gene in output_genes:
            rows.append({
                "CHROM":  row["chromosome"],
                "START":  row["start"],
                "END":    row["end"],
                "GENE":   gene,
                "LOG2":   round(log2, 4),
                "CN":     int(row["cn"]) if not pd.isna(row["cn"]) else "",
                "CALL":   call_type,
                "DEPTH":  round(row["depth"], 1),
                "PROBES": probes,
                "WEIGHT": round(weight, 4),
                "CI_HI":  round(row["ci_hi"], 4),
                "CI_LO":  round(row["ci_lo"], 4),
            })

    # Write output TSV
    sample_name = os.path.basename(args.cns).replace(".tumor_only.call.cns", "")
    tier_suffix = f".tier{args.tier}" if args.tier else ""
    out_tsv     = os.path.join(args.outdir, f"{sample_name}.cnv{tier_suffix}.filtered.tsv")

    out_df = pd.DataFrame(rows)
    if not out_df.empty:
        out_df = out_df.sort_values(["CHROM", "START"])
    out_df.to_csv(out_tsv, sep="\t", index=False)

    # ==========================================================================
    # Summary
    # ==========================================================================

    funnel = {
        "Total segments":                    total,
        f"Pass probes >= {args.min_probes}": pass_probes,
        f"Pass |log2| >= {args.min_log2}":   pass_log2,
        f"Pass weight >= {args.min_weight}": pass_weight,
        "Pass COSMIC gene filter":           pass_gene,
        "Final segments kept":               kept,
    }

    print("=== Filtering Summary ===\n")
    print(f"  Total segments input:            {total:>8,}")
    print(f"  Pass probes >= {args.min_probes}:               {pass_probes:>8,}")
    print(f"  Pass |log2| >= {args.min_log2}:              {pass_log2:>8,}")
    print(f"  Pass weight >= {args.min_weight}:              {pass_weight:>8,}")
    print(f"  Pass COSMIC gene filter:         {pass_gene:>8,}")
    print(f"  Final segments kept:             {kept:>8,}")
    print(f"  Output rows (one per gene):      {len(rows):>8,}")
    print(f"\n  Output TSV: {out_tsv}\n")

    # ==========================================================================
    # QC Report (optional)
    # ==========================================================================

    if args.report:
        report_dir = os.path.join(args.outdir, "reports", "filter_qc")
        os.makedirs(report_dir, exist_ok=True)
        metrics = {
            "funnel":      funnel,
            "segs":        pd.DataFrame(seg_records),
            "gene_counts": gene_counts,
            "gene_tier":   gene_tier or {},
        }
        generate_report(metrics, args, report_dir, sample_name, tier_suffix)


if __name__ == "__main__":
    main()
