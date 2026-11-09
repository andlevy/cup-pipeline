#!/usr/bin/env python3
"""
CUP WES Pipeline -- CNV Segment Filter QC Report
=================================================
Generates diagnostic figures for each filter step applied by
dnaseq-filtercnv.py. Accepts identical arguments so it can be
run alongside the filter script with the same command-line call.

Outputs (under {outdir}/reports/dependantfilter/):
  {sample}.cnv_filter_qc.pdf          -- all figures, one per page
  {sample}.cnv_filter_qc_stats.tsv    -- per-step segment counts
  png/{01..07}_*.png                  -- individual figures
  svg/{01..07}_*.svg                  -- individual figures

Usage:
    python3 dnaseq-filtercnv-report.py \
        --cns <S1.tumor_only.call.cns> \
        --cosmic <Census_all.csv> \
        --outdir <results_dir> \
        [--min_probes 10] \
        [--min_log2 0.3] \
        [--min_weight 0.0] \
        [--tier 1]

Requirements:
    pip install pandas matplotlib numpy --break-system-packages
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


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
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="CNV segment filter QC report")
    p.add_argument("--cns",        required=True,  help="CNVKit tumour-only call file (.cns)")
    p.add_argument("--cosmic",     required=True,  help="COSMIC Cancer Gene Census CSV (Census_all.csv)")
    p.add_argument("--outdir",     required=True,  help="Results root dir (figures go in reports/dependantfilter/)")
    p.add_argument("--min_probes", type=int,   default=10,  help="Minimum probes per segment (default: 10)")
    p.add_argument("--min_log2",   type=float, default=0.3, help="Minimum absolute log2 ratio deviation (default: 0.3)")
    p.add_argument("--min_weight", type=float, default=0.0, help="Minimum segment weight (default: 0.0)")
    p.add_argument("--tier",       type=int,   default=None, help="COSMIC tier filter (1 or 2; default: both)")
    return p.parse_args()


# =============================================================================
# COSMIC loader -- gene -> tier dict
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
# Gene symbol parser (mirrors dnaseq-filtercnv.py)
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
# Metric collection
# =============================================================================

def collect_metrics(args, gene_tier):
    """
    Reads the CNS file once. Applies filter steps in order and records
    per-segment data for every figure.
    """
    df = pd.read_csv(args.cns, sep="\t")
    total = len(df)

    records = []
    gene_counts = {}   # gene -> {"gain": n, "loss": n}  (weight-passing segments)

    n_pass_probes = 0
    n_pass_log2   = 0
    n_pass_weight = 0
    n_pass_gene   = 0

    for _, row in df.iterrows():
        probes = int(row["probes"])
        log2   = float(row["log2"])
        weight = float(row["weight"])
        chrom  = str(row["chromosome"]).replace("chr", "")
        midpt  = (float(row["start"]) + float(row["end"])) / 2.0

        pass_probes = probes >= args.min_probes
        pass_log2   = abs(log2) >= args.min_log2
        pass_weight = weight >= args.min_weight

        gene_symbols = parse_gene_symbols(row["gene"])
        cosmic_genes = [g for g in gene_symbols if g in gene_tier]
        pass_gene    = len(cosmic_genes) > 0

        if not pass_probes:
            stage = "fail_probes"
        elif not pass_log2:
            stage = "fail_log2"
        elif not pass_weight:
            stage = "fail_weight"
        elif not pass_gene:
            stage = "fail_gene"
        else:
            stage = "pass"

        records.append({
            "chrom":        chrom,
            "midpt":        midpt,
            "log2":         log2,
            "probes":       probes,
            "weight":       weight,
            "stage":        stage,
            "cosmic_genes": cosmic_genes,
            "call":         "GAIN" if log2 > 0 else "LOSS",
        })

        if pass_probes:
            n_pass_probes += 1
            if pass_log2:
                n_pass_log2 += 1
                if pass_weight:
                    n_pass_weight += 1
                    call_key = "gain" if log2 > 0 else "loss"
                    for g in gene_symbols:
                        gene_counts.setdefault(g, {"gain": 0, "loss": 0})
                        gene_counts[g][call_key] += 1
                    if pass_gene:
                        n_pass_gene += 1

    funnel = {
        "Total segments":                    total,
        f"Pass probes ≥ {args.min_probes}":  n_pass_probes,
        f"Pass |log2| ≥ {args.min_log2}":    n_pass_log2,
        f"Pass weight ≥ {args.min_weight}":  n_pass_weight,
        "Pass COSMIC gene filter":           n_pass_gene,
    }

    return {
        "funnel":      funnel,
        "segs":        pd.DataFrame(records),
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

    fig, ax = plt.subplots(figsize=(10, 4))
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
    ax.set_xlabel("Segment count")
    ax.set_title("CNV Filtering Funnel — cumulative effect of each step", fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.3 if counts else 1)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 2 -- Probe count distribution
# =============================================================================

def fig_probes(metrics, args):
    probes    = metrics["segs"]["probes"].values.astype(float)
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
    ax.text(0.98, 0.97,
            f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Probes per segment")
    ax.set_ylabel("Segment count")
    ax.set_title("Probe Count Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 3 -- Log2 ratio distribution
# =============================================================================

def fig_log2(metrics, args):
    log2      = metrics["segs"]["log2"].values.astype(float)
    threshold = args.min_log2

    fig, ax = plt.subplots(figsize=(9, 5))
    bound = max(float(np.abs(log2).max()) * 1.05, threshold * 3) if len(log2) else threshold * 3
    bins  = np.linspace(-bound, bound, 61)

    near_diploid = np.abs(log2) < threshold
    gain_mask    = log2 >= threshold
    loss_mask    = log2 <= -threshold

    ax.hist(log2[near_diploid], bins=bins, color=PALETTE["neutral"], alpha=0.6,
            label="Near-diploid (excluded)")
    ax.hist(log2[gain_mask],    bins=bins, color=PALETTE["gain"],    alpha=0.75,
            label="Gain (kept)")
    ax.hist(log2[loss_mask],    bins=bins, color=PALETTE["loss"],    alpha=0.75,
            label="Loss (kept)")

    ax.axvline( threshold, color="black", linestyle="--", linewidth=1.5)
    ax.axvline(-threshold, color="black", linestyle="--", linewidth=1.5,
               label=f"Threshold ±{threshold}")
    ax.axvspan(-threshold, threshold, alpha=0.07, color="grey")

    n_gain = int(gain_mask.sum())
    n_loss = int(loss_mask.sum())
    n_excl = int(near_diploid.sum())
    ax.text(0.98, 0.97,
            f"Gains:    {n_gain:,}\nLosses:   {n_loss:,}\nExcluded: {n_excl:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("log2 ratio")
    ax.set_ylabel("Segment count")
    ax.set_title("Log2 Ratio Distribution — exclusion band shaded", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 4 -- Segment weight distribution
# =============================================================================

def fig_weight(metrics, args):
    weights   = metrics["segs"]["weight"].values.astype(float)
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
    ax.text(0.98, 0.97,
            f"Excluded: {n_excl:,}\nKept:     {n_kept:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="grey", alpha=0.8))

    ax.set_xlabel("Segment weight")
    ax.set_ylabel("Segment count")
    ax.set_title("Segment Weight Distribution", fontweight="bold")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 5 -- Probes vs log2 scatter (joint threshold view)
# =============================================================================

def fig_probes_log2_scatter(metrics, args):
    segs = metrics["segs"]
    if segs.empty:
        return None

    probes = segs["probes"].values.astype(float)
    log2   = segs["log2"].values.astype(float)
    stages = segs["stage"].values

    color_map = [
        PALETTE["pass"]    if s == "pass"         else
        PALETTE["fail"]    if s == "fail_probes"  else
        PALETTE["neutral"] if s == "fail_log2"    else
        PALETTE["orange"]
        for s in stages
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(probes, log2, c=color_map, alpha=0.4, s=8, linewidths=0)

    ax.axvline( args.min_probes, color="black", linestyle="--", linewidth=1.2,
                label=f"min probes = {args.min_probes}")
    ax.axhline( args.min_log2,  color="black", linestyle="-.", linewidth=1.2,
                label=f"|log2| = {args.min_log2}")
    ax.axhline(-args.min_log2,  color="black", linestyle="-.", linewidth=1.2)
    ax.axhline(0, color="grey", linewidth=0.7, linestyle=":")

    legend_handles = [
        mpatches.Patch(color=PALETTE["pass"],    label="Pass all filters"),
        mpatches.Patch(color=PALETTE["fail"],    label="Fail probes"),
        mpatches.Patch(color=PALETTE["neutral"], label="Fail log2"),
        mpatches.Patch(color=PALETTE["orange"],  label="Fail weight / COSMIC"),
    ]
    ax.legend(handles=legend_handles, fontsize=8)
    ax.set_xlabel("Probes per segment")
    ax.set_ylabel("log2 ratio")
    ax.set_title("Probes vs Log2 Ratio — joint threshold view", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# =============================================================================
# Figure 6 -- COSMIC gene coverage
# =============================================================================

def fig_cosmic(metrics, args):
    gene_counts = metrics["gene_counts"]
    gene_tier   = metrics["gene_tier"]

    if not gene_counts:
        return None

    n_tier1 = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if gene_tier.get(g) == 1)
    n_tier2 = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if gene_tier.get(g) == 2)
    n_non   = sum(v["gain"] + v["loss"] for g, v in gene_counts.items() if g not in gene_tier)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: donut
    sizes        = [n_tier1, n_tier2, n_non]
    clrs         = [PALETTE["tier1"], PALETTE["tier2"], PALETTE["neutral"]]
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
    ax1.set_title("COSMIC Gene Coverage\n(weight-passing segments)", fontweight="bold")

    # Panel B: top 20 genes, stacked GAIN / LOSS bars
    top_n       = 20
    total_hits  = {g: v["gain"] + v["loss"] for g, v in gene_counts.items()}
    top_genes   = sorted(total_hits.items(), key=lambda x: x[1], reverse=True)[:top_n]
    gnames      = [g for g, _ in top_genes]
    gains       = [gene_counts[g]["gain"] for g in gnames]
    losses      = [gene_counts[g]["loss"] for g in gnames]

    y = range(len(gnames))
    ax2.barh(list(y), gains,  color=PALETTE["gain"], edgecolor="white", height=0.65, label="Gain")
    ax2.barh(list(y), losses, left=gains, color=PALETTE["loss"],
             edgecolor="white", height=0.65, label="Loss")

    for i, (g_count, l_count, gname) in enumerate(zip(gains, losses, gnames)):
        total = g_count + l_count
        tier  = gene_tier.get(gname)
        marker = " ★" if tier == 1 else " ◆" if tier == 2 else ""
        ax2.text(total, i, f"  {total:,}{marker}", va="center", ha="left", fontsize=8)

    ax2.set_yticks(list(y))
    ax2.set_yticklabels(gnames, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("Segment count")
    ax2.set_title(f"Top {top_n} Genes — Gain/Loss breakdown\n(★ Tier 1  ◆ Tier 2)",
                  fontweight="bold")
    ax2.set_xlim(0, max(g + l for g, l in zip(gains, losses)) * 1.25 if gains else 1)
    ax2.legend(fontsize=8)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    return fig


# =============================================================================
# Figure 7 -- Genome-wide CNV profile
# =============================================================================

def fig_genome(metrics, args):
    segs = metrics["segs"].copy()
    if segs.empty:
        return None

    chrom_order = {c: i for i, c in enumerate(CHR_ORDER)}
    segs["chrom_clean"] = segs["chrom"].str.replace("chr", "", regex=False)
    segs["chrom_idx"]   = segs["chrom_clean"].map(chrom_order)
    segs = segs.dropna(subset=["chrom_idx"])
    segs["chrom_idx"] = segs["chrom_idx"].astype(int)

    # Scale midpoint within each chromosome to [0, 1] then shift by chrom_idx
    chrom_max = segs.groupby("chrom_idx")["midpt"].transform("max").clip(lower=1)
    segs["x"] = segs["chrom_idx"] + segs["midpt"] / chrom_max

    fig, ax = plt.subplots(figsize=(14, 5))

    # Excluded segments (grey, small, behind)
    excl = segs[segs["stage"] != "pass"]
    if not excl.empty:
        ax.scatter(excl["x"], excl["log2"],
                   c=PALETTE["neutral"], alpha=0.2, s=5, linewidths=0)

    # Passing segments coloured by GAIN / LOSS
    passed = segs[segs["stage"] == "pass"]
    if not passed.empty:
        gains  = passed[passed["log2"] > 0]
        losses = passed[passed["log2"] <= 0]
        if not gains.empty:
            ax.scatter(gains["x"],  gains["log2"],
                       c=PALETTE["gain"], alpha=0.7, s=12, linewidths=0, label="Gain")
        if not losses.empty:
            ax.scatter(losses["x"], losses["log2"],
                       c=PALETTE["loss"], alpha=0.7, s=12, linewidths=0, label="Loss")

        # COSMIC-hitting segments overlaid as open circles
        cosmic_hits = passed[passed["cosmic_genes"].apply(len) > 0]
        if not cosmic_hits.empty:
            ax.scatter(cosmic_hits["x"], cosmic_hits["log2"],
                       facecolors="none", edgecolors="black",
                       s=45, linewidths=0.8, zorder=5, label="COSMIC gene")

    # Threshold lines
    ax.axhline( args.min_log2,  color="black", linestyle="--", linewidth=0.8, alpha=0.45)
    ax.axhline(-args.min_log2,  color="black", linestyle="--", linewidth=0.8, alpha=0.45)
    ax.axhline(0, color="grey", linewidth=0.6, linestyle=":")

    # Chromosome boundary lines and x-axis labels
    chrs_present = sorted(segs["chrom_idx"].unique())
    for idx in chrs_present:
        ax.axvline(idx, color="lightgrey", linewidth=0.5, zorder=0)

    tick_pos    = [i + 0.5 for i in chrs_present]
    tick_labels = [CHR_ORDER[i] for i in chrs_present]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=7)

    # Add excluded to legend manually
    handles, lbls = ax.get_legend_handles_labels()
    handles.append(mpatches.Patch(color=PALETTE["neutral"], alpha=0.5, label="Excluded"))
    ax.legend(handles=handles, fontsize=8, loc="upper right")

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("log2 ratio")
    ax.set_title("Genome-wide CNV Profile", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
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
        rows.append({"filter_step": step, "segments": count, "pct_removed_from_prev": pct})
        prev = count
    tsv_path = os.path.join(report_dir, f"{sample_name}.cnv_filter_qc_stats.tsv")
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    print("\n=== CUP WES -- CNV Segment Filter QC Report ===\n")
    print(f"  Input CNS:    {args.cns}")
    print(f"  COSMIC list:  {args.cosmic}")
    print(f"  Min probes:   {args.min_probes}")
    print(f"  Min |log2|:   {args.min_log2}")
    print(f"  Min weight:   {args.min_weight}")
    print(f"  COSMIC tier:  {args.tier if args.tier else 'All'}\n")

    sample_name = os.path.basename(args.cns).replace(".tumor_only.call.cns", "")
    tier_suffix = f".tier{args.tier}" if args.tier else ""
    report_dir  = os.path.join(args.outdir, "reports", "dependantfilter")
    os.makedirs(report_dir, exist_ok=True)

    print("Loading COSMIC Cancer Gene Census...")
    gene_tier = load_cosmic_genes(args.cosmic, args.tier)

    print(f"\nLoading CNS file...")
    metrics = collect_metrics(args, gene_tier)
    print(f"  Done — {metrics['funnel']['Total segments']:,} segments loaded.\n")

    print("=== Filtering Summary ===\n")
    for step, count in metrics["funnel"].items():
        print(f"  {step:<45} {count:>8,}")
    print()

    figure_defs = [
        ("01_funnel",               fig_funnel(metrics, args)),
        ("02_probe_distribution",   fig_probes(metrics, args)),
        ("03_log2_distribution",    fig_log2(metrics, args)),
        ("04_weight_distribution",  fig_weight(metrics, args)),
        ("05_probes_vs_log2",       fig_probes_log2_scatter(metrics, args)),
        ("06_cosmic_genes",         fig_cosmic(metrics, args)),
        ("07_genome_profile",       fig_genome(metrics, args)),
    ]
    figures = [(name, fig) for name, fig in figure_defs if fig is not None]

    print("Saving figures...")
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

    print("Done.\n")


if __name__ == "__main__":
    main()
