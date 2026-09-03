"""
plot_Rg_sasa.py — Plot per-group Rg and SASA distributions from GROMACS xvg output.

Usage:
    python plot_Rg_sasa.py [--region pocket|whole] [seq_ids.txt] [--out-dir DIR]
                            [--structure-source {ngs_observed,designed_assumed,all}]
                            [--structure-guide md_candidate_guide.csv]

Reads Rg_{suffix}.xvg and sasa_{suffix}.xvg from each sequence's run directory,
computes the per-trajectory mean, and saves two comparison plots.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Constants ─────────────────────────────────────────────────────────────────
BASE = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

TYPE_SUBDIR = {
    "Binder": "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry": "neg_fail_gate",
}

GROUP_COLOR = {
    "Binder": "#648FFF",
    "False Positive": "#DC267F",
    "Low Confidence": "#FE6100",
    "Fail Geometry": "#FFB000",
}

GROUP_ORDER = ["Binder", "False Positive", "Low Confidence", "Fail Geometry"]


# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_xvg(path):
    """Reads a GROMACS .xvg file's second column, skipping header lines.

    Args:
        path (str): Path to the .xvg file.

    Returns:
        numpy.ndarray: Values from column 2 of each data line.
    """
    vals = []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                vals.append(float(parts[1]))
    return np.array(vals)


def read_seq_ids(path):
    """Reads a tab-separated seq_ids.txt-style file.

    Args:
        path (str): Path to the file (columns: folder_name, label,
            optional custom_base). Blank lines and lines starting with
            "#" are skipped.

    Returns:
        list[tuple]: (folder_name, label, custom_base) per row; custom_base
        is "" when absent.
    """
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            folder_name = parts[0]
            label = parts[1] if len(parts) > 1 else ""
            custom_base = parts[2] if len(parts) > 2 else ""
            seqs.append((folder_name, label, custom_base))
    return seqs


# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention (pair_XXXX_binder/_nb/_low_pkt/_fail_gate) used
# elsewhere in the repo.
MD_GROUP_SUFFIX = {
    "binder": "binder",
    "non_binder": "nb",
    "negative_low_pocket": "low_pkt",
    "negative_fail_gate": "fail_gate",
}


def load_source_ids(guide_path, source):
    """Loads seq_ids from md_candidate_guide.csv matching one source value.

    Args:
        guide_path (str): Path to md_candidate_guide.csv.
        source (str): Value to match in the `source` column (e.g.
            "ngs_observed" for sequencing-confirmed via Y2H/FACS sort-seq,
            vs. "designed_assumed").

    Returns:
        set: seq_ids (pair_id + mapped md_group suffix) matching `source`.
    """
    guide = pd.read_csv(guide_path)
    matched = guide[guide["source"] == source].copy()
    matched["seq_id"] = matched.apply(
        lambda r: f"{r['pair_id']}_{MD_GROUP_SUFFIX.get(r['md_group'], r['md_group'])}",
        axis=1)
    return set(matched["seq_id"])


def get_rundir(folder_name, label, custom_base):
    """Resolves a sequence's run directory.

    Args:
        folder_name (str): seq_id / folder name.
        label (str): Group label (key into TYPE_SUBDIR), used when
            `custom_base` isn't given.
        custom_base (str): Custom base path override; used as-is if
            non-empty.

    Returns:
        str: Path to the sequence's production run directory.
    """
    if custom_base:
        return os.path.join(custom_base, RUNREL)
    return os.path.join(BASE, TYPE_SUBDIR[label], folder_name, RUNREL)


def collect_means(seqs, xvg_name):
    """Collects each sequence's per-trajectory mean from one .xvg file.

    Args:
        seqs (list[tuple]): (folder_name, label, custom_base) rows (see
            read_seq_ids).
        xvg_name (str): .xvg filename to read from each run directory.

    Returns:
        dict: Maps group label to a list of per-trajectory mean values.
    """
    data = {g: [] for g in GROUP_ORDER}
    for folder_name, label, custom_base in seqs:
        if label not in TYPE_SUBDIR:
            continue
        rundir = get_rundir(folder_name, label, custom_base)
        xvg_path = os.path.join(rundir, xvg_name)
        if not os.path.exists(xvg_path):
            print(f"  WARNING: {xvg_path} not found — skipping")
            continue
        vals = parse_xvg(xvg_path)
        if vals.size > 0:
            data[label].append(vals.mean())
    return data


def make_plot(data, ylabel, title, out_path):
    """Draws a box-and-jitter plot of per-group values and saves it.

    Args:
        data (dict): Maps group label to a list of values (see
            collect_means).
        ylabel (str): Y-axis label.
        title (str): Plot title.
        out_path (str): Path to save the PNG to.
    """
    groups = [g for g in GROUP_ORDER if data.get(g)]
    if not groups:
        print(f"  No data to plot for {title}")
        return

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300, constrained_layout=True)
    rng = np.random.default_rng(42)

    for i, group in enumerate(groups):
        vals = np.array(data[group])
        color = GROUP_COLOR[group]

        ax.boxplot(
            vals,
            positions=[i],
            widths=0.45,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            boxprops=dict(facecolor=color, alpha=0.5),
            whiskerprops=dict(color=color),
            capprops=dict(color=color),
            flierprops=dict(marker="", linestyle="none"),
        )

        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(i + jitter, vals, color=color, alpha=0.85, s=18, zorder=3)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=15, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.4)

    legend_patches = [
        mpatches.Patch(facecolor=GROUP_COLOR[g], label=f"{g} (n={len(data[g])})")
        for g in groups
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=8)

    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    """Loads Rg/SASA data for all sequences and saves the group comparison plots."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--region", choices=["pocket", "whole"], default="whole",
        help="Region to plot: 'pocket' or 'whole' (default: whole)",
    )
    parser.add_argument(
        "seq_ids", nargs="?", default="seq_ids.txt",
        help="Path to seq_ids.txt (default: seq_ids.txt)",
    )
    parser.add_argument(
        "--out-dir", default=".", help="Directory for output PNGs (default: .)"
    )
    parser.add_argument(
        "--structure-source", default="all",
        choices=["ngs_observed", "designed_assumed", "all"],
        help="Filter sequences by md_candidate_guide.csv's source column. "
             "'all' (default) applies no filtering.",
    )
    parser.add_argument(
        "--structure-guide",
        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
        help="Path to md_candidate_guide.csv (default: %(default)s)",
    )
    args = parser.parse_args()

    suffix = "pocket" if args.region == "pocket" else "PL"
    region_label = "Pocket Residues" if args.region == "pocket" else "Protein + Ligand"

    seqs = read_seq_ids(args.seq_ids)
    print(f"Loaded {len(seqs)} sequences from {args.seq_ids}  |  region: {args.region}")

    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        before = len(seqs)
        seqs = [s for s in seqs if s[0] in source_ids]
        print(f"Restricting to source == {args.structure_source}: "
              f"{before} -> {len(seqs)}")

    rg_data = collect_means(seqs, f"Rg_{suffix}.xvg")
    sasa_data = collect_means(seqs, f"sasa_{suffix}.xvg")

    make_plot(
        rg_data,
        ylabel="Mean Rg (nm)",
        title=f"Radius of Gyration by Group — {region_label}",
        out_path=os.path.join(args.out_dir, f"Rg_{suffix}_comparison.png"),
    )

    make_plot(
        sasa_data,
        ylabel="Mean SASA (nm²)",
        title=f"Solvent-Accessible Surface Area by Group — {region_label}",
        out_path=os.path.join(args.out_dir, f"sasa_{suffix}_comparison.png"),
    )


if __name__ == "__main__":
    main()
