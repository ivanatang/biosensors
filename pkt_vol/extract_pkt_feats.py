"""
extract_pocket_features.py

Extracts per-sequence pocket volume summary features from mdpocket
characterization descriptor files for all sequences in seq_ids.txt.

Features extracted per sequence:
    pocket_vol_mean         : mean pocket volume across all frames (Å³)
    pocket_vol_std          : std dev of pocket volume
    pocket_vol_min          : minimum pocket volume
    pocket_vol_max          : maximum pocket volume
    pocket_vol_closed_frac  : fraction of frames below closure threshold

Output: pocket_volume_features.csv

Usage:
    python extract_pocket_features.py [seq_ids.txt] [--threshold 800] [--plot]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── Configurable paths ────────────────────────────────────────────────────────
BASE   = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
# ─────────────────────────────────────────────────────────────────────────────

GROUP_COLORS = {
    "Binder":         "#70AD47",
    "False Positive": "#d62728",
    "Low Confidence": "#e8756a",
    "Fail Geometry":  "#f5b7b1",
}

def get_dir_type(seq_type):
    mapping = {
        "Binder":         "binders",
        "False Positive": "nonbinders",
        "Low Confidence": "neg_low_pkt",
        "Fail Geometry":  "neg_fail_gate",
    }
    return mapping.get(seq_type, seq_type)


def load_descriptors(desc_path):
    """Load mdpocket descriptors file and return DataFrame."""
    df = pd.read_csv(desc_path, sep=r'\s+')
    return df


def extract_features(df, threshold):
    """Compute summary statistics from per-frame pocket volume column."""
    vol = df["pock_volume"].values
    return {
        "pocket_vol_mean":        np.mean(vol),
        "pocket_vol_std":         np.std(vol),
        "pocket_vol_min":         np.min(vol),
        "pocket_vol_max":         np.max(vol),
        "pocket_vol_closed_frac": np.mean(vol < threshold),
        "n_frames":               len(vol),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list",   nargs="?", default="seq_ids.txt")
    parser.add_argument("--threshold", type=float, default=800.0,
                        help="Volume threshold for closed_frac calculation in Å³ (default: 800)")
    parser.add_argument("--plot",      action="store_true",
                        help="Generate summary plots")
    parser.add_argument("--output",    default="pocket_volume_features.csv",
                        help="Output CSV filename (default: pocket_volume_features.csv)")
    args = parser.parse_args()

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    records  = []
    missing  = []

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts       = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            if custom_path:
                run_dir = os.path.join(custom_path, RUNREL)
            else:
                run_dir = os.path.join(BASE, get_dir_type(seq_type), seq_id, RUNREL)

            desc_file = os.path.join(run_dir, f"mdpocket_{seq_id}_descriptors.txt")

            if not os.path.exists(desc_file):
                print(f"MISSING: {seq_id}  [{seq_type}]  →  {desc_file}")
                missing.append(seq_id)
                continue

            try:
                df       = load_descriptors(desc_file)
                features = extract_features(df, args.threshold)
                features["seq_id"]   = seq_id
                features["seq_type"] = seq_type
                records.append(features)
                print(f"OK: {seq_id}  [{seq_type}]  "
                      f"mean={features['pocket_vol_mean']:.1f} Å³  "
                      f"std={features['pocket_vol_std']:.1f}  "
                      f"n={features['n_frames']}")
            except Exception as e:
                print(f"ERROR: {seq_id}  —  {e}")
                missing.append(seq_id)

    if not records:
        print("\nNo descriptors loaded — nothing to write.")
        sys.exit(1)

    # ── Build feature dataframe ───────────────────────────────────────────────
    feat_df = pd.DataFrame(records)
    col_order = ["seq_id", "seq_type", "pocket_vol_mean", "pocket_vol_std",
                 "pocket_vol_min", "pocket_vol_max", "pocket_vol_closed_frac",
                 "n_frames"]
    feat_df = feat_df[col_order]
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    if args.plot:
        # 1. Scatter by group
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        for group, color in GROUP_COLORS.items():
            subset = feat_df[feat_df["seq_type"] == group]
            ax.scatter([group] * len(subset), subset["pocket_vol_mean"],
                       color=color, s=60, zorder=3, alpha=0.8, label=group)
        ax.set_ylabel("Mean pocket volume (Å³)")
        ax.set_title("Binding site pocket volume by group")
        ax.grid(True, alpha=0.4)
        plt.savefig("pocket_vol_by_group.png", dpi=150)
        print("\nSaved: pocket_vol_by_group.png")

        # 2. Boxplot
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        groups  = [g for g in GROUP_COLORS if g in feat_df["seq_type"].values]
        data    = [feat_df[feat_df["seq_type"] == g]["pocket_vol_mean"].values
                   for g in groups]
        colors  = [GROUP_COLORS[g] for g in groups]
        bp = ax.boxplot(data, patch_artist=True, labels=groups)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel("Mean pocket volume (Å³)")
        ax.set_title("Pocket volume distribution by group")
        ax.grid(True, alpha=0.4)
        plt.savefig("pocket_vol_boxplot.png", dpi=150)
        print("Saved: pocket_vol_boxplot.png")

        plt.close("all")


if __name__ == "__main__":
    main()
