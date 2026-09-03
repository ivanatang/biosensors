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
    pocket_vol_early{P}_mean / late{P}_mean / drift{P} / early{P}_closed_frac /
    late{P}_closed_frac : early vs. late trajectory comparison, mirroring
        extract_gate_latch_rmsd_feats.py's early/late RMSD windows. Windows are
        defined as a FRACTION of each sequence's own frame count (default: first/
        last 20%), not an absolute ns window, because many sequences' descriptor
        files cover far fewer frames than the nominal 500ns run (see n_frames) --
        an absolute-ns window like the RMSD script's would silently exclude most
        short trajectories. This makes drift comparable in relative trajectory
        position across sequences, not in absolute simulated time.
    pocket_vol_slope        : linear trend of pocket volume vs. fractional
        trajectory position (0=first frame, 1=last frame), in Å³ per full
        trajectory traversed.

Output: pocket_volume_features.csv

Usage:
    python extract_pocket_features.py [seq_ids_orig.txt] [--threshold 800] [--plot]
                                       [--early-late-frac 0.2]
                                       [--structure-source {ngs_observed,designed_assumed,all}]
                                       [--structure-guide md_candidate_guide.csv]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

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
    """Maps a seq_type label to its data subdirectory.

    Args:
        seq_type (str): One of "Binder", "False Positive", "Low
            Confidence", "Fail Geometry".

    Returns:
        str: Subdirectory name, or `seq_type` unchanged if unrecognized.
    """
    mapping = {
        "Binder":         "binders",
        "False Positive": "nonbinders",
        "Low Confidence": "neg_low_pkt",
        "Fail Geometry":  "neg_fail_gate",
    }
    return mapping.get(seq_type, seq_type)


# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention (pair_XXXX_binder / _nb / _low_pkt / _fail_gate)
# used everywhere else in the repo.
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


def load_descriptors(desc_path):
    """Loads an mdpocket descriptors file.

    Args:
        desc_path (str): Path to a whitespace-delimited descriptors.txt.

    Returns:
        pandas.DataFrame: Parsed descriptor table.
    """
    df = pd.read_csv(desc_path, sep=r'\s+')
    return df


def extract_features(df, threshold, window_frac):
    """Computes summary pocket-volume features, including an early-vs-late
    trajectory drift comparison.

    Early/late windows are defined by fractional position within this
    sequence's own frame count (first/last `window_frac` of frames), not
    absolute time: descriptors.txt has no time column, and sequences vary
    widely in how many frames they actually cover (see n_frames), so a
    fixed-ns window (as used for gate/latch RMSD) would be undefined for
    short trajectories.

    Args:
        df (pandas.DataFrame): Descriptor table with a "pock_volume" column
            (see load_descriptors).
        threshold (float): Volume threshold (A^3) below which a frame
            counts as "closed".
        window_frac (float): Fraction of frames counted as "early"/"late"
            (e.g. 0.2 for first/last 20%).

    Returns:
        dict: Pocket-volume summary stats (mean/std/min/max, closed_frac,
        early/late means and closed_frac, drift, slope, n_frames).
    """
    vol = df["pock_volume"].values
    n   = len(vol)
    pct = int(round(window_frac * 100))

    # Index-sliced (not a position-threshold mask) so early/late windows are
    # never empty, even for very short trajectories where round(window_frac*n)
    # could otherwise land past the last frame satisfying a ">=" cutoff.
    k = max(1, int(round(window_frac * n)))
    early_vol = vol[:k]
    late_vol  = vol[n - k:]
    early_mean = float(early_vol.mean())
    late_mean  = float(late_vol.mean())
    pos        = np.arange(n) / (n - 1) if n > 1 else np.zeros(n)
    slope      = float(linregress(pos, vol).slope) if n > 1 else 0.0

    return {
        "pocket_vol_mean":        np.mean(vol),
        "pocket_vol_std":         np.std(vol),
        "pocket_vol_min":         np.min(vol),
        "pocket_vol_max":         np.max(vol),
        "pocket_vol_closed_frac": np.mean(vol < threshold),
        f"pocket_vol_early{pct}_mean":        early_mean,
        f"pocket_vol_late{pct}_mean":         late_mean,
        f"pocket_vol_drift{pct}":             late_mean - early_mean,
        f"pocket_vol_early{pct}_closed_frac": float(np.mean(early_vol < threshold)),
        f"pocket_vol_late{pct}_closed_frac":  float(np.mean(late_vol < threshold)),
        "pocket_vol_slope":       slope,
        "n_frames":               n,
    }


def main():
    """Extracts pocket-volume features for every sequence in a seq list."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list",   nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_orig.txt")
    parser.add_argument("--threshold", type=float, default=800.0,
                        help="Volume threshold for closed_frac calculation in Å³ (default: 800)")
    parser.add_argument("--plot",      action="store_true",
                        help="Generate summary plots")
    parser.add_argument("--early-late-frac", type=float, default=0.2,
                        help="Fraction of each sequence's own frames counted as "
                             "'early'/'late' for the drift comparison (default: 0.2, "
                             "i.e. first/last 20%%). Proportional rather than a fixed-ns "
                             "window since many sequences cover far fewer frames than "
                             "the nominal 500ns run.")
    parser.add_argument("--output",    default="pocket_volume_features.csv",
                        help="Output CSV filename (default: pocket_volume_features.csv)")
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    parser.add_argument("--suffix", default="",
                        help="Run-directory suffix, e.g. '_qfix' for the "
                             "bond-order/charge-fix systems (default: '', the "
                             "standard production directory). Applies to every "
                             "sequence in seq_list -- run _qfix sequences via a "
                             "separate seq_list from standard ones.")
    args = parser.parse_args()

    global RUNREL
    RUNREL = f"{RUNREL}{args.suffix}"

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    source_ids = None
    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        print(f"Restricting to source == {args.structure_source}: "
              f"{len(source_ids)} in {args.structure_guide}")

    records  = []
    missing  = []
    skipped_wrong_source = []

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts       = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            if source_ids is not None and seq_id not in source_ids:
                skipped_wrong_source.append(seq_id)
                continue

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
                features = extract_features(df, args.threshold, args.early_late_frac)
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
    pct = int(round(args.early_late_frac * 100))
    col_order = ["seq_id", "seq_type", "pocket_vol_mean", "pocket_vol_std",
                 "pocket_vol_min", "pocket_vol_max", "pocket_vol_closed_frac",
                 f"pocket_vol_early{pct}_mean", f"pocket_vol_late{pct}_mean",
                 f"pocket_vol_drift{pct}",
                 f"pocket_vol_early{pct}_closed_frac", f"pocket_vol_late{pct}_closed_frac",
                 "pocket_vol_slope", "n_frames"]
    feat_df = feat_df[col_order]
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")
    if source_ids is not None:
        print(f"  Skipped (source != {args.structure_source}) : {len(skipped_wrong_source)}")

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
