"""Aggregates per-frame hydration counts into per-sequence summary features.

Reads each sequence's {seq_id}_hydration_{region}_{TAG}.csv (produced by
hydration_calc.py, one row per frame: time_ns plus one hydration_count_*
column per cutoff radius) and computes summary features for every
hydration_count_* column present, mirroring
pkt_vol/extract_pkt_feats.py's early/late/drift/slope column family for
consistency with the rest of the feature table.

Early/late windows are fractional (first/last 20% of frames, matching
extract_pkt_feats.py) rather than a fixed-ns window, since sequences can
vary in how many frames they cover after striding. Because
hydration_calc.py records real simulation time (unlike mdpocket's
descriptors.txt), both a fractional-position slope and an absolute-ns
slope are emitted.

Output: water_density_feats.csv

Usage:
    python aggregate_hydration_feats.py [seq_ids_ngs_observed.txt] [--start-ns 40] [--end-ns 500]
                                         [--reference-region {ligand,pocket_residues}]
                                         [--early-late-frac 0.2] [--out water_density_feats.csv]
                                         [--structure-source {ngs_observed,designed_assumed,all}]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import linregress

# ── Configurable paths ────────────────────────────────────────────────────────
BASE = "/scratch/alpine/ivta1597/LCA_boltz_models"
# ─────────────────────────────────────────────────────────────────────────────

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

# md_candidate_guide.csv's md_group values use different suffixes than the
# seq_id naming convention used everywhere else in the repo.
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


def extract_features(df, window_frac):
    """Computes summary statistics for every hydration_count_* column.

    One column per cutoff radius written by hydration_calc.py. Includes
    an early-vs-late trajectory comparison over a fractional window (since
    frame counts can vary across sequences after striding).

    Args:
        df: Per-frame hydration DataFrame (time_ns plus hydration_count_*
            columns; see hydration_calc.py).
        window_frac (float): Fraction of frames counted as "early"/"late"
            for the drift comparison (e.g. 0.2 for the first/last 20%).

    Returns:
        tuple[dict, list[str]]: (features, count_cols) -- summary features
        (mean/std/min/max, early/late means, drift, fractional and
        per-ns slope, per hydration_count_* column, plus n_frames) and
        the list of hydration_count_* column names found.
    """
    t = df["time_ns"].values
    n = len(t)
    pct = int(round(window_frac * 100))
    k = max(1, int(round(window_frac * n)))
    pos = np.arange(n) / (n - 1) if n > 1 else np.zeros(n)

    count_cols = [c for c in df.columns if c.startswith("hydration_count_")]
    features = {}
    for col in count_cols:
        cnt = df[col].values
        early_mean = float(cnt[:k].mean())
        late_mean  = float(cnt[n - k:].mean())
        frac_slope = float(linregress(pos, cnt).slope) if n > 1 else 0.0
        ns_slope   = float(linregress(t, cnt).slope) if n > 1 else 0.0

        features.update({
            f"{col}_mean":        float(np.mean(cnt)),
            f"{col}_std":         float(np.std(cnt)),
            f"{col}_min":         float(np.min(cnt)),
            f"{col}_max":         float(np.max(cnt)),
            f"{col}_early{pct}_mean": early_mean,
            f"{col}_late{pct}_mean":  late_mean,
            f"{col}_drift{pct}":      late_mean - early_mean,
            f"{col}_slope":       frac_slope,
            f"{col}_slope_per_ns": ns_slope,
        })

    features["n_frames"] = n
    return features, count_cols


def main():
    """Aggregates hydration CSVs across a seq list into one feature table."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?",
                        default="/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt")
    parser.add_argument("--start-ns", type=float, default=40.0)
    parser.add_argument("--end-ns",   type=float, default=500.0)
    parser.add_argument("--reference-region", choices=["ligand", "pocket_residues"],
                        default="ligand",
                        help="Which hydration_calc.py reference-region output to read "
                             "(default: ligand)")
    parser.add_argument("--early-late-frac", type=float, default=0.2,
                        help="Fraction of each sequence's own frames counted as "
                             "'early'/'late' for the drift comparison (default: 0.2).")
    parser.add_argument("--output", default=None,
                        help="Output CSV filename (default: water_density_feats.csv, "
                             "or water_density_feats{suffix}.csv if --suffix is set)")
    parser.add_argument("--base", default=BASE)
    parser.add_argument("--suffix", default="",
                        help="Run-directory suffix, e.g. '_qfix', matching whatever "
                             "hydration_calc.py was run with for these sequences "
                             "(default: '', the standard directory). Also folded into "
                             "the default --output filename so a _qfix run never "
                             "overwrites the standard one.")
    parser.add_argument("--structure-source", default="all",
                        choices=["ngs_observed", "designed_assumed", "all"],
                        help="Filter sequences by md_candidate_guide.csv's source column. "
                             "'all' (default) applies no filtering.")
    parser.add_argument("--structure-guide",
                        default="/projects/ivta1597/biosensors/md_candidate_guide.csv",
                        help="Path to md_candidate_guide.csv (default: %(default)s)")
    args = parser.parse_args()
    if args.output is None:
        args.output = f"water_density_feats{args.suffix}.csv"

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    source_ids = None
    if args.structure_source != "all":
        source_ids = load_source_ids(args.structure_guide, args.structure_source)
        print(f"Restricting to source == {args.structure_source}: "
              f"{len(source_ids)} in {args.structure_guide}")

    TAG = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
    REGION_TAG = {"ligand": "ligand", "pocket_residues": "pocket"}[args.reference_region]

    records = []
    missing = []
    skipped_wrong_source = []
    all_count_cols = []   # preserves discovery order across sequences

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
                run_dir = os.path.join(custom_path, f"water_density_{TAG}{args.suffix}")
            else:
                dir_type = TYPE_SUBDIR.get(seq_type, seq_type)
                run_dir  = os.path.join(args.base, dir_type, seq_id, f"water_density_{TAG}{args.suffix}")

            csv_file = os.path.join(run_dir, f"{seq_id}_hydration_{REGION_TAG}_{TAG}.csv")

            if not os.path.exists(csv_file):
                print(f"MISSING: {seq_id}  [{seq_type}]  ->  {csv_file}")
                missing.append(seq_id)
                continue

            try:
                df = pd.read_csv(csv_file)
                features, count_cols = extract_features(df, args.early_late_frac)
                features["seq_id"]   = seq_id
                features["seq_type"] = seq_type
                records.append(features)
                for c in count_cols:
                    if c not in all_count_cols:
                        all_count_cols.append(c)
                headline = count_cols[0] if count_cols else None
                headline_mean = features.get(f"{headline}_mean") if headline else None
                if headline_mean is not None:
                    print(f"OK: {seq_id}  [{seq_type}]  "
                          f"{headline}_mean={headline_mean:.2f}  n={features['n_frames']}")
                else:
                    print(f"OK: {seq_id}  [{seq_type}]  n={features['n_frames']}")
            except Exception as e:
                print(f"ERROR: {seq_id}  -  {e}")
                missing.append(seq_id)

    if not records:
        print("\nNo sequences loaded - nothing to write.")
        sys.exit(1)

    feat_df = pd.DataFrame(records)
    pct = int(round(args.early_late_frac * 100))
    col_order = ["seq_id", "seq_type"]
    for col in all_count_cols:
        col_order += [f"{col}_mean", f"{col}_std", f"{col}_min", f"{col}_max",
                      f"{col}_early{pct}_mean", f"{col}_late{pct}_mean",
                      f"{col}_drift{pct}", f"{col}_slope", f"{col}_slope_per_ns"]
    col_order.append("n_frames")
    feat_df = feat_df[col_order]
    feat_df.to_csv(args.output, index=False)

    print(f"\nFeatures written to: {args.output}")
    print(f"  Sequences processed : {len(records)}")
    print(f"  Sequences missing   : {len(missing)}")
    if missing:
        print(f"  Missing seq_ids     : {', '.join(missing)}")
    if source_ids is not None:
        print(f"  Skipped (source != {args.structure_source}) : {len(skipped_wrong_source)}")


if __name__ == "__main__":
    main()
