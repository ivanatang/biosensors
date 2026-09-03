"""Aggregates per-sequence salt-bridge occupancy tables into one feature CSV.

Aggregates per-sequence salt-bridge occupancy tables (output of
salt_bridge_analysis.py) into a single wide-format CSV, one row per
sequence, ready to merge into the ML feature table on seq_id. Mirrors the
seq_ids.txt-driven path-construction / config.yaml conventions used by
water_analysis/aggregate_r_scores.py.

Per-sequence output columns:
    max_saltbridge_occupancy_pct: occupancy of the single most-contacted
        basic residue.
    n_saltbridges_gt50pct: count of residues with occupancy >= 50%.
    mean_top3_occupancy_pct: mean occupancy of the top-3 most-contacted
        residues (missing slots treated as 0% occupancy, so this stays
        comparable across sequences with fewer than 3 detected salt
        bridges).

Usage:
    python aggregate_salt_bridge_feats.py --config config.yaml --seq_list seq_ids.txt \
        --out_csv salt_bridge/saltbridge_features_all_seqs.csv
"""
import os
import argparse
import numpy as np
import pandas as pd
import yaml


def parse_args():
    """Parses CLI arguments for the salt-bridge feature aggregation run.

    Returns:
        argparse.Namespace: Parsed arguments (config, seq_list, out_csv,
        start_ns, end_ns).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default="config.yaml",
                        help='Path to config.yaml (default: config.yaml in repo root)')
    parser.add_argument('--seq_list', default="seq_ids.txt",
                        help='Text file with seq_id, seq_type, and optional '
                             'custom base path, one per line')
    parser.add_argument('--out_csv', required=True)
    parser.add_argument('--start-ns', type=float, default=40.0,
                        help='Start of analysis window in ns (default: 40), '
                             'must match the value used for salt_bridge_analysis.py')
    parser.add_argument('--end-ns', type=float, default=500.0,
                        help='End of analysis window in ns (default: 500), '
                             'must match the value used for salt_bridge_analysis.py')
    return parser.parse_args()


def load_seq_list(seq_list_path):
    """Loads (seq_id, seq_type, custom_path) entries from a seq_ids file.

    Args:
        seq_list_path (str): Path to a seq_ids.txt-style file. Accepts
            tab- or comma-separated columns. Blank lines and lines
            starting with "#" are skipped.

    Returns:
        list[tuple[str, str, str | None]]: (seq_id, seq_type, custom_path)
        per entry; seq_type defaults to "Unknown" and custom_path to None
        when absent.
    """
    entries = []
    with open(seq_list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.replace(',', '\t').split('\t')]
            seq_id   = parts[0]
            seq_type = parts[1] if len(parts) > 1 else 'Unknown'
            custom   = parts[2] if len(parts) > 2 else None
            entries.append((seq_id, seq_type, custom))
    return entries


def get_occupancy_path(seq_id, seq_type, custom_dir, cfg, tag):
    """Resolves the path to a sequence's salt-bridge occupancy CSV.

    Mirrors salt_bridge_analysis.py: the default window (40-500ns) reads
    from the original, untagged directory; any other window reads from its
    own tagged subdirectory.

    Args:
        seq_id (str): Sequence identifier.
        seq_type (str): seq_type label, used only if `seq_id` doesn't carry
            a recognized suffix.
        custom_dir (str | None): Custom base path from seq_ids.txt, if any.
        cfg (dict): Parsed config.yaml.
        tag (str): Window tag, e.g. "40_500ns".

    Returns:
        str: Path to the occupancy CSV for this sequence/window.
    """
    sb_cfg = cfg["salt_bridge"]
    occ_file = sb_cfg["output_files"]["occupancy"]
    subdir = sb_cfg["output_subdir"] if tag == "40_500ns" else f"{sb_cfg['output_subdir']}_{tag}"

    if custom_dir:
        return os.path.join(custom_dir, subdir, occ_file)

    base   = os.path.expandvars(cfg["paths"]["base"])
    runrel = cfg["paths"]["runrel"]

    if seq_id.endswith('_binder'):
        type_subdir = cfg["paths"]["type_subdir"]["binder"]
    elif seq_id.endswith('_nb'):
        type_subdir = cfg["paths"]["type_subdir"]["nb"]
    elif seq_id.endswith('_low_pkt'):
        type_subdir = cfg["paths"]["type_subdir"]["low_pkt"]
    elif seq_id.endswith('_fail_gate'):
        type_subdir = cfg["paths"]["type_subdir"]["fail_gate"]
    else:
        type_subdir = cfg["paths"]["type_subdir"].get(seq_type, "")

    return os.path.join(base, type_subdir, seq_id, runrel, subdir, occ_file)


def summarize(occ_df):
    """Summarizes a sequence's per-residue salt-bridge occupancy table.

    Args:
        occ_df (pandas.DataFrame): Occupancy table with an
            "occupancy_pct" column (one row per residue).

    Returns:
        tuple[float, int, float]: (max occupancy %, count of residues with
        occupancy >= 50%, mean occupancy % of the top 3 residues, missing
        slots treated as 0%). (0.0, 0, 0.0) if `occ_df` is empty.
    """
    if occ_df.empty:
        return 0.0, 0, 0.0

    occ = occ_df["occupancy_pct"].to_numpy()
    top3 = np.sort(occ)[::-1][:3]
    top3_padded = np.pad(top3, (0, max(0, 3 - len(top3))), constant_values=0.0)

    max_occ    = occ.max()
    n_gt50     = int((occ >= 50).sum())
    mean_top3  = float(top3_padded.mean())
    return float(max_occ), n_gt50, mean_top3


if __name__ == "__main__":
    args = parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    seq_list = load_seq_list(args.seq_list)
    tag = f"{int(args.start_ns)}_{int(args.end_ns)}ns"

    rows    = []
    missing = []

    for seq_id, seq_type, custom_dir in seq_list:
        path = get_occupancy_path(seq_id, seq_type, custom_dir, cfg, tag)

        if not os.path.exists(path):
            print(f"  MISSING: {path}")
            missing.append(seq_id)
            continue

        occ_df = pd.read_csv(path)
        max_occ, n_gt50, mean_top3 = summarize(occ_df)

        rows.append({
            "seq_id":                        seq_id,
            "seq_type":                      seq_type,
            "max_saltbridge_occupancy_pct":  round(max_occ, 2),
            "n_saltbridges_gt50pct":         n_gt50,
            "mean_top3_occupancy_pct":       round(mean_top3, 2),
        })
        print(f"  Loaded {seq_id}  [{seq_type}]  "
              f"max={max_occ:.1f}%  n>=50%={n_gt50}  top3_mean={mean_top3:.1f}%")

    if not rows:
        raise RuntimeError("No sequences loaded. Check --config/--seq_list paths.")

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    out_df.to_csv(args.out_csv, index=False)
    print(f"\nSalt-bridge features saved -> {args.out_csv}")
    print(f"  Shape: {out_df.shape}  ({len(out_df)} sequences x 3 feature columns)")

    if missing:
        print(f"\nMissing ({len(missing)}): {missing}")
