"""
Aggregate contact features into wide-format CSVs
=================================================
Reads {seq_id}_R_scores_{TAG}.csv for each sequence and produces two output CSVs:

  r_scores_all_sequences_{TAG}.csv
      Rows: one per sequence
      Columns: seq_id, seq_type, R_{resSeq} for all 181 residues
      NaN (no contact) written as empty cell.

  dw_scores_all_sequences_{TAG}.csv
      Rows: one per sequence
      Columns: seq_id, seq_type, D_{resSeq}, W_{resSeq} for all 181 residues
      Residues with no contact (I == 0) are written as 0.0 since zero
      occupancy is an unambiguous measurement, not a missing value.

Usage
-----
    python aggregate_r_scores.py --seq_list seq_ids.txt --out_dir /path/to/output
    python aggregate_r_scores.py --seq_list seq_ids.txt --out_dir /path/to/output --start-ns 40 --end-ns 250

seq_ids.txt -- tab or comma separated, one sequence per line.
Two-column format (default path construction):
    pair_3069_binder    Binder
    pair_3070_binder    Binder

Three-column format (custom base directory — parent of water_contacts_{TAG}):
    seq14_binder    Binder    /scratch/alpine/ivta1597/LCA_boltz_models/binders/seq14_binder/HMR/dodecahedron
    seq1_nb         False Positive    /scratch/alpine/ivta1597/LCA_boltz_models/nonbinders/seq1_nb/HMR/dodecahedron

NOTE: the third column should be the directory CONTAINING water_contacts_{TAG}/,
not the water_contacts directory itself. The script appends water_contacts_{TAG}
automatically so the same seq_ids.txt works for any time window.
"""

import os
import argparse
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
BASE = "/scratch/alpine/ivta1597/LCA_boltz_models"


# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq_list', default=None,
                        help='Text file with seq_id, seq_type, and optional '
                             'custom base path, one per line')
    parser.add_argument('--seq_ids',   nargs='+', default=None)
    parser.add_argument('--seq_types', nargs='+', default=None)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--base', default=BASE)
    parser.add_argument('--start-ns', type=float, default=40.0,
                        help='Start of analysis window in ns (default: 40)')
    parser.add_argument('--end-ns',   type=float, default=500.0,
                        help='End of analysis window in ns (default: 500)')
    return parser.parse_args()


def load_seq_list(args):
    """
    Returns a list of (seq_id, seq_type, custom_dir_or_None) tuples.
    custom_dir is the full path to the water_contacts directory, or None
    to use the default path construction.
    """
    if args.seq_list:
        entries = []
        with open(args.seq_list) as f:
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
    if args.seq_ids:
        types = args.seq_types or ['Unknown'] * len(args.seq_ids)
        return [(sid, stype, None) for sid, stype in zip(args.seq_ids, types)]
    raise ValueError("Provide either --seq_list or --seq_ids")


def get_csv_path(seq_id, base, custom_dir, tag):
    """
    Return the full path to the R_scores CSV for this sequence.

    For default paths, infers subdirectory from seq_id suffix and appends
    water_contacts_{tag}/{seq_id}_R_scores_{tag}.csv.

    For custom paths, custom_dir is the parent directory that CONTAINS
    water_contacts_{tag}/ (i.e. the sequence-level or HMR/dodecahedron dir).
    """
    if custom_dir:
        return os.path.join(custom_dir, f"water_contacts_{tag}",
                            f"{seq_id}_R_scores_{tag}.csv")

    if seq_id.endswith('_binder'):
        subdir = 'binders'
    elif seq_id.endswith('_nb'):
        subdir = 'nonbinders'
    elif seq_id.endswith('_low_pkt'):
        subdir = 'neg_low_pkt'
    elif seq_id.endswith('_fail_gate'):
        subdir = 'neg_fail_gate'
    else:
        subdir = None

    if subdir:
        return os.path.join(base, subdir, seq_id, f"water_contacts_{tag}",
                            f"{seq_id}_R_scores_{tag}.csv")

    return os.path.join(base, seq_id, f"water_contacts_{tag}",
                        f"{seq_id}_R_scores_{tag}.csv")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args     = parse_args()
    seq_list = load_seq_list(args)
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Derive time window tag ─────────────────────────────────────────────
    TAG = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
    print(f"Analysis window: {args.start_ns:.0f}–{args.end_ns:.0f} ns  "
          f"(tag: {TAG})")

    r_rows       = []   # R-score rows (NaN for no contact)
    dw_rows      = []   # D and W rows (0.0 for no contact)
    all_residues = None
    missing      = []

    for seq_id, seq_type, custom_dir in seq_list:
        path = get_csv_path(seq_id, args.base, custom_dir, TAG)

        if not os.path.exists(path):
            print(f"  MISSING: {path}")
            missing.append(seq_id)
            continue

        df = pd.read_csv(path)

        # Set residue order from the first file loaded
        if all_residues is None:
            all_residues = df['resSeq'].tolist()

        df_idx = df.set_index('resSeq')

        # ── R-score row ───────────────────────────────────────────────────
        r_row = {'seq_id': seq_id, 'seq_type': seq_type}
        for res in all_residues:
            if res in df_idx.index:
                val = df_idx.loc[res, 'R']
                r_row[f'R_{res}'] = val   # NaN preserved as-is
            else:
                r_row[f'R_{res}'] = np.nan

        # ── D and W rows ──────────────────────────────────────────────────
        dw_row = {'seq_id': seq_id, 'seq_type': seq_type}
        for res in all_residues:
            if res in df_idx.index:
                dw_row[f'D_{res}'] = df_idx.loc[res, 'D']
                dw_row[f'W_{res}'] = df_idx.loc[res, 'W']
            else:
                dw_row[f'D_{res}'] = 0.0
                dw_row[f'W_{res}'] = 0.0

        r_rows.append(r_row)
        dw_rows.append(dw_row)

        path_note = " [custom path]" if custom_dir else ""
        print(f"  Loaded {seq_id}  [{seq_type}]{path_note}")

    if not r_rows:
        raise RuntimeError("No sequences loaded. Check --base path and seq_ids.")

    n_seq = len(r_rows)
    n_res = len(all_residues)

    # ── Save R-score CSV ──────────────────────────────────────────────────
    r_df   = pd.DataFrame(r_rows)
    r_path = os.path.join(args.out_dir, f"r_scores_all_sequences_{TAG}.csv")
    r_df.to_csv(r_path, index=False, na_rep='')
    print(f"\nR-scores saved  -> {r_path}")
    print(f"  Shape: {r_df.shape}  "
          f"({n_seq} sequences x {n_res} R columns + 2 id columns)")

    # ── Save D/W CSV ──────────────────────────────────────────────────────
    dw_df   = pd.DataFrame(dw_rows)
    dw_path = os.path.join(args.out_dir, f"dw_scores_all_sequences_{TAG}.csv")
    dw_df.to_csv(dw_path, index=False)
    print(f"D/W scores saved -> {dw_path}")
    print(f"  Shape: {dw_df.shape}  "
          f"({n_seq} sequences x {n_res * 2} D/W columns + 2 id columns)")

    if missing:
        print(f"\nMissing ({len(missing)}): {missing}")
