"""
aggregate_contact_features.py
------------------------------
After all SLURM jobs finish, run this on the login node to collect
per-sequence *_contact_summary_{TAG}.csv files into one table ready to
merge with feat_table.xlsx.

Usage:
    python aggregate_contact_features.py                  # full 40-500 ns
    python aggregate_contact_features.py --start-ns 40 --end-ns 250
"""

import os
import glob
import argparse
import pandas as pd

# ─────────────────────────────────────────────
# PATHS — mirror contact_type_analysis.py
# ─────────────────────────────────────────────
base_contacts = "/scratch/alpine/ivta1597/LCA_boltz_models/LIG_contacts_flex"
seq_ids_file  = os.path.join(base_contacts, "seq_ids.txt")

# ─────────────────────────────────────────────
# ARGUMENTS
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--start-ns', type=float, default=40.0,
                    help='Start of analysis window in ns (default: 40)')
parser.add_argument('--end-ns',   type=float, default=500.0,
                    help='End of analysis window in ns (default: 500)')
args = parser.parse_args()

TAG         = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
results_dir = os.path.join(base_contacts, f"contact_type_results_{TAG}")
out_path    = os.path.join(results_dir, f"contact_features_all_{TAG}.csv")

print(f"Window      : {args.start_ns:.0f}-{args.end_ns:.0f} ns  (tag: {TAG})")
print(f"Results dir : {results_dir}")

# ─────────────────────────────────────────────
# LOAD
# ─────────────────────────────────────────────
summary_files = sorted(glob.glob(
    os.path.join(results_dir, f"*_contact_summary_{TAG}.csv")
))
print(f"Found {len(summary_files)} summary files")

if not summary_files:
    raise FileNotFoundError(
        f"No summary CSVs found in: {results_dir}\n"
        f"Expected pattern: *_contact_summary_{TAG}.csv"
    )

dfs      = [pd.read_csv(f) for f in summary_files]
combined = pd.concat(dfs, ignore_index=True)

# ─────────────────────────────────────────────
# CHECK COVERAGE vs seq_ids.txt
# ─────────────────────────────────────────────
if os.path.exists(seq_ids_file):
    with open(seq_ids_file) as fh:
        all_ids = [l.split()[0] for l in fh if l.strip()]
    missing = set(all_ids) - set(combined["seq_id"])
    if missing:
        print(f"WARNING: {len(missing)} sequences missing results: {missing}")
    else:
        print("All sequences accounted for.")

print(f"\nFeature table shape: {combined.shape}")
print(combined.to_string())

combined.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
