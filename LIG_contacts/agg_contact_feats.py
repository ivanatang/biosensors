"""Collects per-sequence contact-summary CSVs into one merged feature table.

Run on the login node after all SLURM jobs finish. Merges
*_contact_summary_{TAG}.csv files into one table ready to merge with
feat_table.xlsx.

Usage:
    python aggregate_contact_features.py                  # full 40-500 ns
    python aggregate_contact_features.py --start-ns 40 --end-ns 250
"""

import os
import argparse
import pandas as pd

# ─────────────────────────────────────────────
# PATHS — mirror contact_type_analysis.py
# ─────────────────────────────────────────────
# Both per-sequence results and the combined output CSV live in the
# persistent repo location (not scratch), so they survive scratch's
# 90-day auto-deletion.
out_dir = "/projects/ivta1597/biosensors/LIG_contacts"

# ─────────────────────────────────────────────
# ARGUMENTS
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--seq_list', default="/projects/ivta1597/biosensors/seq_ids_orig.txt",
                    help='Path to seq_ids.txt-style sequence list, used only for the '
                         'coverage check (default: seq_ids_orig.txt in the repo root)')
parser.add_argument('--start-ns', type=float, default=40.0,
                    help='Start of analysis window in ns (default: 40)')
parser.add_argument('--end-ns',   type=float, default=500.0,
                    help='End of analysis window in ns (default: 500)')
parser.add_argument('--ligand-region', choices=['whole', 'core', 'tail'], default='whole',
                    help='Ligand region the underlying contact_type_analysis.py runs '
                         'were restricted to (default: whole)')
parser.add_argument('--suffix', default='',
                    help="Run-directory/output suffix, e.g. '_qfix', matching "
                         "whatever contact_type_analysis.py was run with for these "
                         "sequences (default: '', the standard directory). Also "
                         "appended to this script's own output filename so a _qfix "
                         "run never overwrites the standard one.")
args = parser.parse_args()
seq_ids_file = args.seq_list

TAG         = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
REGION_TAG  = "" if args.ligand_region == "whole" else f"_{args.ligand_region}"
results_dir = os.path.join(out_dir, f"contact_type_results_{TAG}{REGION_TAG}{args.suffix}")
out_path    = os.path.join(out_dir, f"contact_features_all_{TAG}{REGION_TAG}{args.suffix}.csv")
os.makedirs(out_dir, exist_ok=True)

print(f"Window      : {args.start_ns:.0f}-{args.end_ns:.0f} ns  (tag: {TAG})")
print(f"Region      : {args.ligand_region}")
print(f"Results dir : {results_dir}")
print(f"Output dir  : {out_dir}")

# ─────────────────────────────────────────────
# LOAD -- explicit per-sequence paths from seq_list, not a glob over the
# whole shared results dir. A glob also picks up leftover output from any
# sequence that ever had this step run, including ones no longer in the
# current cohort (e.g. seq_ids_orig.txt's superseded 200-sequence list),
# inflating the combined table with unlabeled rows. See the equivalent fix
# in agg_gate_latch_water_bridge.py / agg_residue_atom_split.py.
# ─────────────────────────────────────────────
with open(seq_ids_file) as fh:
    all_ids = [l.split()[0] for l in fh if l.strip() and not l.startswith('#')]

dfs     = []
missing = []
for seq_id in all_ids:
    path = os.path.join(results_dir, f"{seq_id}_contact_summary_{TAG}{REGION_TAG}.csv")
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        missing.append(seq_id)
        continue
    dfs.append(pd.read_csv(path))

print(f"Found {len(dfs)} of {len(all_ids)} summary files")
if not dfs:
    raise FileNotFoundError(
        f"No summary CSVs found for any sequence in {seq_ids_file}.\n"
        f"Expected pattern: {results_dir}/<seq_id>_contact_summary_{TAG}{REGION_TAG}.csv"
    )

combined = pd.concat(dfs, ignore_index=True)

if missing:
    print(f"\nWARNING: {len(missing)} sequences missing results: {missing}")
else:
    print("All sequences accounted for.")

print(f"\nFeature table shape: {combined.shape}")
print(combined.to_string())

combined.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
