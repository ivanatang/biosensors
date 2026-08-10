"""
extract_rmsf_feats.py

Extracts per-sequence RMSF features for the regions of interest (gate, latch,
Lb7a5 loop, C-terminal recoil helix, plus individual pocket residues) from
gmx rmsf outputs, for every sequence in seq_ids.txt. Reads directly from the
PetaLibrary archive (BASE, matching config.yaml's paths.base and every
sibling extraction script) and writes into REPO_DIR/analysis/, so this runs
entirely on Alpine's login node -- no local Mac / OneDrive sync needed.

Two tables are produced:
    rmsf_single_residues_per_seq{TAG}.csv
        Per-residue RMSF for single pocket residues of interest.
        Source: rmsf_PL.xvg (full per-residue RMSF)

    rmsf_ca_per_seq_summary{TAG}.csv
        Mean +/- SD RMSF per structural region (Ca atoms only).
        Source: rmsf_PL_ca_{gate,latch,Lb7a5,recoil}.xvg

TAG defaults to "_500ns" and is set once at the top of the CONFIG block
(or via --tag), so re-running against a different analysis window only
requires changing that one value. --end-ns controls which windowed
subdirectory (e.g. "250ns") the source .xvg files are read from, and is
independent of --tag (which only affects output filenames).

Usage:
    python extract_rmsf_feats.py [seq_ids.txt] [--tag _500ns] [--end-ns 500]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

# ── Configurable paths / window tag ───────────────────────────────────────────
# BASE mirrors config.yaml's paths.base and every sibling script
# (contact_type_analysis.py, extract_gate_latch_rmsd_feats.py,
# salt_bridge_analysis.py): rmsf_PL*.xvg is written here by
# post_processing_pipeline_worker.sh, on the durable PetaLibrary archive, not
# scratch (auto-deletes after 90 days) or OneDrive. Runs on Alpine's login
# node directly, no local Mac / OneDrive sync step needed.
BASE     = "/pl/active/shirts_archive/IvanaTang/biosensors"
REPO_DIR = "/projects/ivta1597/biosensors"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
TAG    = "_500ns"   # appended to output CSV filenames only
NM_TO_ANG = 10.0
# ─────────────────────────────────────────────────────────────────────────────

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

REGIONS_CA = {
    "Gate (r84-90)":             "rmsf_PL_ca_gate.xvg",
    "Latch (r114-118)":          "rmsf_PL_ca_latch.xvg",
    "Lb7a5 (r148-155)":          "rmsf_PL_ca_Lb7a5.xvg",
    "C-term helix (r154-166)":   "rmsf_PL_ca_recoil.xvg",
}

SINGLE_RESIDUES = {
    "Q69":  69,
    "I134": 134,
    "Y23":  23,
    "K59":  59,
    "R79":  79,
    "I110": 110,
    "G163": 163,
}


def rmsf_run_dir(seq_id, group_label, end_ns):
    """Directory containing a sequence's rmsf_*.xvg outputs. Newer pipeline
    runs write outputs under runrel/{end_ns}ns/; older ones write directly
    into runrel/. Prefer whichever location actually has rmsf_PL.xvg rather
    than guessing from seq_id, since the two layouts don't map cleanly onto
    naming prefix or pair ID (e.g. pair_0482_low_pkt uses the flat layout
    despite matching the "resubmitted with new pipeline" ID list)."""
    run_dir = os.path.join(BASE, TYPE_SUBDIR[group_label], seq_id, RUNREL)
    nested_dir = os.path.join(run_dir, f"{int(end_ns)}ns")
    if os.path.exists(os.path.join(nested_dir, "rmsf_PL.xvg")):
        return nested_dir
    return run_dir


def get_data(filepath):
    """Parse a GROMACS .xvg file, skipping comment/annotation lines."""
    x_data, y_data = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            cols = line.split()
            if len(cols) >= 2:
                x_data.append(float(cols[0]))
                y_data.append(float(cols[1]))
    return np.array(x_data), np.array(y_data)


def load_seq_ids(seq_list_path):
    """Yields (seq_id, group_label) tuples from a seq_ids.txt-style file."""
    with open(seq_list_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            seq_id     = parts[0].strip()
            group_label = parts[1].strip() if len(parts) > 1 else ""
            yield seq_id, group_label


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?", default="seq_ids.txt")
    parser.add_argument("--tag", default=TAG,
                        help=f"Suffix appended to output CSV filenames (default: {TAG})")
    parser.add_argument("--end-ns", type=float, default=500.0,
                        help="End of analysis window in ns; selects which "
                             "windowed subdirectory to read .xvg files from "
                             "(default: 500)")
    args = parser.parse_args()
    tag = args.tag
    end_ns = args.end_ns

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    all_systems = list(load_seq_ids(args.seq_list))

    # ── 1. Single-residue RMSF (from full per-residue rmsf_PL.xvg) ───────────
    rows_single, missing_single = [], []
    for seq_id, group_label in all_systems:
        xvg = os.path.join(rmsf_run_dir(seq_id, group_label, end_ns), "rmsf_PL.xvg")
        if not os.path.exists(xvg):
            missing_single.append(seq_id)
            continue

        res_id, rmsf_nm = get_data(xvg)
        res_id  = res_id.astype(int)
        rmsf_ang = rmsf_nm * NM_TO_ANG

        row = {"Sequence": seq_id, "Group": group_label}
        for label, resnum in SINGLE_RESIDUES.items():
            vals = rmsf_ang[res_id == resnum]
            row[f"{label} RMSF (A)"] = float(vals[0]) if vals.size > 0 else np.nan
        rows_single.append(row)

    if missing_single:
        print(f"rmsf_PL.xvg missing for {len(missing_single)} sequences: "
              f"{missing_single[:5]}{'...' if len(missing_single) > 5 else ''}")

    if not rows_single:
        print("\nNo rmsf_PL.xvg files found - nothing to write.")
        sys.exit(1)

    single_res_df = pd.DataFrame(rows_single).set_index("Sequence")
    single_res_df = single_res_df.round(
        {f"{label} RMSF (A)": 5 for label in SINGLE_RESIDUES}
    )

    out_dir = os.path.join(REPO_DIR, "analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path_single = os.path.join(out_dir, f"rmsf_single_residues_per_seq{tag}.csv")
    single_res_df.to_csv(csv_path_single)
    print(f"Saved ({len(single_res_df)} rows) -> {csv_path_single}")

    # ── 2. Per-region Ca RMSF summary (mean +/- SD) ──────────────────────────
    rows_ca, missing_ca, empty_ca = [], [], []
    for seq_id, group_label in all_systems:
        run_dir = rmsf_run_dir(seq_id, group_label, end_ns)
        row = {"Sequence": seq_id, "Group": group_label}
        any_found = False
        for region_name, region_fname in REGIONS_CA.items():
            xvg = os.path.join(run_dir, region_fname)
            if not os.path.exists(xvg):
                row[f"{region_name} mean (A)"] = np.nan
                row[f"{region_name} SD (A)"]   = np.nan
                continue
            any_found = True
            _, rmsf_nm = get_data(xvg)
            if rmsf_nm.size == 0:
                # File exists (gmx wrote it) but has zero data rows, i.e. the
                # underlying index group had zero atoms for this sequence --
                # NaN, not an np.mean()-on-empty-array warning with no context.
                empty_ca.append((seq_id, region_name, xvg))
                row[f"{region_name} mean (A)"] = np.nan
                row[f"{region_name} SD (A)"]   = np.nan
                continue
            vals = rmsf_nm * NM_TO_ANG
            row[f"{region_name} mean (A)"] = round(float(np.mean(vals)), 5)
            row[f"{region_name} SD (A)"]   = (
                round(float(np.std(vals, ddof=1)), 5) if len(vals) > 1 else np.nan
            )
        if not any_found:
            missing_ca.append(seq_id)
            continue
        rows_ca.append(row)

    if missing_ca:
        print(f"rmsf_PL_ca_*.xvg missing for {len(missing_ca)} sequences: "
              f"{missing_ca[:5]}{'...' if len(missing_ca) > 5 else ''}")

    if empty_ca:
        print(f"\nrmsf_PL_ca_*.xvg present but EMPTY (zero data rows) for "
              f"{len(empty_ca)} (sequence, region) pair(s):")
        for seq_id, region_name, xvg in empty_ca:
            print(f"  {seq_id:<28} {region_name:<22} {xvg}")

    if not rows_ca:
        print("\nNo rmsf_PL_ca_*.xvg files found - nothing to write.")
        sys.exit(1)

    df_ca = pd.DataFrame(rows_ca).set_index("Sequence")
    csv_path_ca = os.path.join(out_dir, f"rmsf_ca_per_seq_summary{tag}.csv")
    df_ca.to_csv(csv_path_ca)
    print(f"Saved ({len(df_ca)} rows) -> {csv_path_ca}")


if __name__ == "__main__":
    main()
