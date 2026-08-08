#!/usr/bin/env python3
"""
parse_gate_latch_hbond.py
----------------------------
Cross-references gmx hbond's real (distance + angle) H-bond counts, from
gate_latch_hbond_gmx.sh, against the specific frames where
gate_latch_water_bridge.py (--dump-frames) found the literal single-water
gate-latch-ligand triple bridge active.

This is the join that actually answers "when the documented network is
forming, is a genuine hydrogen bond present, and does it land on the
backbone or side chain" -- neither script alone answers it: the heavy-atom
4 A test in gate_latch_water_bridge.py doesn't check H-bond geometry, and
gmx hbond alone can't identify a single specific bridging water so its raw
counts also include ordinary surface hydration unrelated to the network.

For each of gate_backbone / gate_sidechain / latch_backbone / latch_sidechain,
this reports, separately for triple-bridge-active frames and inactive frames:
  frac_hbond_present : fraction of frames with >=1 real H-bond (gmx hbond
                        count > 0) between water_sol and that residue-part
  mean_hbond_count    : mean H-bond count over those frames

A large gap between the triple-bridge and non-triple-bridge fractions for a
given part (backbone or side chain) is the actual evidence for which part
the network's H-bond is going through -- not just which atom happens to be
geometrically closest in a single frame.

Usage
-----
    conda activate biosensors
    python parse_gate_latch_hbond.py --seq_id pair_3085_binder --seq_type binders \
        --start-ns 40 --end-ns 500
"""

import os
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--seq_id',   required=True)
parser.add_argument('--seq_type', required=True)
parser.add_argument('--start-ns', type=float, default=40.0)
parser.add_argument('--end-ns',   type=float, default=500.0)
parser.add_argument('--frames-start-ns', type=float, default=0.0,
                    help="start-ns used when the --dump-frames CSV from "
                         "gate_latch_water_bridge.py was generated (that "
                         "script defaults to 0, not 40 -- must match the "
                         "TAG in the frames CSV filename actually on disk).")
parser.add_argument('--frames-end-ns', type=float, default=500.0)
parser.add_argument('--ligand-region', choices=['whole', 'core', 'tail'], default='core',
                    help="Must match the --ligand-region the --dump-frames "
                         "CSV was generated with (affects its filename).")
parser.add_argument('--join-tol-ps', type=float, default=200.0,
                    help="Max time gap (ps) allowed when matching a gmx "
                         "hbond frame to a triple-bridge frame -- generous "
                         "default since the two scripts may use different "
                         "strides (gmx hbond typically runs every raw frame "
                         "at 37.5 ps spacing; gate_latch_water_bridge.py "
                         "uses STRIDE=10, i.e. 375 ps spacing).")
args = parser.parse_args()
seq_id, seq_type = args.seq_id, args.seq_type
ligand_region = args.ligand_region

HBOND_TAG  = f"{int(args.start_ns)}_{int(args.end_ns)}ns"
FRAMES_TAG = f"{int(args.frames_start_ns)}_{int(args.frames_end_ns)}ns"
REGION_TAG = "" if ligand_region == "whole" else f"_{ligand_region}"

base = "/scratch/alpine/ivta1597/LCA_boltz_models"
prod = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
rundir = os.path.join(base, seq_type, seq_id, prod)

frames_csv = os.path.join(
    base, seq_type, seq_id,
    f"gate_latch_water_bridge_{FRAMES_TAG}{REGION_TAG}",
    f"{seq_id}_gate_latch_bridge_frames_{FRAMES_TAG}{REGION_TAG}.csv",
)

GROUPS = ["gate_backbone", "gate_sidechain", "latch_backbone", "latch_sidechain"]


def read_hbond_num_xvg(path):
    """gmx hbond -num output: time (ps) in col 0, H-bond count in col 1.
    Skips comment/metadata lines per this repo's standard .xvg convention."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(('#', '@')):
                continue
            parts = line.split()
            rows.append((float(parts[0]), float(parts[1])))
    df = pd.DataFrame(rows, columns=["time_ps", "hbond_count"])
    df["time_ns"] = df["time_ps"] / 1000.0
    return df[["time_ns", "hbond_count"]]


def main():
    if not os.path.exists(frames_csv):
        raise FileNotFoundError(
            f"Frames CSV not found: {frames_csv}\n"
            f"Re-run gate_latch_water_bridge.py with --dump-frames "
            f"--start-ns {args.frames_start_ns} --end-ns {args.frames_end_ns} "
            f"--ligand-region {ligand_region} first."
        )
    frames = pd.read_csv(frames_csv).sort_values("time_ns")
    print(f"[{seq_id}] Loaded {len(frames)} triple-bridge frames "
          f"({frames['triple_bridge'].sum()} active)")

    summary_rows = []

    for grp in GROUPS:
        xvg_path = os.path.join(
            rundir, f"hbond_{grp}_{int(args.start_ns)}_{int(args.end_ns)}ns_num.xvg"
        )
        if not os.path.exists(xvg_path):
            print(f"[{seq_id}] WARNING: missing {xvg_path}, skipping {grp}")
            continue

        hb = read_hbond_num_xvg(xvg_path).sort_values("time_ns")

        joined = pd.merge_asof(
            frames, hb, on="time_ns", direction="nearest",
            tolerance=args.join_tol_ps / 1000.0,
        )
        n_unmatched = joined["hbond_count"].isna().sum()
        if n_unmatched:
            print(f"[{seq_id}] {grp}: {n_unmatched} frames had no gmx hbond "
                  f"sample within {args.join_tol_ps} ps -- dropped")
        joined = joined.dropna(subset=["hbond_count"])

        for label, sub in [("triple_bridge_active", joined[joined["triple_bridge"]]),
                            ("triple_bridge_inactive", joined[~joined["triple_bridge"]])]:
            if len(sub) == 0:
                continue
            summary_rows.append(dict(
                seq_id=seq_id,
                group=grp,
                condition=label,
                n_frames=len(sub),
                frac_hbond_present=round(float((sub["hbond_count"] > 0).mean()), 4),
                mean_hbond_count=round(float(sub["hbond_count"].mean()), 4),
            ))

    df_out = pd.DataFrame(summary_rows)
    out_path = os.path.join(rundir, f"{seq_id}_gate_latch_hbond_crossref_{HBOND_TAG}.csv")
    df_out.to_csv(out_path, index=False)

    print(f"\n-- {seq_id}: real H-bond presence, backbone vs side chain -----------")
    print(df_out.to_string(index=False))
    print(f"\nWrote: {out_path}")

    active = df_out[df_out["condition"] == "triple_bridge_active"]
    if not active.empty:
        gate = active[active["group"].str.startswith("gate_")]
        latch = active[active["group"].str.startswith("latch_")]
        for label, sub in [("Gate", gate), ("Latch", latch)]:
            if len(sub) < 2:
                continue
            bb = sub[sub["group"].str.endswith("backbone")]["frac_hbond_present"]
            sc = sub[sub["group"].str.endswith("sidechain")]["frac_hbond_present"]
            if len(bb) and len(sc):
                print(f"{label}: backbone H-bond present in {float(bb.iloc[0])*100:.1f}% "
                      f"of active frames, side chain in {float(sc.iloc[0])*100:.1f}%")


if __name__ == "__main__":
    main()
