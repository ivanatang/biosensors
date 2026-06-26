"""
contact_type_analysis.py
------------------------
Compute protein–ligand contact type features (hydrophobic, polar, charged)
for a single PYR1 sequence trajectory on Alpine HPC.

seq_id is always passed as a positional CLI argument — either directly
for local testing, or by contact_type_worker.sh for SLURM array runs.

Usage:
    python contact_type_analysis.py <seq_id>

Output:
    <output_dir>/<seq_id>_contact_perframe.csv    — per-frame contact counts
    <output_dir>/<seq_id>_residue_occupancy.csv   — per-residue occupancy
    <output_dir>/<seq_id>_contact_summary.csv     — scalar features for feat_table
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

# ─────────────────────────────────────────────
# PATHS  — edit these for your environment
# ─────────────────────────────────────────────
base    = "/scratch/alpine/ivta1597/LCA_boltz_models"
runrel  = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

output_dir = os.path.join(base, "LIG_contacts_flex/contact_type_results")   # default; overridden by TAG in main()

# ─────────────────────────────────────────────
# TYPE SUBDIRECTORY MAP
# ─────────────────────────────────────────────
TYPE_SUBDIR = {
    "binder":    "binders",
    "nb":        "nonbinders",
    "low_pkt":   "neg_low_pkt",
    "fail_gate": "neg_fail_gate",
}

def get_type_subdir(seq_id):
    """
    Extract the type suffix from seq_id and return the corresponding
    subdirectory. Checks two-token suffixes (low_pkt, fail_gate) before
    single-token ones (binder, nb) to avoid partial matches.

    Examples:
        pair_3059_binder    -> binders
        pair_0664_nb        -> nonbinders
        pair_0008_low_pkt   -> neg_low_pkt
        pair_0715_fail_gate -> neg_fail_gate
    """
    tokens = seq_id.split("_")

    # Check two-token suffix first
    if len(tokens) >= 2:
        two_token = "_".join(tokens[-2:])
        if two_token in TYPE_SUBDIR:
            return TYPE_SUBDIR[two_token]

    # Fall back to single-token suffix
    if tokens[-1] in TYPE_SUBDIR:
        return TYPE_SUBDIR[tokens[-1]]

    raise ValueError(
        f"Cannot determine type subdirectory from seq_id '{seq_id}'. "
        f"Expected suffix from: {list(TYPE_SUBDIR.keys())}"
    )

# ─────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────
LIGAND_RESNAME  = "LIG"
CUTOFF_NM       = 0.45     # 4.5 Å heavy-atom contact cutoff
STRIDE          = 1        # set >1 only for quick debugging; use 1 for publication

# ─────────────────────────────────────────────
# RESIDUE CHEMISTRY MAP
# ─────────────────────────────────────────────
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR       = {"SER", "THR", "CYS", "TYR", "ASN", "GLN"}
POS_CHARGED = {"ARG", "LYS", "HIS"}
NEG_CHARGED = {"ASP", "GLU"}

def classify_residue(resname):
    if resname in HYDROPHOBIC:  return "hydrophobic"
    if resname in POLAR:        return "polar"
    if resname in POS_CHARGED:  return "pos_charged"
    if resname in NEG_CHARGED:  return "neg_charged"
    return "other"


# ─────────────────────────────────────────────
# ARGUMENTS
# ─────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute protein–ligand contact type features for one sequence."
    )
    parser.add_argument("seq_id",
                        help="Sequence identifier (e.g. pair_3059_binder)")
    parser.add_argument("--start-ns", type=float, default=40.0,
                        help="Start of analysis window in ns (default: 40)")
    parser.add_argument("--end-ns",   type=float, default=500.0,
                        help="End of analysis window in ns (default: 500)")
    return parser.parse_args()


# ─────────────────────────────────────────────
# LOAD TRAJECTORY
# ─────────────────────────────────────────────
def load_trajectory(seq_id, start_ps, end_ps):
    type_subdir = get_type_subdir(seq_id)
    seq_dir  = os.path.join(base, type_subdir, seq_id, runrel)
    top_path = os.path.join(seq_dir, "medoid_PL.pdb")
    xtc_path = os.path.join(seq_dir, "PL_only_40_500ns.xtc")

    for path in (top_path, xtc_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Expected file not found: {path}")

    print(f"[{seq_id}] Loading topology: {top_path}")
    print(f"[{seq_id}] Loading trajectory: {xtc_path}  (stride={STRIDE})")

    traj = md.load(xtc_path, top=top_path, stride=STRIDE)
    print(f"[{seq_id}] Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")

    # ── Restrict to analysis window ────────────────────────────────────
    mask = (traj.time >= start_ps) & (traj.time <= end_ps)
    traj = traj[mask]
    print(f"[{seq_id}] Time window: {start_ps/1000:.0f}–{end_ps/1000:.0f} ns  "
          f"({mask.sum()} frames retained)")
    # ──────────────────────────────────────────────────────────────────

    return traj


# ─────────────────────────────────────────────
# CONTACT TYPE COMPUTATION
# ─────────────────────────────────────────────
def compute_contact_type_features(traj, seq_id):
    """
    For each protein residue within CUTOFF_NM of LCA (any frame),
    compute per-frame binary contact and return:
      - per_frame_df : (n_frames) counts of each contact type per frame
      - residue_df   : (n_residues) occupancy per residue, for inspection
    """
    top = traj.topology

    # ── Ligand heavy atoms ──────────────────────────────────────────────
    lig_atoms = [
        a.index for a in top.atoms
        if a.residue.name == LIGAND_RESNAME and a.element.symbol != "H"
    ]
    if not lig_atoms:
        raise ValueError(
            f"[{seq_id}] No atoms found with residue name '{LIGAND_RESNAME}'. "
            f"Check LIGAND_RESNAME. Available residue names: "
            f"{sorted({r.name for r in top.residues})}"
        )
    print(f"[{seq_id}] Ligand '{LIGAND_RESNAME}': {len(lig_atoms)} heavy atoms")

    # ── Protein residues ─────────────────────────────────────────────────
    prot_residues = [
        (r.index, r.name, classify_residue(r.name))
        for r in top.residues if r.is_protein
    ]
    print(f"[{seq_id}] Protein residues: {len(prot_residues)}")

    # ── Per-residue contact calculation ──────────────────────────────────
    records = []
    for res_idx, res_name, res_class in prot_residues:
        prot_atoms = [
            a.index for a in top.residue(res_idx).atoms
            if a.element.symbol != "H"
        ]
        if not prot_atoms:
            continue

        # All heavy-atom pairs between this residue and LCA
        pairs = np.array(
            [[pa, la] for pa in prot_atoms for la in lig_atoms],
            dtype=int
        )

        # Minimum distance to LCA per frame  →  shape (n_frames,)
        dists    = md.compute_distances(traj, pairs)   # (n_frames, n_pairs)
        min_dist = dists.min(axis=1)
        in_contact = (min_dist < CUTOFF_NM).astype(np.int8)

        records.append({
            "res_idx":    res_idx,
            "res_name":   res_name,
            "res_class":  res_class,
            "occupancy":  in_contact.mean(),
            "in_contact": in_contact,        # (n_frames,) — kept for summing
        })

    residue_df = pd.DataFrame([
        {k: v for k, v in r.items() if k != "in_contact"}
        for r in records
    ])

    # ── Per-frame type counts ─────────────────────────────────────────────
    contact_matrix = np.vstack([r["in_contact"] for r in records])  # (n_res, n_frames)
    type_labels    = np.array([r["res_class"] for r in records])

    contact_types = ["hydrophobic", "polar", "pos_charged", "neg_charged", "other"]
    per_frame = {}
    for t in contact_types:
        mask = (type_labels == t)
        per_frame[f"n_{t}"] = contact_matrix[mask].sum(axis=0).astype(int)

    total = sum(per_frame[f"n_{t}"] for t in contact_types)
    per_frame["n_total"] = total
    per_frame["frac_hydrophobic"] = np.where(
        total > 0, per_frame["n_hydrophobic"] / total, 0.0
    )

    per_frame_df = pd.DataFrame(per_frame)
    per_frame_df.index.name = "frame"

    return per_frame_df, residue_df


# ─────────────────────────────────────────────
# SUMMARISE → SCALAR FEATURES FOR feat_table
# ─────────────────────────────────────────────
def summarise_features(per_frame_df, seq_id):
    """
    Collapse per-frame arrays into scalar features.
    Occupancy columns: fraction of frames with ≥1 contact of that type.
    """
    rows = {"seq_id": seq_id}

    numeric_cols = [c for c in per_frame_df.columns
                    if c not in ("frame",)]

    for col in numeric_cols:
        arr = per_frame_df[col].values
        rows[f"mean_{col}"] = arr.mean()
        rows[f"std_{col}"]  = arr.std()
        if col.startswith("n_"):
            rows[f"occ_{col}_gt0"] = (arr > 0).mean()

    return pd.DataFrame([rows])


# ─────────────────────────────────────────────
# DIAGNOSTIC PRINT
# ─────────────────────────────────────────────
def print_diagnostics(per_frame_df, residue_df, seq_id):
    print(f"\n── Contact summary for {seq_id} ─────────────────────────")
    print(f"  Frames analysed : {len(per_frame_df)}")
    print(f"  Mean total contacts     : {per_frame_df['n_total'].mean():.2f}")
    print(f"  Mean hydrophobic        : {per_frame_df['n_hydrophobic'].mean():.2f}")
    print(f"  Mean polar              : {per_frame_df['n_polar'].mean():.2f}")
    print(f"  Mean pos_charged        : {per_frame_df['n_pos_charged'].mean():.2f}")
    print(f"  Mean neg_charged        : {per_frame_df['n_neg_charged'].mean():.2f}")
    print(f"  Mean frac_hydrophobic   : {per_frame_df['frac_hydrophobic'].mean():.3f}")

    top_contacts = (
        residue_df[residue_df["occupancy"] > 0.10]
        .sort_values("occupancy", ascending=False)
        .head(10)
    )
    print(f"\n  Top residues by occupancy (>10%):")
    print(top_contacts[["res_name", "res_class", "occupancy"]].to_string(index=False))
    print("─" * 55)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    args     = parse_args()
    seq_id   = args.seq_id
    start_ns = args.start_ns
    end_ns   = args.end_ns
    start_ps = int(start_ns * 1000)
    end_ps   = int(end_ns   * 1000)
    TAG      = f"{int(start_ns)}_{int(end_ns)}ns"

    # Tagged output directory — one per time window, consistent with
    # water_contacts_{TAG}/ naming used by the water analysis pipeline
    tagged_output_dir = os.path.join(base, f"LIG_contacts_flex/contact_type_results_{TAG}")
    os.makedirs(tagged_output_dir, exist_ok=True)

    print(f"\n=== contact_type_analysis.py  |  seq_id={seq_id}  "
          f"window={start_ns:.0f}–{end_ns:.0f} ns ===")

    # Check if already done (safe to re-queue without reprocessing)
    summary_out = os.path.join(tagged_output_dir,
                               f"{seq_id}_contact_summary_{TAG}.csv")
    if os.path.exists(summary_out):
        print(f"[{seq_id}] Output already exists, skipping: {summary_out}")
        sys.exit(0)

    traj = load_trajectory(seq_id, start_ps, end_ps)

    per_frame_df, residue_df = compute_contact_type_features(traj, seq_id)

    print_diagnostics(per_frame_df, residue_df, seq_id)

    # ── Write outputs ─────────────────────────────────────────────────────
    per_frame_out = os.path.join(tagged_output_dir,
                                 f"{seq_id}_contact_perframe_{TAG}.csv")
    residue_out   = os.path.join(tagged_output_dir,
                                 f"{seq_id}_residue_occupancy_{TAG}.csv")

    per_frame_df.to_csv(per_frame_out)
    residue_df.to_csv(residue_out, index=False)
    summarise_features(per_frame_df, seq_id).to_csv(summary_out, index=False)

    print(f"[{seq_id}] Wrote: {per_frame_out}")
    print(f"[{seq_id}] Wrote: {residue_out}")
    print(f"[{seq_id}] Wrote: {summary_out}")
    print(f"[{seq_id}] Done.")


if __name__ == "__main__":
    main()
