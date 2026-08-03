"""
residue_atom_split_contact.py
------------------------------
For a single protein residue (default: 116, the Arg on the latch flagged
by Ryan), split its heavy atoms into backbone (N, CA, C, O) vs side chain
and compute each group's contact occupancy with the ligand separately.

Motivation: contact_type_analysis.py and R_score_calc.py both test contact
at the whole-residue level (any heavy atom on the residue within cutoff of
the ligand). That can't distinguish "this residue's backbone sits near the
ligand" from "this residue's side chain reaches out to the ligand", which
matters for residue 116: if its Arg side chain points away from the ligand
in the Hab1-bound state (per Ryan), then a core-preferential contact
signal at 116 is more likely a backbone/latch-positioning effect than a
side-chain-ligand interaction.

seq_id is always passed as a positional CLI argument, either directly for
local testing or by residue_atom_split_worker.sh for SLURM array runs.

Usage:
    python residue_atom_split_contact.py <seq_id> [--resseq 116]
        [--ligand-region core] [--start-ns 40] [--end-ns 500]

Output:
    <output_dir>/<seq_id>_res<RESSEQ>_atomsplit_{TAG}{REGION_TAG}.csv
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

# ─────────────────────────────────────────────
# PATHS  — mirror contact_type_analysis.py
# ─────────────────────────────────────────────
base    = "/pl/active/shirts_archive/IvanaTang/biosensors"
runrel  = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
output_base = "/projects/ivta1597/biosensors/LIG_contacts"

TYPE_SUBDIR = {
    "binder":    "binders",
    "nb":        "nonbinders",
    "low_pkt":   "neg_low_pkt",
    "fail_gate": "neg_fail_gate",
}

def get_type_subdir(seq_id):
    tokens = seq_id.split("_")
    if len(tokens) >= 2:
        two_token = "_".join(tokens[-2:])
        if two_token in TYPE_SUBDIR:
            return TYPE_SUBDIR[two_token]
    if tokens[-1] in TYPE_SUBDIR:
        return TYPE_SUBDIR[tokens[-1]]
    raise ValueError(
        f"Cannot determine type subdirectory from seq_id '{seq_id}'. "
        f"Expected suffix from: {list(TYPE_SUBDIR.keys())}"
    )


def run_dir(seq_id):
    """Mirrors contact_type_analysis.py's run_dir(): prefer the nested
    500ns/ subdirectory if that's where medoid_PL.pdb actually lives."""
    flat_dir   = os.path.join(base, get_type_subdir(seq_id), seq_id, runrel)
    nested_dir = os.path.join(flat_dir, "500ns")
    if os.path.exists(os.path.join(nested_dir, "medoid_PL.pdb")):
        return nested_dir
    return flat_dir


# ─────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────
LIGAND_RESNAME  = "LIG"
CUTOFF_NM       = 0.45     # 4.5 Å heavy-atom contact cutoff, matches contact_type_analysis.py
STRIDE          = 1

BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O"}

# Same core/tail ligand-atom split as contact_type_analysis.py / R_score_calc.py.
LIGAND_TAIL_ORDINALS     = {0, 1, 3, 4, 7, 9, 19}   # O41,O42,C44,C45,C48,C50,C60
LIGAND_EXPECTED_ELEMENTS = ['O', 'O', 'O'] + ['C'] * 24


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split one residue's ligand contact into backbone vs side chain."
    )
    parser.add_argument("seq_id", help="Sequence identifier (e.g. pair_3059_binder)")
    parser.add_argument("--resseq", type=int, default=116,
                        help="PDB residue number to analyze (default: 116)")
    parser.add_argument("--start-ns", type=float, default=40.0)
    parser.add_argument("--end-ns",   type=float, default=500.0)
    parser.add_argument("--ligand-region", choices=["whole", "core", "tail"], default="core",
                        help="Ligand atoms to test contact against (default: core, "
                             "since this is testing the core-preferential finding at 116)")
    return parser.parse_args()


def load_trajectory(seq_id, start_ps, end_ps):
    seq_dir  = run_dir(seq_id)
    top_path = os.path.join(seq_dir, "medoid_PL.pdb")
    xtc_path = os.path.join(seq_dir, "PL_only_40_500ns.xtc")

    for path in (top_path, xtc_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Expected file not found: {path}")

    print(f"[{seq_id}] Loading topology: {top_path}")
    print(f"[{seq_id}] Loading trajectory: {xtc_path}  (stride={STRIDE})")

    traj = md.load(xtc_path, top=top_path, stride=STRIDE)
    mask = (traj.time >= start_ps) & (traj.time <= end_ps)
    traj = traj[mask]
    print(f"[{seq_id}] Time window: {start_ps/1000:.0f}-{end_ps/1000:.0f} ns  "
          f"({mask.sum()} frames retained)")
    return traj


def get_ligand_atoms(top, seq_id, ligand_region):
    lig_heavy_atoms_all = [
        a for a in top.atoms
        if a.residue.name == LIGAND_RESNAME and a.element.symbol != "H"
    ]
    if not lig_heavy_atoms_all:
        raise ValueError(
            f"[{seq_id}] No atoms found with residue name '{LIGAND_RESNAME}'. "
            f"Available residue names: {sorted({r.name for r in top.residues})}"
        )

    if ligand_region != "whole":
        got_elements = [a.element.symbol for a in lig_heavy_atoms_all]
        if got_elements != LIGAND_EXPECTED_ELEMENTS:
            raise ValueError(
                f"[{seq_id}] Ligand heavy-atom element order doesn't match the "
                f"expected LCA pattern (3xO then 24xC). Got: {got_elements}"
            )
        if ligand_region == "core":
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i not in LIGAND_TAIL_ORDINALS]
        elif ligand_region == "tail":
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i in LIGAND_TAIL_ORDINALS]

    return [a.index for a in lig_heavy_atoms_all]


def get_target_residue(top, seq_id, resseq):
    matches = [r for r in top.residues if r.is_protein and r.resSeq == resseq]
    if not matches:
        raise ValueError(f"[{seq_id}] No protein residue with resSeq={resseq} found.")
    if len(matches) > 1:
        raise ValueError(
            f"[{seq_id}] Multiple protein residues with resSeq={resseq} found "
            f"(multi-chain topology?): {[(r.chain.index, r.name) for r in matches]}. "
            f"This script assumes a single-chain protein."
        )
    return matches[0]


def min_dist_occupancy(traj, atom_indices, lig_atoms):
    """Per-frame min distance from atom_indices to lig_atoms, plus the
    fraction of frames within CUTOFF_NM (occupancy)."""
    if not atom_indices:
        return None, None
    pairs = np.array([[a, l] for a in atom_indices for l in lig_atoms], dtype=int)
    dists = md.compute_distances(traj, pairs)      # (n_frames, n_pairs)
    min_dist = dists.min(axis=1)
    occupancy = float((min_dist < CUTOFF_NM).mean())
    return occupancy, float(min_dist.mean())


def main():
    args     = parse_args()
    seq_id   = args.seq_id
    resseq   = args.resseq
    start_ns = args.start_ns
    end_ns   = args.end_ns
    start_ps = int(start_ns * 1000)
    end_ps   = int(end_ns   * 1000)
    TAG      = f"{int(start_ns)}_{int(end_ns)}ns"
    ligand_region = args.ligand_region
    REGION_TAG = "" if ligand_region == "whole" else f"_{ligand_region}"

    tagged_output_dir = os.path.join(
        output_base, f"residue_atomsplit_results_{TAG}{REGION_TAG}"
    )
    os.makedirs(tagged_output_dir, exist_ok=True)

    print(f"\n=== residue_atom_split_contact.py  |  seq_id={seq_id}  resseq={resseq}  "
          f"window={start_ns:.0f}-{end_ns:.0f} ns  region={ligand_region} ===")

    summary_out = os.path.join(
        tagged_output_dir, f"{seq_id}_res{resseq}_atomsplit_{TAG}{REGION_TAG}.csv"
    )
    if os.path.exists(summary_out):
        print(f"[{seq_id}] Output already exists, skipping: {summary_out}")
        sys.exit(0)

    traj = load_trajectory(seq_id, start_ps, end_ps)
    top  = traj.topology

    lig_atoms = get_ligand_atoms(top, seq_id, ligand_region)
    print(f"[{seq_id}] Ligand '{LIGAND_RESNAME}': {len(lig_atoms)} heavy atoms "
          f"(region={ligand_region})")

    residue = get_target_residue(top, seq_id, resseq)
    print(f"[{seq_id}] Target residue: {residue.name}{residue.resSeq} "
          f"(topology index {residue.index}, chain {residue.chain.index})")

    backbone_atoms  = [a.index for a in residue.atoms if a.name in BACKBONE_ATOM_NAMES]
    sidechain_atoms = [a.index for a in residue.atoms
                        if a.element.symbol != "H" and a.name not in BACKBONE_ATOM_NAMES]

    bb_occ, bb_mindist = min_dist_occupancy(traj, backbone_atoms, lig_atoms)
    sc_occ, sc_mindist = min_dist_occupancy(traj, sidechain_atoms, lig_atoms)

    if sc_occ is None:
        print(f"[{seq_id}] WARNING: residue {residue.name}{residue.resSeq} has no "
              f"side-chain heavy atoms (Gly?) — sidechain columns will be NaN.")

    row = dict(
        seq_id               = seq_id,
        resSeq                = residue.resSeq,
        resName                = residue.name,
        ligand_region          = ligand_region,
        n_frames                = traj.n_frames,
        backbone_n_atoms        = len(backbone_atoms),
        sidechain_n_atoms       = len(sidechain_atoms),
        backbone_occupancy      = bb_occ,
        sidechain_occupancy     = sc_occ,
        backbone_mean_mindist_nm  = bb_mindist,
        sidechain_mean_mindist_nm = sc_mindist,
    )

    print(f"\n── Atom-split contact summary for {seq_id}, "
          f"{residue.name}{residue.resSeq} ─────────────────────────")
    print(f"  Backbone  : occupancy={bb_occ:.3f}  mean_min_dist={bb_mindist:.3f} nm  "
          f"(n_atoms={len(backbone_atoms)})")
    if sc_occ is not None:
        print(f"  Side chain: occupancy={sc_occ:.3f}  mean_min_dist={sc_mindist:.3f} nm  "
              f"(n_atoms={len(sidechain_atoms)})")
    else:
        print(f"  Side chain: n/a (no side-chain heavy atoms)")
    print("─" * 60)

    pd.DataFrame([row]).to_csv(summary_out, index=False)
    print(f"[{seq_id}] Wrote: {summary_out}")
    print(f"[{seq_id}] Done.")


if __name__ == "__main__":
    main()
