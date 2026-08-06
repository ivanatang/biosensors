#!/usr/bin/env python3
"""
gate_latch_water_bridge.py
----------------------------
Tests whether a water-mediated bridge connecting the gate loop (r84-90),
the latch loop (r114-118), and the ligand forms in this apo (no HAB1)
simulation. Motivated by published PYR1-HAB1-ligand structures showing a
water network linking gate, latch, ligand, and a Trp on HAB1 (see Leonard
et al., ACS Chem. Biol. 2024, the same paper R_score_calc.py's D/W/R
framework is drawn from). HAB1 is not present in these simulations, so
this can only test whether the gate-latch-ligand sub-network assembles on
its own, not whether it reproduces the full published network -- state
that limitation explicitly wherever these results are reported.

Two related but distinct bridge definitions are computed per frame, both
using the same 4 Angstrom heavy-atom distance cutoff R_score_calc.py uses
for its D/W terms:

  triple_bridge : at least one SINGLE water molecule is simultaneously
                   within the cutoff of the ligand, a gate heavy atom, AND
                   a latch heavy atom at the same time -- the strict,
                   literal "one water touches all three" reading of "the
                   network is forming."

  co_occurrence : a ligand-gate bridging water AND a ligand-latch bridging
                   water are BOTH present in the same frame, but not
                   necessarily the same water molecule -- a looser signal
                   that both sides of the network are water-engaged at
                   once, which could represent a network still assembling
                   (two separate bridging waters that haven't merged into
                   one, or never will without HAB1 present to help
                   organize them).

The simpler two-way ligand-gate and ligand-latch bridge occupancies (the
same logic as R_score_calc.py's W term, generalized from one residue to
the pooled gate/latch region) are also reported, so triple_bridge and
co_occurrence can be read against their baseline components rather than
in isolation.

Output
------
  {out_dir}/{seq_id}_gate_latch_bridge_{TAG}{REGION_TAG}.csv  -- one-row summary

Usage
-----
    conda activate biosensors
    python gate_latch_water_bridge.py --seq_id pair_3059_binder --seq_type binders
    python gate_latch_water_bridge.py --seq_id pair_3059_binder --seq_type binders \
        --start-ns 40 --end-ns 500 --ligand-region core
"""

import os
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

parser = argparse.ArgumentParser()
parser.add_argument('--seq_id',   required=True)
parser.add_argument('--seq_type', required=True,
                    help="Directory name: binders | nonbinders | neg_low_pkt | neg_fail_gate")
parser.add_argument('--start-ns', type=float, default=40.0)
parser.add_argument('--end-ns',   type=float, default=500.0)
parser.add_argument('--ligand-region', choices=['whole', 'core', 'tail'], default='core',
                    help="Ligand atoms to test bridging against (default: core, since "
                         "gate/latch close over the ligand's steroid-ring end, not the "
                         "solvent-facing carboxylate tail -- see analysis/core_vs_tail).")
args = parser.parse_args()
seq_id, seq_type = args.seq_id, args.seq_type
ligand_region = args.ligand_region

START_NS, END_NS = args.start_ns, args.end_ns
START_PS, END_PS = int(START_NS * 1000), int(END_NS * 1000)
TAG = f"{int(START_NS)}_{int(END_NS)}ns"
REGION_TAG = "" if ligand_region == "whole" else f"_{ligand_region}"

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG -- same paths/conventions as R_score_calc.py
# ─────────────────────────────────────────────────────────────────────────────
input_base  = "/pl/active/shirts_archive/IvanaTang/biosensors"
output_base = "/scratch/alpine/ivta1597/LCA_boltz_models"
prod = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.gro")

out_dir = os.path.join(output_base, seq_type, seq_id, f"gate_latch_water_bridge_{TAG}{REGION_TAG}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "WAT", "SOL"}

GATE_RESIDUES  = set(range(84, 91))    # CLAUDE.md: Gate loop, residues 84-90
LATCH_RESIDUES = set(range(114, 119))  # CLAUDE.md: Latch, residues 114-118

# Same core/tail ligand-atom split as R_score_calc.py / contact_type_analysis.py.
LIGAND_TAIL_ORDINALS     = {0, 1, 3, 4, 7, 9, 19}   # O41,O42,C44,C45,C48,C50,C60
LIGAND_EXPECTED_ELEMENTS = ['O', 'O', 'O'] + ['C'] * 24

HEAVY_CUT = 0.40   # nm (= 4.0 A) -- same cutoff as R_score_calc.py's D/W terms
STRIDE    = 10     # set to 1 for publication quality


def main():
    summary_out = os.path.join(
        out_dir, f"{seq_id}_gate_latch_bridge_{TAG}{REGION_TAG}.csv"
    )
    if os.path.exists(summary_out):
        print(f"[{seq_id}] Output already exists, skipping: {summary_out}")
        return

    for path in (traj_path, top_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[{seq_id}] Expected file not found: {path}")

    # ── 1. Load trajectory ──────────────────────────────────────────────
    traj = md.load(traj_path, top=top_path, stride=STRIDE)
    mask = (traj.time >= START_PS) & (traj.time <= END_PS)
    traj = traj[mask]
    top = traj.topology
    nf = traj.n_frames
    print(f"[{seq_id}] {nf} frames  |  {traj.n_atoms} atoms  |  "
          f"window {START_NS:.0f}-{END_NS:.0f} ns")

    # ── 2. Parse topology ───────────────────────────────────────────────
    def _heavy(atoms):
        return [a.index for a in atoms if a.element.symbol != 'H']

    lig_res = [r for r in top.residues if r.name == LIG_RESNAME]
    if not lig_res:
        raise ValueError(f"[{seq_id}] No residue named '{LIG_RESNAME}' found.")
    lig_heavy_atoms_all = [a for r in lig_res for a in r.atoms if a.element.symbol != 'H']

    if ligand_region != "whole":
        got = [a.element.symbol for a in lig_heavy_atoms_all]
        if got != LIGAND_EXPECTED_ELEMENTS:
            raise ValueError(
                f"[{seq_id}] Ligand heavy-atom element order doesn't match the "
                f"expected LCA pattern (3xO then 24xC). Got: {got}"
            )
        if ligand_region == "core":
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i not in LIGAND_TAIL_ORDINALS]
        else:
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i in LIGAND_TAIL_ORDINALS]
    lig_heavy = np.array([a.index for a in lig_heavy_atoms_all], dtype=int)

    wat_O = []
    for r in top.residues:
        if r.name in WATER_RESNAMES:
            O = [a.index for a in r.atoms if a.element.symbol == 'O']
            H = [a.index for a in r.atoms if a.element.symbol == 'H']
            if O and len(H) == 2:
                wat_O.append(O[0])
    wat_O = np.array(wat_O, dtype=int)

    gate_heavy = np.array(_heavy(
        [a for r in top.residues if r.is_protein and r.resSeq in GATE_RESIDUES
         for a in r.atoms]
    ), dtype=int)
    latch_heavy = np.array(_heavy(
        [a for r in top.residues if r.is_protein and r.resSeq in LATCH_RESIDUES
         for a in r.atoms]
    ), dtype=int)

    if gate_heavy.size == 0:
        raise ValueError(f"[{seq_id}] No gate heavy atoms found (resSeq 84-90) -- "
                          f"check topology residue numbering.")
    if latch_heavy.size == 0:
        raise ValueError(f"[{seq_id}] No latch heavy atoms found (resSeq 114-118) -- "
                          f"check topology residue numbering.")

    print(f"[{seq_id}] Ligand: {len(lig_heavy)} heavy atoms (region={ligand_region})  |  "
          f"Waters: {len(wat_O)}  |  Gate: {len(gate_heavy)} heavy atoms  |  "
          f"Latch: {len(latch_heavy)} heavy atoms")

    # ── 3. Per-frame bridge detection ───────────────────────────────────
    lig_xyz   = traj.xyz[:, lig_heavy, :]
    watO_xyz  = traj.xyz[:, wat_O, :]
    gate_xyz  = traj.xyz[:, gate_heavy, :]
    latch_xyz = traj.xyz[:, latch_heavy, :]

    gate_bridge      = np.zeros(nf, dtype=bool)
    latch_bridge     = np.zeros(nf, dtype=bool)
    triple_bridge    = np.zeros(nf, dtype=bool)
    n_triple_waters  = np.zeros(nf, dtype=int)

    print("Computing per-frame water bridges...")
    for f in range(nf):
        if f % 200 == 0:
            print(f"  frame {f:>5d}/{nf}", end='\r', flush=True)

        wp = watO_xyz[f]
        lp = lig_xyz[f]
        gp = gate_xyz[f]
        tp = latch_xyz[f]

        near_lig   = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1).min(axis=1) < HEAVY_CUT
        near_gate  = np.linalg.norm(wp[:, None, :] - gp[None, :, :], axis=-1).min(axis=1) < HEAVY_CUT
        near_latch = np.linalg.norm(wp[:, None, :] - tp[None, :, :], axis=-1).min(axis=1) < HEAVY_CUT

        gl  = near_lig & near_gate      # waters bridging ligand<->gate
        ll  = near_lig & near_latch     # waters bridging ligand<->latch
        glt = gl & near_latch           # waters bridging ligand<->gate<->latch, all at once

        gate_bridge[f]     = gl.any()
        latch_bridge[f]    = ll.any()
        triple_bridge[f]   = glt.any()
        n_triple_waters[f] = int(glt.sum())

    print()

    # ── 4. Aggregate to occupancies ─────────────────────────────────────
    row = dict(
        seq_id                      = seq_id,
        ligand_region                = ligand_region,
        n_frames                      = nf,
        gate_bridge_occupancy          = round(float(gate_bridge.mean()), 4),
        latch_bridge_occupancy         = round(float(latch_bridge.mean()), 4),
        co_occurrence_occupancy        = round(float((gate_bridge & latch_bridge).mean()), 4),
        triple_bridge_occupancy        = round(float(triple_bridge.mean()), 4),
        mean_n_triple_bridge_waters    = round(float(n_triple_waters.mean()), 4),
    )

    print(f"\n-- Gate-latch-ligand water bridge summary for {seq_id} --------------")
    for k, v in row.items():
        print(f"  {k}: {v}")
    print("-" * 60)

    pd.DataFrame([row]).to_csv(summary_out, index=False)
    print(f"[{seq_id}] Wrote: {summary_out}")
    print(f"[{seq_id}] Done.")


if __name__ == "__main__":
    main()
