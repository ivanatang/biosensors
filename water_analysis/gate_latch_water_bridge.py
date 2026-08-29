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

PERIODIC-BOUNDARY FIX (important, read this): an earlier version of this
script computed distances with plain `np.linalg.norm` on raw coordinates,
with no minimum-image correction. Verified directly on pair_3085_binder
(40-500ns window): that method gives triple_bridge_occupancy=0.7832
(exactly what's in the already-recorded lab notebook entry and the
water_bridge ML feature group); the periodic-aware method below gives
0.8704. 91.3% frame-by-frame agreement between the two methods, and of
the disagreements, ALL are false negatives in the old method (107 frames
where a real bridge existed but raw-coordinate distance missed it because
the water was stored in a different periodic image) and there are ZERO
false positives. This is a one-directional systematic UNDERCOUNT, not
random noise -- it would affect every sequence in the cohort the same
way. All distance calculations here now use mdtraj's compute_distances
with periodic=True, which handles this repo's triclinic (dodecahedron)
boxes correctly. Existing occupancy numbers computed before this fix
should be treated as underestimates until recomputed; the direction of
any Binder-vs-FP comparison is plausibly still valid since both groups
would be similarly affected, but that has NOT been separately confirmed.

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
                   necessarily the same water molecule.

DYNAMICS METRICS (new): beyond the aggregate occupancy, this now also
tracks, per sequence:

  first_appearance_ns : simulation time of the first frame where a triple
                         bridge is present (NaN if it never appears).

  n_runs               : number of separate contiguous stretches where a
                         triple bridge is continuously present (i.e. how
                         many times it forms/breaks over the trajectory --
                         "is it always there, or does it flicker").

  n_distinct_waters     : number of DIFFERENT water molecules that ever
                         fill the bridging role over the trajectory (asks
                         "is it the same water the whole time, or does the
                         specific molecule doing the bridging keep
                         changing").

  mean_run_duration_ns,
  median_run_duration_ns : how long a given water molecule typically holds
                         the bridging position before being replaced,
                         estimated from run lengths (a lower bound on true
                         residence time -- limited by STRIDE's temporal
                         resolution, see note below).

Because first_appearance_ns is only meaningful if the window actually
starts at t=0, START_NS defaults to 0.0 here (NOT the 40 ns equilibration
cutoff used by convention elsewhere in this repo) -- override with
--start-ns if you specifically want the old 40-500ns occupancy-only
comparison window.

Output
------
  {out_dir}/{seq_id}_gate_latch_bridge_{TAG}{REGION_TAG}.csv  -- one-row summary

Usage
-----
    conda activate biosensors
    python gate_latch_water_bridge.py --seq_id pair_3059_binder --seq_type binders
    python gate_latch_water_bridge.py --seq_id pair_3059_binder --seq_type binders \
        --start-ns 0 --end-ns 500 --ligand-region core
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
parser.add_argument('--start-ns', type=float, default=0.0,
                    help="Default 0 (not the usual 40 ns cutoff) so first_appearance_ns "
                         "is measured from the true start of the production run.")
parser.add_argument('--end-ns',   type=float, default=500.0)
parser.add_argument('--ligand-region', choices=['whole', 'core', 'tail'], default='core',
                    help="Ligand atoms to test bridging against (default: core, since "
                         "gate/latch close over the ligand's steroid-ring end, not the "
                         "solvent-facing carboxylate tail -- see analysis/core_vs_tail).")
parser.add_argument('--suffix', default='',
                    help="Run-directory/output suffix, e.g. '_qfix' for the "
                         "bond-order/charge-fix systems (default: '', the "
                         "standard production directory)")
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
prod = f"prod_md_0p9_cutoff_3dt_64x1_16PME_642dd{args.suffix}"

traj_path = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.gro")

out_dir = os.path.join(output_base, seq_type, seq_id, f"gate_latch_water_bridge_{TAG}{REGION_TAG}{args.suffix}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "WAT", "SOL"}

GATE_RESIDUES  = set(range(84, 91))    # CLAUDE.md: Gate loop, residues 84-90
LATCH_RESIDUES = set(range(114, 119))  # CLAUDE.md: Latch, residues 114-118

# Same core/tail ligand-atom split as R_score_calc.py / contact_type_analysis.py.
LIGAND_TAIL_ORDINALS     = {0, 1, 3, 4, 7, 9, 19}   # O41,O42,C44,C45,C48,C50,C60
LIGAND_EXPECTED_ELEMENTS = ['O', 'O', 'O'] + ['C'] * 24

HEAVY_CUT = 0.40   # nm (= 4.0 A) -- same cutoff as R_score_calc.py's D/W terms
STRIDE    = 10     # 37.5 ps * 10 = 375 ps spacing. Run-duration metrics are a
                    # LOWER BOUND on true water residence time at this
                    # resolution -- individual waters may exchange faster
                    # than this stride can resolve. Set to 1 for a precise
                    # answer on residence time specifically (much slower).


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

    wat_res_list = [r for r in top.residues if r.name in WATER_RESNAMES]
    wat_O = []
    wat_res_of_O = []
    for r in wat_res_list:
        O = [a.index for a in r.atoms if a.element.symbol == 'O']
        H = [a.index for a in r.atoms if a.element.symbol == 'H']
        if O and len(H) == 2:
            wat_O.append(O[0])
            wat_res_of_O.append(r)
    wat_O = np.array(wat_O, dtype=int)
    wat_resSeq = np.array([r.resSeq for r in wat_res_of_O], dtype=int)

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

    # ── 3. Per-frame bridge detection, periodic-aware, staged for speed ──
    # Stage 1: which waters are EVER near the ligand core at all (periodic
    # minimum-image distance) -- cuts the ~9000-water pool down to only the
    # ones worth checking against gate/latch in stage 2, since computing
    # gate/latch distances for every water in every frame is far more
    # expensive than needed (most waters never come near the pocket at all).
    lig_pairs = np.array([[w, l] for w in wat_O for l in lig_heavy], dtype=int)
    d_lig = md.compute_distances(traj, lig_pairs, periodic=True).reshape(
        nf, len(wat_O), len(lig_heavy))
    near_lig_all = d_lig.min(axis=2) < HEAVY_CUT
    cand_idx = np.where(near_lig_all.any(axis=0))[0]
    print(f"[{seq_id}] candidate waters ever near ligand core: {len(cand_idx)} of {len(wat_O)}")

    if len(cand_idx) == 0:
        near_lig_cand = np.zeros((nf, 0), dtype=bool)
        near_gate_cand = np.zeros((nf, 0), dtype=bool)
        near_latch_cand = np.zeros((nf, 0), dtype=bool)
        cand_resSeq = np.array([], dtype=int)
    else:
        wat_O_cand = wat_O[cand_idx]
        cand_resSeq = wat_resSeq[cand_idx]
        gate_pairs = np.array([[w, g] for w in wat_O_cand for g in gate_heavy], dtype=int)
        latch_pairs = np.array([[w, l] for w in wat_O_cand for l in latch_heavy], dtype=int)
        d_gate = md.compute_distances(traj, gate_pairs, periodic=True).reshape(
            nf, len(wat_O_cand), len(gate_heavy))
        d_latch = md.compute_distances(traj, latch_pairs, periodic=True).reshape(
            nf, len(wat_O_cand), len(latch_heavy))
        near_gate_cand = d_gate.min(axis=2) < HEAVY_CUT
        near_latch_cand = d_latch.min(axis=2) < HEAVY_CUT
        near_lig_cand = near_lig_all[:, cand_idx]

    gate_bridge_frame  = near_lig_cand.any(axis=1) if near_lig_cand.shape[1] else np.zeros(nf, bool)
    # Recompute gate/latch bridge (ligand<->gate, ligand<->latch) using the
    # SAME candidate-water restriction for consistency across all metrics.
    gl = near_lig_cand & near_gate_cand
    ll = near_lig_cand & near_latch_cand
    glt = gl & near_latch_cand   # (nf, n_cand) -- triple bridge per candidate water

    gate_bridge  = gl.any(axis=1)
    latch_bridge = ll.any(axis=1)
    triple_bridge = glt.any(axis=1)
    n_triple_waters = glt.sum(axis=1)

    # ── 4. Dynamics metrics ─────────────────────────────────────────────
    if triple_bridge.any():
        first_idx = int(np.argmax(triple_bridge))
        first_appearance_ns = float(traj.time[first_idx]) / 1000.0
    else:
        first_appearance_ns = float("nan")

    # identity per frame: which candidate water(s) satisfy the triple
    # condition (resSeq), used for run-counting and distinct-water-counting
    frame_waters = []
    for f in range(nf):
        qualifying = np.where(glt[f])[0]
        frame_waters.append(tuple(sorted(int(cand_resSeq[i]) for i in qualifying))
                             if len(qualifying) else None)

    distinct_waters = set()
    for fw in frame_waters:
        if fw is not None:
            distinct_waters.update(fw)
    n_distinct_waters = len(distinct_waters)

    # contiguous runs of "a bridge is present" (regardless of which water),
    # to count how many separate formation events there are
    runs = []
    prev_present = False
    run_len = 0
    for present in triple_bridge:
        if present and not prev_present:
            run_len = 1
        elif present and prev_present:
            run_len += 1
        elif not present and prev_present:
            runs.append(run_len)
        prev_present = bool(present)
    if prev_present:
        runs.append(run_len)

    n_runs = len(runs)
    time_per_frame_ns = STRIDE * 0.0375  # 37.5 ps raw spacing * stride
    if runs:
        mean_run_duration_ns = float(np.mean(runs)) * time_per_frame_ns
        median_run_duration_ns = float(np.median(runs)) * time_per_frame_ns
    else:
        mean_run_duration_ns = float("nan")
        median_run_duration_ns = float("nan")

    # ── 5. Aggregate to occupancies + write ─────────────────────────────
    row = dict(
        seq_id                      = seq_id,
        ligand_region                = ligand_region,
        n_frames                      = nf,
        window_start_ns                = START_NS,
        window_end_ns                  = END_NS,
        gate_bridge_occupancy          = round(float(gate_bridge.mean()), 4),
        latch_bridge_occupancy         = round(float(latch_bridge.mean()), 4),
        co_occurrence_occupancy        = round(float((gate_bridge & latch_bridge).mean()), 4),
        triple_bridge_occupancy        = round(float(triple_bridge.mean()), 4),
        mean_n_triple_bridge_waters    = round(float(n_triple_waters.mean()), 4),
        first_appearance_ns            = round(first_appearance_ns, 4) if first_appearance_ns == first_appearance_ns else float("nan"),
        n_runs                          = n_runs,
        n_distinct_waters               = n_distinct_waters,
        mean_run_duration_ns            = round(mean_run_duration_ns, 4) if mean_run_duration_ns == mean_run_duration_ns else float("nan"),
        median_run_duration_ns          = round(median_run_duration_ns, 4) if median_run_duration_ns == median_run_duration_ns else float("nan"),
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
