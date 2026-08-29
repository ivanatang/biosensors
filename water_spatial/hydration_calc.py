#!/usr/bin/env python3
"""
Computes per-frame water-oxygen counts within several first-hydration-shell
distances of a reference atom set, for a single sequence.

Reference region: `--reference-region ligand` (default) uses LCA's own heavy
atoms; `--reference-region pocket_residues` uses the 27 consensus
pocket-lining protein residues from Leonard et al. 2024 (`POCKET_RESIDS`,
kept in sync with the same list in compute_Rg_sasa.sh). Both are genuine
water-density measurements -- raw water-oxygen counts near real atoms, never
counting positions inside the reference atoms themselves, so there's no
atom-occlusion confound the way there is for a flat geometric sphere average
(see water_spatial/extract_water_spatial_feats.py's pocket_density()). Since
the reference atom set is the same molecule/residues in every sequence, raw
counts across sequences are density-comparable (same effective denominator).

Rationale
---------
LCA is strongly hydrophobic. If binders settle the ligand into the pocket in
a way that excludes water (the classic hydrophobic-effect driving force for
binding), their hydration-shell counts should trend lower than nonbinders',
whose pockets may stay more solvent-exposed. Complements R_score_calc.py's
per-residue water-mediated-CONTACT score (which only counts water bridging a
specific residue-ligand contact) with a bulk measure of how much water sits
near the reference region overall, whether or not it's bridging any
particular contact.

Cutoffs: multiple radii computed in one pass (default 3/4/5/6/8 A) from a
single per-frame minimum-distance array -- 3.5 A is the standard
first-hydration-shell O...O radius; the wider cutoffs give a coarse shell
profile. No normalization is applied: the reference region is the same
across all sequences (same ligand, or the same pocket-residue set), so raw
per-frame counts are already comparable across sequences.

Distances use `mdtraj.compute_distances(..., periodic=True)`, which applies
the minimum-image convention correctly for triclinic boxes (this system uses
a rhombic dodecahedral box) -- a prior version of this script used a raw
`np.linalg.norm` on absolute coordinates with no periodic correction, which
could silently undercount hydration for waters/ligand near a periodic
boundary. To keep memory bounded, distances are computed one reference atom
at a time (looping over the reference set, which is small: ~27 ligand heavy
atoms or ~150-200 pocket-residue heavy atoms) rather than building the full
water x reference-atom pairs array up front, which would be tens of GB at
full trajectory length -- same class of memory blowup fixed in
extract_water_spatial_feats.py's pocket_density() earlier this project.

Output
------
  {out_dir}/{seq_id}_hydration_{region}_{TAG}.csv
    columns: time_ns, hydration_count_{region}_{cutoff}A (one per cutoff)

Usage
-----
    conda activate biosensors
    python hydration_calc.py --seq_id pair_3059_binder --seq_type binders
    python hydration_calc.py --seq_id pair_3059_binder --seq_type binders \
        --reference-region pocket_residues
"""

import os
import warnings
import numpy as np
import pandas as pd
import mdtraj as md
import argparse

warnings.filterwarnings("ignore", category=DeprecationWarning)

# Kept in sync with compute_Rg_sasa.sh's POCKET_RESIDS -- the 27 consensus
# pocket-lining residue positions (Leonard et al. 2024). If that list ever
# changes, update both.
POCKET_RESIDS = {59, 60, 61, 62, 79, 81, 83, 87, 88, 89, 91, 92, 94, 108, 109,
                  110, 115, 117, 120, 122, 141, 158, 159, 160, 163, 164, 167}

parser = argparse.ArgumentParser()
parser.add_argument('--seq_id',    required=True)
parser.add_argument('--seq_type',  required=True)
parser.add_argument('--start-ns',  type=float, default=40.0,
                    help='Start of analysis window in ns (default: 40)')
parser.add_argument('--end-ns',    type=float, default=500.0,
                    help='End of analysis window in ns (default: 500)')
parser.add_argument('--reference-region', choices=['ligand', 'pocket_residues'],
                    default='ligand',
                    help="Reference atom set for distance measurement "
                         "(default: ligand)")
parser.add_argument('--cutoffs', type=str, default='3,4,5,6,8',
                    help='Comma-separated hydration-shell cutoffs in Angstrom '
                         '(default: 3,4,5,6,8)')
parser.add_argument('--stride',    type=int, default=10,
                    help='Frame stride (default: 10; set to 1 for publication quality)')
parser.add_argument('--suffix', default='',
                    help="Run-directory/output suffix, e.g. '_qfix' for the "
                         "bond-order/charge-fix systems (default: '', the "
                         "standard production directory)")
args = parser.parse_args()
seq_id   = args.seq_id
seq_type = args.seq_type

# ── Analysis window ────────────────────────────────────────────────
START_NS = args.start_ns
END_NS   = args.end_ns
START_PS = int(START_NS * 1000)
END_PS   = int(END_NS   * 1000)
TAG      = f"{int(START_NS)}_{int(END_NS)}ns"   # e.g. "40_500ns"
REGION_TAG = {'ligand': 'ligand', 'pocket_residues': 'pocket'}[args.reference_region]

CUTOFFS_ANG = [float(x) for x in args.cutoffs.split(',')]

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG  ← edit these before running
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory inputs are read from the PetaLibrary archive, not scratch --
# scratch auto-deletes after 90 days and older runs' xtc/gro may already be
# gone, and this script needs the raw solvated trajectory (water retained),
# which the trimmed/protein-only scratch intermediates used elsewhere don't
# have. Mirrors R_score_calc.py's input_base.
input_base  = "/pl/active/shirts_archive/IvanaTang/biosensors"
output_base = "/scratch/alpine/ivta1597/LCA_boltz_models"
prod   = f"prod_md_0p9_cutoff_3dt_64x1_16PME_642dd{args.suffix}"

traj_path = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(input_base, seq_type, seq_id, prod, "prod_md_500ns.gro")

out_dir = os.path.join(output_base, seq_type, seq_id, f"water_density_{TAG}{args.suffix}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "WAT", "SOL"}

STRIDE = args.stride


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD TRAJECTORY
# ─────────────────────────────────────────────────────────────────────────────
traj = md.load(traj_path, top=top_path, stride=STRIDE)

mask = (traj.time >= START_PS) & (traj.time <= END_PS)
traj = traj[mask]
print(f"  Time window: {START_NS:.0f}-{END_NS:.0f} ns  "
      f"({mask.sum()} of {mask.shape[0]} frames retained after striding)")

top = traj.topology
nf  = traj.n_frames
print(f"  {nf} frames  |  {traj.n_atoms} atoms")

# ─────────────────────────────────────────────────────────────────────────────
# 2. PARSE TOPOLOGY — reference atom set and water oxygens
# ─────────────────────────────────────────────────────────────────────────────
all_res = list(top.residues)
wat_res = [r for r in all_res if r.name in WATER_RESNAMES]

if args.reference_region == 'ligand':
    lig_res = [r for r in all_res if r.name == LIG_RESNAME]
    if not lig_res:
        raise ValueError(
            f"No residue named '{LIG_RESNAME}' found in topology.\n"
            "Check LIG_RESNAME in the CONFIG block.")
    ref_atoms = [a for r in lig_res for a in r.atoms if a.element.symbol != 'H']
else:
    pocket_res = [r for r in all_res if r.is_protein and r.resSeq in POCKET_RESIDS]
    found_resids = {r.resSeq for r in pocket_res}
    missing = POCKET_RESIDS - found_resids
    if missing:
        print(f"  WARNING: pocket resids not found in topology: {sorted(missing)}")
    if not pocket_res:
        raise ValueError(
            f"No protein residues matching POCKET_RESIDS found in topology.")
    ref_atoms = [a for r in pocket_res for a in r.atoms if a.element.symbol != 'H']

ref_idx = np.array([a.index for a in ref_atoms], dtype=int)

# Water O atoms — require exactly 2 H atoms in the residue to confirm it's a
# real water molecule, not a stray crystallographic oxygen. Mirrors
# R_score_calc.py's water-oxygen extraction.
wat_O = []
for r in wat_res:
    O = [a.index for a in r.atoms if a.element.symbol == 'O']
    H = [a.index for a in r.atoms if a.element.symbol == 'H']
    if O and len(H) == 2:
        wat_O.append(O[0])
wat_O = np.array(wat_O, dtype=int)

print(f"  Reference region    : {args.reference_region}")
print(f"  Reference atoms     : {len(ref_idx)}")
print(f"  Water oxygens       : {len(wat_O)}")

if len(wat_O) == 0:
    raise ValueError(
        f"No water oxygens found (checked WATER_RESNAMES={WATER_RESNAMES}).\n"
        "Check the water residue naming in this topology.")

# ─────────────────────────────────────────────────────────────────────────────
# 3. PER-FRAME MINIMUM WATER-TO-REFERENCE DISTANCE
# ─────────────────────────────────────────────────────────────────────────────
# One reference atom at a time, rather than the full (water x reference)
# pairs array up front -- keeps peak memory at O(frames x n_waters) instead
# of O(frames x n_waters x n_reference_atoms), which would be tens of GB at
# full trajectory length for the pocket_residues region. periodic=True
# applies the minimum-image convention correctly for this system's triclinic
# (rhombic dodecahedral) box.
min_dist_nm = np.full((nf, len(wat_O)), np.inf, dtype=np.float32)

print("Computing minimum water-to-reference distances...")
for i, r_idx in enumerate(ref_idx):
    if i % 20 == 0:
        print(f"  reference atom {i:>4d}/{len(ref_idx)}", end='\r', flush=True)
    pairs = np.column_stack([wat_O, np.full(len(wat_O), r_idx, dtype=int)])
    d = md.compute_distances(traj, pairs, periodic=True)   # (F, W)
    np.minimum(min_dist_nm, d, out=min_dist_nm)

print()

# ─────────────────────────────────────────────────────────────────────────────
# 4. HYDRATION COUNTS AT EACH CUTOFF
# ─────────────────────────────────────────────────────────────────────────────
count_cols = {}
for cutoff_ang in CUTOFFS_ANG:
    cutoff_nm = cutoff_ang / 10.0
    col = f"hydration_count_{REGION_TAG}_{f'{cutoff_ang:g}'.replace('.', 'p')}A"
    counts = (min_dist_nm < cutoff_nm).sum(axis=1)
    count_cols[col] = counts
    print(f"  {col}: mean={counts.mean():.2f}  std={counts.std():.2f}  "
          f"min={counts.min()}  max={counts.max()}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. SAVE CSV
# ─────────────────────────────────────────────────────────────────────────────
time_ns = traj.time / 1000.0   # mdtraj reports ps
df = pd.DataFrame({"time_ns": time_ns, **count_cols})

csv_path = os.path.join(out_dir, f"{seq_id}_hydration_{REGION_TAG}_{TAG}.csv")
df.to_csv(csv_path, index=False)
print(f"Hydration count table saved -> {csv_path}")
