#!/usr/bin/env python3
"""
Computes the R-score for every protein residue relative to the ligand,
following Leonard et al., ACS Chem. Biol. 2024, 19, 1757-1772.

What is computed
----------------
For each protein residue and each trajectory frame, two boolean contact
types are recorded using a 4 Å heavy-atom distance threshold:

  Direct (D)        : any heavy atom on the residue is within 4 Å of any
                      heavy atom on the ligand — the residue touches the
                      ligand without an intervening water.

  Water-mediated (W): a water oxygen is simultaneously within 4 Å of the
                      ligand AND within 4 Å of the residue — the water acts
                      as a bridge between residue and ligand.

  Total (I)         : union of D and W — was either contact present?
                      Scales R toward zero for residues with rare contact,
                      avoiding false signal from noise.

These per-frame booleans are averaged to scalar occupancies, then combined:

  R = (D - W) * I          (Leonard et al., Methods)

  R = +1 : purely direct contact
  R =  0 : equal direct/water-mediated, or no contact at all
  R = -1 : purely water-mediated contact

Residues with R < -0.7 are classified as having "dominant" water-mediated
interactions (Figure 5A/B, Figure S3).  These are the inputs to Step 2
(hbond_threshold_calibration.py).

Output
------
  {out_dir}/{seq_id}_R_scores.csv   — full per-residue table
  {out_dir}/{seq_id}_R_scores.png   — bar chart coloured by R-score

Usage
-----
    conda activate IS_env
    python compute_r_scores.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
import argparse

warnings.filterwarnings("ignore", category=DeprecationWarning)

parser = argparse.ArgumentParser()
parser.add_argument('--seq_id',    required=True)
parser.add_argument('--seq_type',  required=True)
parser.add_argument('--start-ns',  type=float, default=40.0,
                    help='Start of analysis window in ns (default: 40)')
parser.add_argument('--end-ns',    type=float, default=500.0,
                    help='End of analysis window in ns (default: 500)')
args = parser.parse_args()
seq_id   = args.seq_id
seq_type = args.seq_type

# ── Analysis window ────────────────────────────────────────────────
START_NS = args.start_ns
END_NS   = args.end_ns
START_PS = int(START_NS * 1000)
END_PS   = int(END_NS   * 1000)
TAG      = f"{int(START_NS)}_{int(END_NS)}ns"   # e.g. "40_250ns", "40_500ns"

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG  ← edit these before running
# ─────────────────────────────────────────────────────────────────────────────
base   = "/scratch/alpine/ivta1597/LCA_boltz_models"
ext    = "HMR/dodecahedron"
prod   = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path = os.path.join(base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(base, seq_type, seq_id, prod, "prod_md_500ns.gro")
out_dir   = os.path.join(base, seq_type, seq_id, f"water_contacts_{TAG}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "WAT", "SOL"}
ION_RESNAMES   = {"NA", "CL", "NA+", "CL-"}

HEAVY_CUT = 0.40   # nm (= 4.0 Å) — heavy-atom distance threshold
R_DOM     = -0.70  # residues with R below this are "dominant" water-mediated
STRIDE    = 10     # set to 1 for publication quality


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD TRAJECTORY
# ─────────────────────────────────────────────────────────────────────────────
traj = md.load(traj_path, top=top_path, stride=STRIDE)

# ── Restrict to analysis window ────────────────────────────────────
mask = (traj.time >= START_PS) & (traj.time <= END_PS)
traj = traj[mask]
print(f"  Time window: {START_NS:.0f}–{END_NS:.0f} ns  "
      f"({mask.sum()} of {mask.shape[0]} frames retained after striding)")
# ──────────────────────────────────────────────────────────────────

top  = traj.topology
nf   = traj.n_frames
print(f"  {nf} frames  |  {traj.n_atoms} atoms")

# ─────────────────────────────────────────────────────────────────────────────
# 2. PARSE TOPOLOGY — separate ligand, waters, and protein residues
# ─────────────────────────────────────────────────────────────────────────────
def _heavy(indices):
    return [i for i in indices if top.atom(i).element.symbol != 'H']

all_res  = list(top.residues)
skip     = WATER_RESNAMES | ION_RESNAMES | {LIG_RESNAME}

lig_res  = [r for r in all_res if r.name == LIG_RESNAME]
wat_res  = [r for r in all_res if r.name in WATER_RESNAMES]
prot_res = [r for r in all_res if r.name not in skip]

if not lig_res:
    raise ValueError(
        f"No residue named '{LIG_RESNAME}' found in topology.\n"
        "Check LIG_RESNAME in the CONFIG block.")

# Ligand heavy atoms
lig_heavy = np.array(
    _heavy([a.index for r in lig_res for a in r.atoms]), dtype=int)

# Water O atoms (parallel array to wat_H pairs — not needed here but kept
# for consistency with downstream scripts)
wat_O = []
for r in wat_res:
    O = [a.index for a in r.atoms if a.element.symbol == 'O']
    H = [a.index for a in r.atoms if a.element.symbol == 'H']
    if O and len(H) == 2:
        wat_O.append(O[0])
wat_O = np.array(wat_O, dtype=int)

# Protein: heavy-atom index arrays keyed by residue topology index
prot_heavy_by_res = {}
for r in prot_res:
    hv = _heavy([a.index for a in r.atoms])
    if hv:
        prot_heavy_by_res[r.index] = np.array(hv, dtype=int)

print(f"  Ligand  : {len(lig_heavy)} heavy atoms")
print(f"  Waters  : {len(wat_O)}")
print(f"  Protein : {len(prot_heavy_by_res)} residues with heavy atoms")

# ─────────────────────────────────────────────────────────────────────────────
# 3. COMPUTE PER-FRAME DIRECT AND WATER-MEDIATED CONTACT BOOLEANS
# ─────────────────────────────────────────────────────────────────────────────
# Pre-extract coordinate slices once to avoid repeated indexing in the loop
lig_xyz  = traj.xyz[:, lig_heavy, :]   # (F, L, 3)
watO_xyz = traj.xyz[:, wat_O,     :]   # (F, W, 3)

res_indices = [r.index for r in prot_res if r.index in prot_heavy_by_res]
res_heavy   = [prot_heavy_by_res[ri] for ri in res_indices]

direct_occ = {ri: np.zeros(nf, dtype=bool) for ri in res_indices}
wmed_occ   = {ri: np.zeros(nf, dtype=bool) for ri in res_indices}

print("Computing contact occupancies…")
for f in range(nf):
    if f % 200 == 0:
        print(f"  frame {f:>5d}/{nf}", end='\r', flush=True)

    lp = lig_xyz[f]    # (L, 3)
    wp = watO_xyz[f]   # (W, 3)

    # Waters within HEAVY_CUT of any ligand heavy atom — potential bridges
    dWL = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)  # (W, L)
    near_lig_wpos = wp[dWL.min(axis=1) < HEAVY_CUT]                  # (W', 3)
    has_bridging  = near_lig_wpos.shape[0] > 0

    for ri, rh in zip(res_indices, res_heavy):
        rp = traj.xyz[f, rh, :]   # (R, 3)

        # Direct: any residue heavy atom within HEAVY_CUT of any ligand heavy atom
        dRL = np.linalg.norm(rp[:, None, :] - lp[None, :, :], axis=-1)  # (R, L)
        direct_occ[ri][f] = dRL.min() < HEAVY_CUT

        # Water-mediated: a bridging water is also within HEAVY_CUT of the residue
        if has_bridging:
            dRW = np.linalg.norm(
                rp[:, None, :] - near_lig_wpos[None, :, :], axis=-1)    # (R, W')
            wmed_occ[ri][f] = dRW.min() < HEAVY_CUT


# ─────────────────────────────────────────────────────────────────────────────
# 4. COMPUTE R-SCORES
#
#   R = (D - W) * I

# All three terms are scalar occupancies (mean of per-frame booleans):
#   D   = mean(direct_occ)
#   W   = mean(wmed_occ)
#   I   = mean(direct_occ OR wmed_occ)   ← union, not sum
#
# When D+W = 0 (residue never contacts the ligand), R is undefined → NaN.
# ─────────────────────────────────────────────────────────────────────────────
records = []
for r in prot_res:
    ri = r.index
    if ri not in direct_occ:
        continue

    D = direct_occ[ri].astype(float).mean()
    W = wmed_occ[ri].astype(float).mean()
    I = (direct_occ[ri] | wmed_occ[ri]).astype(float).mean()
    
    R = (D - W) * I if (D + W) > 1e-9 else np.nan

    records.append(dict(
        res_index = ri,
        chain     = r.chain.index,
        resSeq    = r.resSeq,
        resname   = r.name,
        D         = round(D, 4),
        W         = round(W, 4),
        I         = round(I, 4),
        R         = round(R, 4) if not np.isnan(R) else np.nan,
    ))

df = pd.DataFrame(records)

# ─────────────────────────────────────────────────────────────────────────────
# 5. SAVE CSV
# ─────────────────────────────────────────────────────────────────────────────
csv_path = os.path.join(out_dir, f"{seq_id}_R_scores_{TAG}.csv")
df.to_csv(csv_path, index=False)
print(f"R-score table saved → {csv_path}")

# Print dominant residues to stdout
dom = df[df['R'] < R_DOM].dropna(subset=['R']).sort_values('R')
print(f"\nDominant water-mediated residues (R < {R_DOM}):  {len(dom)} found")
if len(dom):
    print(dom[['resname', 'resSeq', 'chain', 'D', 'W', 'I', 'R']].to_string(index=False))
else:
    print("  None — check LIG_RESNAME or consider lowering R_DOM.")

# ─────────────────────────────────────────────────────────────────────────────
# 6. PLOT — bar chart of R-scores for residues with any contact (I > 0.02)
#
# Bars are coloured on a diverging scale:
#   red   (R ≈ -1) : dominant water-mediated
#   white (R ≈  0) : mixed or negligible
#   blue  (R ≈ +1) : dominant direct
#
# The R_DOM threshold is marked with a dashed line.
# ─────────────────────────────────────────────────────────────────────────────
df_plot = df[df['I'] > 0.02].dropna(subset=['R']).sort_values('resSeq')

if df_plot.empty:
    print("No residues with I > 0.02 — skipping plot.")
else:
    # Map R ∈ [-1, +1] to a diverging colormap
    cmap   = plt.cm.RdBu
    colors = [cmap((r + 1) / 2) for r in df_plot['R']]
    labels = [f"{row['resname']}{row['resSeq']}" for _, row in df_plot.iterrows()]

    fig, ax = plt.subplots(figsize=(max(10, len(df_plot) * 0.45), 4),
                           constrained_layout=True)

    ax.bar(range(len(df_plot)), df_plot['R'], color=colors,
           edgecolor='k', linewidth=0.35)
    ax.axhline(R_DOM, color='k', linestyle='--', linewidth=1.0,
               label=f'R_DOM = {R_DOM}  (dominant water-mediated)')
    ax.axhline(0, color='grey', linewidth=0.6, alpha=0.5)

    ax.set_xticks(range(len(df_plot)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel('R-score', fontsize=10)
    ax.set_xlabel('Residue', fontsize=10)
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(f'{seq_id} — Water-mediated contact R-scores\n'
                 f'(residues with I > 0.02,  n={len(df_plot)})',
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    # Annotate dominant residues directly on the plot
    for i, (_, row) in enumerate(df_plot.iterrows()):
        if row['R'] < R_DOM:
            ax.text(i, row['R'] - 0.05, f"{row['resname']}{row['resSeq']}",
                    ha='center', va='top', fontsize=6, rotation=90,
                    color='#8B0000')

    png_path = os.path.join(out_dir, f"{seq_id}_R_scores_{TAG}.png")
    fig.savefig(png_path, dpi=300)
    plt.show()
