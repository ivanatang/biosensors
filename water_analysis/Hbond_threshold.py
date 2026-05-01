"""
H-bond threshold calibration — Figure S4 reproduction
======================================================
Reproduces the control-residue analysis from Leonard et al. (ACS Chem. Biol.
2024) Supplementary Figure S4, used to validate and derive the thresholds:

    "Strong" : mean H–acceptor distance < STRONG_MAX  Å
    "Stable" : std  H–acceptor distance < STABLE_STD  Å

For each dominant residue (R < -0.7 from water_hbond_stability.py), two
control residues of the **same amino acid type** are identified:

  Solvent-exposed control
    — Same AA type as reference residue
    — High relative SASA (> SASA_EXPOSED_FRAC of the max observed for that type)
    — Zero ligand contact occupancy (I == 0 in R-score table)

  Pocket control
    — Same AA type as reference residue
    — Located near the binding site (within POCKET_RADIUS of ligand centroid)
    — Very low ligand contact occupancy (I < POCKET_I_MAX, i.e. <10%)

Water–H-acceptor distances are computed for each control residue across the
trajectory, then KDE curves are plotted per reference residue (mirroring
Figure S4 panels).  Thresholds are derived from the combined distribution of
all control bonds:

    STRONG_MAX = percentile at COVERAGE_PCT of combined per-residue means
    STABLE_STD = percentile at COVERAGE_PCT of combined per-residue std devs

Usage
-----
1. Run water_hbond_stability.py first to produce {seq_id}_R_scores.csv.
2. Adjust CONFIG below (paths, seq_id, LIG_RESNAME).
3. Run:  conda activate IS_env && python hbond_threshold_calibration.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG  (match settings in water_hbond_stability.py)
# ─────────────────────────────────────────────────────────────────────────────
base   = "/scratch/alpine/ivta1597/LCA_boltz_models"
seq_type = "binders"
seq_id = "pair_XXXX"
prod   = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path    = os.path.join(base, seq_id, prod, "prod_md_500ns.xtc")
top_path     = os.path.join(base, seq_id, prod, "prod_md_500ns.gro")
rscores_csv  = os.path.join(base, seq_id, "water_contacts", f"{seq_id}_R_scores.csv")
out_dir      = os.path.join(base, seq_id, "water_contacts")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"
WATER_RESNAMES = {"HOH", "SOL"}
ION_RESNAMES   = {"NA", "CL", "NA+", "CL-"}

R_DOM           = -0.70   # dominant water-mediated threshold
HEAVY_CUT       = 0.40    # nm  (= 4.0 Å)
POCKET_RADIUS   = 1.00    # nm  (= 10 Å)  — residues within this of ligand centroid
POCKET_I_MAX    = 0.10    # max combined ligand contact occupancy for pocket control
SASA_EXPOSED_FRAC = 0.40  # min fraction of max-per-type SASA for exposed control

COVERAGE_PCT = 95.0       # percentile used to set thresholds
STRIDE       = 10         # trajectory stride for H-bond distance computation
CALIB_STRIDE = 10         # stride for SASA averaging


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def heavy(indices, top):
    return [i for i in indices if top.atom(i).element.symbol != 'H']

def acceptors(indices, top):
    return [i for i in indices if top.atom(i).element.symbol in ('O', 'N')]


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD
# ─────────────────────────────────────────────────────────────────────────────
print(f"Loading trajectory (stride={STRIDE})…")
traj = md.load(traj_path, top=top_path, stride=STRIDE)
top  = traj.topology
nf   = traj.n_frames
print(f"  {nf} frames | {traj.n_atoms} atoms")

df_R = pd.read_csv(rscores_csv)
print(f"  R-score table: {len(df_R)} residues")


# ─────────────────────────────────────────────────────────────────────────────
# 2. PARSE TOPOLOGY
# ─────────────────────────────────────────────────────────────────────────────
all_res     = list(top.residues)
lig_res     = [r for r in all_res if r.name == LIG_RESNAME]
wat_res     = [r for r in all_res if r.name in WATER_RESNAMES]
skip        = WATER_RESNAMES | ION_RESNAMES | {LIG_RESNAME}
prot_res    = [r for r in all_res if r.name not in skip]

if not lig_res:
    raise ValueError(f"No residue named '{LIG_RESNAME}' in topology.")

lig_all_idx = [a.index for r in lig_res for a in r.atoms]
lig_heavy   = np.array(heavy(lig_all_idx, top), dtype=int)

wat_O_idx, wat_H_idx = [], []
for r in wat_res:
    O = [a.index for a in r.atoms if a.element.symbol == 'O']
    H = [a.index for a in r.atoms if a.element.symbol == 'H']
    if O and len(H) == 2:
        wat_O_idx.append(O[0])
        wat_H_idx.append((H[0], H[1]))
wat_O_idx = np.array(wat_O_idx, dtype=int)
H_by_O    = {wat_O_idx[i]: wat_H_idx[i] for i in range(len(wat_O_idx))}

prot_heavy_by_res = {}
prot_acc_by_res   = {}
for r in prot_res:
    all_idx = [a.index for a in r.atoms]
    hv  = heavy(all_idx, top)
    acc = acceptors(all_idx, top)
    if hv:
        prot_heavy_by_res[r.index] = np.array(hv,  dtype=int)
        prot_acc_by_res[r.index]   = np.array(acc, dtype=int)

res_by_index = {r.index: r for r in prot_res}
print(f"  Ligand heavy: {len(lig_heavy)} | Waters: {len(wat_O_idx)} "
      f"| Protein residues: {len(prot_res)}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. DOMINANT RESIDUES  (from R-score CSV)
# ─────────────────────────────────────────────────────────────────────────────
dom = df_R[df_R['R'] < R_DOM].dropna(subset=['R']).copy()
print(f"\nDominant residues (R < {R_DOM}):")
print(dom[['resname', 'resSeq', 'R', 'I']].to_string(index=False))


# ─────────────────────────────────────────────────────────────────────────────
# 4. COMPUTE MEAN SASA PER RESIDUE  (averaged over sub-sampled frames)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\nComputing SASA (calib_stride={CALIB_STRIDE})…")
sub_sasa = traj[::max(1, CALIB_STRIDE // STRIDE)]  # relative to loaded stride
sasa_arr = md.shrake_rupley(sub_sasa, mode='residue')   # (F, n_res)
mean_sasa = sasa_arr.mean(axis=0)   # (n_res,)  in nm²

# Topology residue order matches sasa_arr columns
topo_res_list = list(top.residues)
sasa_by_res_idx = {topo_res_list[i].index: mean_sasa[i]
                   for i in range(len(topo_res_list))}


# ─────────────────────────────────────────────────────────────────────────────
# 5. LIGAND CENTROID  (mean position across frames)
# ─────────────────────────────────────────────────────────────────────────────
lig_centroid = traj.xyz[:, lig_heavy, :].mean(axis=(0, 1))   # (3,) nm


# ─────────────────────────────────────────────────────────────────────────────
# 6. FIND CONTROL RESIDUES
# ─────────────────────────────────────────────────────────────────────────────
def residue_centroid(res, traj):
    """Mean position of heavy atoms for this residue, averaged over frames."""
    hv = heavy([a.index for a in res.atoms], top)
    if not hv:
        return None
    return traj.xyz[:, hv, :].mean(axis=(0, 1))   # (3,) nm


def find_controls(ref_resname, ref_res_idx, df_R,
                  prot_res, sasa_by_res_idx,
                  lig_centroid, traj):
    """
    For a given reference residue type, find:
      exposed : same AA, high SASA, I == 0 (no ligand contact), not the ref
      pocket  : same AA, near ligand centroid, I < POCKET_I_MAX, not the ref

    Returns
    -------
    dict with keys 'exposed' and 'pocket', each a list of residue objects
    (may be empty if no suitable candidate found).
    """
    # Build lookup: res_index → I occupancy
    I_by_res = dict(zip(df_R['res_index'], df_R['I'].fillna(0.0)))

    # Max SASA for this residue type (to compute relative SASA)
    type_sasa = [sasa_by_res_idx.get(r.index, 0.0)
                 for r in prot_res if r.name == ref_resname]
    if not type_sasa:
        return {'exposed': [], 'pocket': []}
    max_sasa = max(type_sasa) if max(type_sasa) > 0 else 1.0

    exposed_cands, pocket_cands = [], []

    for r in prot_res:
        if r.name != ref_resname or r.index == ref_res_idx:
            continue

        I    = I_by_res.get(r.index, 0.0)
        sasa = sasa_by_res_idx.get(r.index, 0.0)
        rel_sasa = sasa / max_sasa

        # Distance from residue centroid to ligand centroid
        cen = residue_centroid(r, traj)
        if cen is None:
            continue
        dist_to_lig = np.linalg.norm(cen - lig_centroid)

        # Solvent-exposed: high SASA, no ligand contact
        if rel_sasa >= SASA_EXPOSED_FRAC and I == 0.0:
            exposed_cands.append((rel_sasa, r))

        # Pocket: near ligand, low contact
        if dist_to_lig <= POCKET_RADIUS and I < POCKET_I_MAX and I > 0:
            pocket_cands.append((dist_to_lig, r))

    # Pick best candidate for each type
    exposed = [r for _, r in sorted(exposed_cands, reverse=True)[:3]]
    pocket  = [r for _, r in sorted(pocket_cands)[:3]]

    return {'exposed': exposed, 'pocket': pocket}


print("\nFinding control residues…")
controls = {}
for _, row in dom.iterrows():
    ri    = int(row['res_index'])
    rname = row['resname']
    rseq  = row['resSeq']
    label = f"{rname}{rseq}"

    ctrls = find_controls(rname, ri, df_R, prot_res,
                          sasa_by_res_idx, lig_centroid, traj)
    controls[label] = ctrls

    exp_labels = [f"{r.name}{r.resSeq}" for r in ctrls['exposed']]
    pkt_labels = [f"{r.name}{r.resSeq}" for r in ctrls['pocket']]
    print(f"  Ref {label:10s} → exposed: {exp_labels}  pocket: {pkt_labels}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. COMPUTE H-ACCEPTOR DISTANCES FOR CONTROL RESIDUES
# ─────────────────────────────────────────────────────────────────────────────
watO_xyz = traj.xyz[:, wat_O_idx, :]   # (F, W, 3)


def ha_distances_for_residue(res, traj, wat_O_idx, watO_xyz, H_by_O, HEAVY_CUT):
    """
    For a single residue, collect H-acceptor distances for all frames where
    any water is within HEAVY_CUT of the residue.  This captures typical
    water-protein H-bonds at that residue, independent of ligand.

    Returns np.ndarray of distances in Å.
    """
    rh_idx = prot_heavy_by_res.get(res.index, np.array([], dtype=int))
    ra_idx = prot_acc_by_res.get(res.index,   np.array([], dtype=int))
    if len(rh_idx) == 0 or len(ra_idx) == 0:
        return np.array([])

    ha_dists = []

    for f in range(traj.n_frames):
        rp = traj.xyz[f, rh_idx, :]   # (R, 3)
        wp = watO_xyz[f]              # (W, 3)

        # Waters within HEAVY_CUT of any residue heavy atom
        dWR = np.linalg.norm(wp[:, None, :] - rp[None, :, :], axis=-1)
        near_mask = dWR.min(axis=1) < HEAVY_CUT
        if not near_mask.any():
            continue

        near_O_global = wat_O_idx[near_mask]

        min_d = np.inf
        for O_idx in near_O_global:
            H_pair = H_by_O.get(O_idx)
            if H_pair is None:
                continue
            for H_idx in H_pair:
                H_pos = traj.xyz[f, H_idx, :]
                for acc_idx in ra_idx:
                    d = np.linalg.norm(H_pos - traj.xyz[f, acc_idx, :])
                    if d < min_d:
                        min_d = d

        if min_d < np.inf:
            ha_dists.append(min_d * 10.0)   # nm → Å

    return np.array(ha_dists)


print("\nComputing H-acceptor distances for control residues…")
ctrl_data = {}   # label → {'exposed': [...arrays], 'pocket': [...arrays]}

for ref_label, ctrls in controls.items():
    ctrl_data[ref_label] = {'exposed': [], 'pocket': []}

    for kind in ('exposed', 'pocket'):
        for res in ctrls[kind]:
            print(f"  {ref_label} [{kind}] {res.name}{res.resSeq}…", end=' ')
            arr = ha_distances_for_residue(res, traj, wat_O_idx, watO_xyz,
                                           H_by_O, HEAVY_CUT)
            print(f"{len(arr)} contact frames")
            ctrl_data[ref_label][kind].append(arr)


# ─────────────────────────────────────────────────────────────────────────────
# 8. DERIVE THRESHOLDS FROM COMBINED CONTROL DISTRIBUTIONS
# ─────────────────────────────────────────────────────────────────────────────
all_means, all_stds = [], []

for ref_label, kinds in ctrl_data.items():
    for kind in ('exposed', 'pocket'):
        for arr in kinds[kind]:
            if len(arr) >= 5:
                all_means.append(arr.mean())
                all_stds.append(arr.std())

if len(all_means) < 3:
    print("\nWARNING: Too few control residues found. "
          "Try lowering SASA_EXPOSED_FRAC or POCKET_I_MAX.")
    STRONG_MAX = 2.5
    STABLE_STD = 0.45
else:
    all_means = np.array(all_means)
    all_stds  = np.array(all_stds)
    STRONG_MAX = float(np.percentile(all_means, COVERAGE_PCT))
    STABLE_STD = float(np.percentile(all_stds,  COVERAGE_PCT))

    # Fraction of controls within thresholds (validation check)
    frac_strong = (all_means < STRONG_MAX).mean() * 100
    frac_stable = (all_stds  < STABLE_STD).mean() * 100

    print(f"\n{'─'*55}")
    print(f"  Derived thresholds  ({COVERAGE_PCT}th percentile of {len(all_means)} controls):")
    print(f"    STRONG_MAX = {STRONG_MAX:.3f} Å  ({frac_strong:.0f}% of controls have mean < this)")
    print(f"    STABLE_STD = {STABLE_STD:.3f} Å  ({frac_stable:.0f}% of controls have std  < this)")
    print(f"{'─'*55}")


# ─────────────────────────────────────────────────────────────────────────────
# 9. PLOT — KDE distributions per reference residue (Figure S4 style)
# ─────────────────────────────────────────────────────────────────────────────
COLOR_EXP = '#5B9BD5'   # light blue  (solvent exposed)
COLOR_PKT = '#1F4E79'   # dark navy   (pocket)

n_panels = len(ctrl_data)
if n_panels == 0:
    print("No dominant residues found — nothing to plot.")
else:
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(2.8 * n_panels, 3.8),
                             constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    x_grid = np.linspace(0, 5, 500)

    for ax, (ref_label, kinds) in zip(axes, ctrl_data.items()):

        plotted_any = False

        for kind, color, zorder in [('exposed', COLOR_EXP, 2),
                                     ('pocket',  COLOR_PKT, 3)]:
            arrays = [a for a in kinds[kind] if len(a) >= 10]
            if not arrays:
                continue

            # Pool all distances from all controls of this kind
            pooled = np.concatenate(arrays)
            try:
                kde = gaussian_kde(pooled, bw_method='scott')
                y   = kde(x_grid)
                ax.fill_between(x_grid, y, alpha=0.35, color=color, zorder=zorder)
                ax.plot(x_grid, y, color=color, linewidth=1.5,
                        label='Solvent Exposed' if kind == 'exposed' else 'Pocket',
                        zorder=zorder)
                plotted_any = True
            except Exception:
                pass

        # Threshold line
        ax.axvline(STRONG_MAX, color='#d62728', linestyle='--',
                   linewidth=1.0, alpha=0.8,
                   label=f'STRONG_MAX = {STRONG_MAX:.2f} Å')

        ax.set_title(f'Ref = {ref_label}', fontsize=8, fontweight='bold')
        ax.set_xlabel('H–Acceptor Distance (Å)', fontsize=7)
        ax.set_ylabel('Kernel Density', fontsize=7)
        ax.set_xlim(0, 5)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)
        if ax is axes[0] and plotted_any:
            ax.legend(fontsize=6, loc='upper right')

    fig.suptitle(
        f'{seq_id} — Control H-bond distributions (Figure S4 style)\n'
        f'STRONG_MAX = {STRONG_MAX:.2f} Å  |  STABLE_STD = {STABLE_STD:.2f} Å  '
        f'({COVERAGE_PCT}th pctile, n={len(all_means)} controls)',
        fontsize=9)

    out = os.path.join(out_dir, f"{seq_id}_hbond_calibration_S4.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"\nFigure saved → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 10. SUPPLEMENTARY: mean / std summary for each control residue
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print(f"{'Reference':<12} {'Kind':<10} {'Control':<12} "
      f"{'n':>6} {'mean (Å)':>10} {'std (Å)':>9} "
      f"{'<STRONG_MAX':>12} {'<STABLE_STD':>12}")
print(f"{'─'*70}")

for ref_label, ctrls in controls.items():
    for kind in ('exposed', 'pocket'):
        for res, arr in zip(ctrls[kind], ctrl_data[ref_label][kind]):
            ctrl_label = f"{res.name}{res.resSeq}"
            if len(arr) < 5:
                print(f"{ref_label:<12} {kind:<10} {ctrl_label:<12} "
                      f"{'<5':>6}  (insufficient data)")
                continue
            m = arr.mean()
            s = arr.std()
            in_strong = '✓' if m < STRONG_MAX else '✗'
            in_stable = '✓' if s < STABLE_STD else '✗'
            print(f"{ref_label:<12} {kind:<10} {ctrl_label:<12} "
                  f"{len(arr):>6} {m:>10.3f} {s:>9.3f} "
                  f"{in_strong:>12} {in_stable:>12}")

print(f"{'─'*70}")

# ─────────────────────────────────────────────────────────────────────────────
# 11. SAVE THRESHOLDS TO JSON  (read by water_hbond_stability.py)
# ─────────────────────────────────────────────────────────────────────────────
import json

thresholds = {
    "STRONG_MIN":   0.0,
    "STRONG_MAX":   round(STRONG_MAX, 4),
    "STABLE_STD":   round(STABLE_STD, 4),
    "coverage_pct": COVERAGE_PCT,
    "n_controls":   int(len(all_means)),
    "seq_id":       seq_id,
}
thresh_path = os.path.join(out_dir, f"{seq_id}_thresholds.json")
with open(thresh_path, 'w') as fh:
    json.dump(thresholds, fh, indent=2)

print(f"  STRONG_MIN = {thresholds['STRONG_MIN']}")
print(f"  STRONG_MAX = {thresholds['STRONG_MAX']} Å")
print(f"  STABLE_STD = {thresholds['STABLE_STD']} Å")