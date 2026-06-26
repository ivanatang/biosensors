#!/usr/bin/env python3
"""
Water-mediated H-bond stability analysis
=========================================
Implements the R-score and H-bond stability methodology from:
  Leonard et al., ACS Chem. Biol. 2024, 19, 1757-1772.

Key quantities
--------------
R-score  = (D - W) * I
  D  = per-residue occupancy of *direct*        heavy-atom contacts  (<4 Å)
  W  = per-residue occupancy of *water-mediated* heavy-atom contacts  (<4 Å water bridge)
  I  = occupancy of either contact type (union, not sum — avoids double-counting
       frames where both contact types co-occur simultaneously)
  → R = +1  : purely direct;   R = -1 : purely water-mediated
  → R =  0  : equal direct and water-mediated (D == W), or no contact at all
  Threshold: R < -0.7 → "dominant" water-mediated interaction (Figure 5A/B)

H-bond stability (for dominant residues, Figure 5 C/D)
  "Strong" : STRONG_MIN <= mean H-acceptor distance <= STRONG_MAX
  "Stable" : std  H-acceptor distance  < STABLE_STD
  Thresholds are derived from control residue distributions
  (hbond_threshold_calibration.py, Figure S4 workflow).
  Both criteria must hold for a water-mediated H-bond to be
  classified as functionally significant (Methods, paper).

Usage  (two-step workflow matching the paper)
-----
Step 1 — R-scores:
    Set seq_id and paths in CONFIG, then:
        conda activate IS_env && python water_hbond_stability.py
    This produces {seq_id}_R_scores.csv identifying dominant residues.

Step 2 — Calibrate thresholds (Figure S4):
        python hbond_threshold_calibration.py
    This finds control residues, plots KDE distributions, and writes
    {seq_id}_thresholds.json.

Step 3 — H-bond stability classification (Figure 5C/D):
    Re-run water_hbond_stability.py (it reads the thresholds JSON automatically)
    and produces the final stability results and plots.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from scipy.stats import gaussian_kde
from collections import defaultdict
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
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
base   = "/scratch/alpine/ivta1597/LCA_boltz_models"
ext    = "HMR/dodecahedron"
prod   = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path    = os.path.join(base, seq_type, seq_id, prod, "prod_md_500ns.xtc")
top_path     = os.path.join(base, seq_type, seq_id, prod, "prod_md_500ns.gro")
out_dir   = os.path.join(base, seq_type, seq_id, f"water_contacts_{TAG}")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"          # residue name of your ligand in the .gro/.pdb
WATER_RESNAMES = {"HOH", "WAT", "SOL"}
ION_RESNAMES   = {"NA", "CL", "NA+", "CL-"}

HEAVY_CUT    = 0.40    # nm  (= 4.0 Å)  heavy-atom contact threshold
R_DOM        = -0.70   # R-score threshold for dominant water-mediated

thresh_path  = os.path.join(base, seq_type, seq_id, f"water_contacts_{TAG}", f"{seq_id}_thresholds_{TAG}.json")
rscores_path = os.path.join(base, seq_type, seq_id, f"water_contacts_{TAG}", f"{seq_id}_R_scores_{TAG}.csv")

STRIDE = 10    # stride for main analysis; use 1 for publication quality

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD TRAJECTORY
# ─────────────────────────────────────────────────────────────────────────────
def load_trajectory(traj_path, top_path, stride):
    print(f"Loading trajectory (stride={stride})…")
    traj = md.load(traj_path, top=top_path, stride=stride)
    print(f"  {traj.n_frames} frames  |  {traj.n_atoms} atoms")
    return traj


# ─────────────────────────────────────────────────────────────────────────────
# 2. ATOM-INDEX HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def heavy(indices, top):
    """Return subset of indices that are heavy (non-H) atoms."""
    return [i for i in indices if top.atom(i).element.symbol != 'H']

def acceptors(indices, top):
    """Return subset of indices that are H-bond acceptors (O or N)."""
    return [i for i in indices if top.atom(i).element.symbol in ('O', 'N')]

def parse_topology(top, LIG_RESNAME, WATER_RESNAMES, ION_RESNAMES):
    """Extract ligand, water, and protein atom index arrays."""
    all_residues = list(top.residues)

    lig_res   = [r for r in all_residues if r.name == LIG_RESNAME]
    wat_res   = [r for r in all_residues if r.name in WATER_RESNAMES]
    skip      = WATER_RESNAMES | ION_RESNAMES | {LIG_RESNAME}
    prot_res  = [r for r in all_residues if r.name not in skip]

    if not lig_res:
        raise ValueError(f"No residue named '{LIG_RESNAME}' found in topology.")

    # Ligand
    lig_all_idx  = [a.index for r in lig_res for a in r.atoms]
    lig_heavy_idx = heavy(lig_all_idx, top)
    lig_acc_idx   = acceptors(lig_all_idx, top)

    # Waters: parallel arrays for O and (H1, H2)
    wat_O_idx  = []
    wat_H_idx  = []           # list of tuples (H1, H2)
    for r in wat_res:
        O = [a.index for a in r.atoms if a.element.symbol == 'O']
        H = [a.index for a in r.atoms if a.element.symbol == 'H']
        if O and len(H) == 2:
            wat_O_idx.append(O[0])
            wat_H_idx.append((H[0], H[1]))
    wat_O_idx = np.array(wat_O_idx, dtype=int)

    # Protein residues: store heavy indices and acceptor indices
    prot_heavy_by_res = {}
    prot_acc_by_res   = {}
    for r in prot_res:
        all_idx = [a.index for a in r.atoms]
        hv  = heavy(all_idx, top)
        acc = acceptors(all_idx, top)
        if hv:
            prot_heavy_by_res[r.index] = np.array(hv, dtype=int)
            prot_acc_by_res[r.index]   = np.array(acc, dtype=int)

    print(f"  Ligand  : {len(lig_heavy_idx)} heavy atoms")
    print(f"  Waters  : {len(wat_O_idx)}")
    print(f"  Protein : {len(prot_heavy_by_res)} residues with heavy atoms")

    return (lig_heavy_idx, lig_acc_idx,
            wat_O_idx, wat_H_idx,
            prot_res, prot_heavy_by_res, prot_acc_by_res)


# ─────────────────────────────────────────────────────────────────────────────
# 3. COMPUTE DIRECT AND WATER-MEDIATED CONTACT OCCUPANCIES
# ─────────────────────────────────────────────────────────────────────────────
def compute_contact_occupancies(traj, lig_heavy_idx, wat_O_idx,
                                prot_res, prot_heavy_by_res,
                                HEAVY_CUT):
    """
    For each protein residue, compute per-frame boolean arrays for:
      direct[res_idx] → bool (F,) : direct heavy-atom contact with ligand
      wmed[res_idx]   → bool (F,) : water-mediated contact with ligand

    Strategy (vectorised per frame):
      1. Find waters whose O is within HEAVY_CUT of ANY ligand heavy atom.
      2. For each protein residue, check:
           direct : any res heavy atom within HEAVY_CUT of any lig heavy atom
           wmed   : any of the above bridging waters within HEAVY_CUT of any
                    res heavy atom
    """
    nf = traj.n_frames
    n_lig = len(lig_heavy_idx)

    # Pre-extract coordinate arrays  (n_frames, n_atoms, 3)
    lig_xyz  = traj.xyz[:, lig_heavy_idx, :]   # (F, L, 3)
    watO_xyz = traj.xyz[:, wat_O_idx,     :]   # (F, W, 3)

    direct_occ = {r.index: np.zeros(nf, dtype=bool) for r in prot_res}
    wmed_occ   = {r.index: np.zeros(nf, dtype=bool) for r in prot_res}

    # Collect residue arrays for batch processing
    res_indices  = [r.index for r in prot_res if r.index in prot_heavy_by_res]
    res_heavy    = [prot_heavy_by_res[ri] for ri in res_indices]

    print("Computing contact occupancies…")
    for f in range(nf):
        if f % 200 == 0:
            print(f"  frame {f:>5d}/{nf}", end='\r', flush=True)

        lp = lig_xyz[f]    # (L, 3)
        wp = watO_xyz[f]   # (W, 3)

        # ── Waters near ligand ──────────────────────────────────────────────
        # Broadcast: (W, L, 3) → (W, L) distances
        dWL = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)  # (W, L)
        near_lig_mask = dWL.min(axis=1) < HEAVY_CUT                      # (W,)
        near_lig_wpos = wp[near_lig_mask]                                 # (W', 3)
        has_bridging  = near_lig_wpos.shape[0] > 0

        # ── Per-residue contacts ────────────────────────────────────────────
        for ri, rh in zip(res_indices, res_heavy):
            rp = traj.xyz[f, rh, :]  # (R, 3)

            # Direct contact
            dRL = np.linalg.norm(rp[:, None, :] - lp[None, :, :], axis=-1)  # (R, L)
            direct_occ[ri][f] = dRL.min() < HEAVY_CUT

            # Water-mediated: bridging water also within HEAVY_CUT of residue
            if has_bridging:
                dRW = np.linalg.norm(rp[:, None, :] - near_lig_wpos[None, :, :], axis=-1)  # (R, W')
                wmed_occ[ri][f] = dRW.min() < HEAVY_CUT

    print(f"\n  Done.")
    return direct_occ, wmed_occ


# ─────────────────────────────────────────────────────────────────────────────
# 3b. COMPUTE WATER-MEDIATED OCCUPANCY FOR SPECIFIC RESIDUES ONLY
# ─────────────────────────────────────────────────────────────────────────────
def compute_wmed_occ_for_residues(traj, lig_heavy_idx, wat_O_idx,
                                   target_res_indices, prot_heavy_by_res,
                                   HEAVY_CUT):
    """Compute water-mediated contact booleans only for a specific list of
    residue topology indices.  Used when dominant residues are already known
    from R_scores.csv, avoiding a full scan of all residues and skipping
    direct contact computation since only wmed_occ is needed downstream.
    """
    nf       = traj.n_frames
    lig_xyz  = traj.xyz[:, lig_heavy_idx, :]   # (F, L, 3)
    watO_xyz = traj.xyz[:, wat_O_idx,     :]   # (F, W, 3)

    wmed_occ  = {ri: np.zeros(nf, dtype=bool) for ri in target_res_indices}
    res_heavy = [(ri, prot_heavy_by_res[ri]) for ri in target_res_indices
                 if ri in prot_heavy_by_res]

    print(f"Computing wmed_occ for {len(res_heavy)} dominant residues...")
    for f in range(nf):
        if f % 200 == 0:
            print(f'  frame {f:>5d}/{nf}', end='\r', flush=True)

        lp = lig_xyz[f]
        wp = watO_xyz[f]

        dWL           = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)
        near_lig_wpos = wp[dWL.min(axis=1) < HEAVY_CUT]
        if near_lig_wpos.shape[0] == 0:
            continue

        for ri, rh in res_heavy:
            rp  = traj.xyz[f, rh, :]
            dRW = np.linalg.norm(rp[:, None, :] - near_lig_wpos[None, :, :], axis=-1)
            wmed_occ[ri][f] = dRW.min() < HEAVY_CUT

    print("\n  Done.")
    return wmed_occ


# ─────────────────────────────────────────────────────────────────────────────
# 4. COMPUTE R-SCORES
# ─────────────────────────────────────────────────────────────────────────────
def compute_R_scores(prot_res, direct_occ, wmed_occ):
    """
    R = (D - W) * I

    As stated in the paper (Leonard et al. ACS Chem. Biol. 2024):
      D = occupancy of direct nonbonded interactions (heavy atom < 4 Å)
      W = occupancy of water-mediated nonbonded interactions
      I = total occupancy of either direct OR water-mediated (union, not sum)

    R returns 0 when D == W because the numerator (D - W) = 0 directly,
    which is the correct behaviour for balanced direct/water-mediated contact.
    I is the boolean union across frames to avoid double-counting frames
    where both contact types are simultaneously present.
    """
    records = []
    for r in prot_res:
        ri = r.index
        D_f = direct_occ[ri].astype(float)
        W_f = wmed_occ[ri].astype(float)
        I_f = (direct_occ[ri] | wmed_occ[ri]).astype(float)

        D = D_f.mean()
        W = W_f.mean()
        I = I_f.mean()

        R = (D - W) * I if (D + W) > 1e-9 else np.nan

        records.append(dict(
            res_index=ri, chain=r.chain.index,
            resSeq=r.resSeq, resname=r.name,
            D=round(D, 4), W=round(W, 4), I=round(I, 4),
            R=round(R, 4) if not np.isnan(R) else np.nan,
        ))

    df = pd.DataFrame(records)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 5a. LOAD THRESHOLDS  (derived by hbond_threshold_calibration.py / Figure S4)
# ─────────────────────────────────────────────────────────────────────────────
def load_thresholds(thresh_path):
    """
    Read STRONG_MIN, STRONG_MAX, STABLE_STD from the JSON file written by
    hbond_threshold_calibration.py.

    Raises a descriptive RuntimeError if the file is missing, so the user
    knows they need to run the calibration script first — matching the paper's
    workflow where Figure S4 (control residue distributions) is generated
    before the water-contact classifications in Figure 5C/D.
    """
    import json
    if not os.path.exists(thresh_path):
        raise RuntimeError(
            f"Threshold file not found:\n  {thresh_path}\n\n"
            "Run hbond_threshold_calibration.py first to generate Figure S4 "
            "and derive thresholds from control residue H-bond distributions.\n"
            "That script writes the required JSON automatically."
        )
    with open(thresh_path) as fh:
        d = json.load(fh)
    print(f"  Loaded thresholds from {thresh_path}")
    print(f"    STRONG_MIN = {d['STRONG_MIN']} Å")
    print(f"    STRONG_MAX = {d['STRONG_MAX']} Å  (from {d['n_controls']} control residues)")
    print(f"    STABLE_STD = {d['STABLE_STD']} Å  ({d['coverage_pct']}th percentile)")
    return float(d['STRONG_MIN']), float(d['STRONG_MAX']), float(d['STABLE_STD'])


# ─────────────────────────────────────────────────────────────────────────────
# 5b. H-ACCEPTOR DISTANCE DISTRIBUTIONS (dominant residues)
# ─────────────────────────────────────────────────────────────────────────────
def compute_hbond_stability(traj, df_dominant, wmed_occ,
                            lig_acc_idx,
                            wat_O_idx, wat_H_idx,
                            prot_heavy_by_res, prot_acc_by_res,
                            lig_heavy_idx,
                            HEAVY_CUT,
                            STRONG_MIN, STRONG_MAX, STABLE_STD):
    """
    For each dominant residue (R < R_DOM), compute H-acceptor distance
    distributions separately for the ligand and the protein residue,
    matching Figure 5C/D of the paper which shows two curves per panel.

    For each bridging frame:
      ha_dists_lig : min H-acceptor distance from bridging water H atoms
                     to any acceptor atom on the LIGAND
      ha_dists_res : min H-acceptor distance from bridging water H atoms
                     to any acceptor atom on the PROTEIN RESIDUE

    This separation reveals which end of the water bridge is the stronger
    and more stable interaction — information lost when the two are combined
    into a single minimum.

    Classification uses the combined (ligand+residue) minimum, consistent
    with the paper's Methods description, but both curves are stored for
    plotting.

    Classification:
      Strong : STRONG_MIN <= mean H-acceptor dist <= STRONG_MAX
      Stable : std  H-acceptor dist < STABLE_STD
      Significant : strong AND stable
    """
    nf = traj.n_frames
    lig_xyz  = traj.xyz[:, lig_heavy_idx, :]   # (F, L, 3)
    watO_xyz = traj.xyz[:, wat_O_idx,     :]   # (F, W, 3)

    H_by_O = {wat_O_idx[i]: wat_H_idx[i] for i in range(len(wat_O_idx))}

    results = {}

    for _, row in df_dominant.iterrows():
        ri     = int(row['res_index'])
        label  = f"{row['resname']}{row['resSeq']}"
        rh_idx = prot_heavy_by_res.get(ri, np.array([], dtype=int))
        ra_idx = prot_acc_by_res.get(ri, np.array([], dtype=int))

        ha_dists_lig = []   # water-H → ligand acceptor
        ha_dists_res = []   # water-H → residue acceptor
        ha_dists_min = []   # min of both (used for classification)

        for f in range(nf):
            if not wmed_occ[ri][f]:
                continue

            lp = lig_xyz[f]   # (L, 3)
            wp = watO_xyz[f]  # (W, 3)
            rp = traj.xyz[f, rh_idx, :] if len(rh_idx) else np.empty((0, 3))

            # ── Find bridging waters ────────────────────────────────────────
            dWL = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)
            near_lig = dWL.min(axis=1) < HEAVY_CUT

            if len(rp) == 0 or near_lig.sum() == 0:
                continue

            dWR = np.linalg.norm(
                wp[near_lig][:, None, :] - rp[None, :, :], axis=-1)
            near_res = dWR.min(axis=1) < HEAVY_CUT

            bridging_global_O = wat_O_idx[near_lig][near_res]
            if len(bridging_global_O) == 0:
                continue

            # ── H-acceptor distances — ligand and residue separately ────────
            min_lig = np.inf
            min_res = np.inf

            for O_idx in bridging_global_O:
                H_pair = H_by_O.get(O_idx)
                if H_pair is None:
                    continue
                for H_idx in H_pair:
                    H_pos = traj.xyz[f, H_idx, :]

                    # Distance to each ligand acceptor
                    for acc_idx in lig_acc_idx:
                        d = np.linalg.norm(H_pos - traj.xyz[f, acc_idx, :])
                        if d < min_lig:
                            min_lig = d

                    # Distance to each residue acceptor
                    for acc_idx in ra_idx:
                        d = np.linalg.norm(H_pos - traj.xyz[f, acc_idx, :])
                        if d < min_res:
                            min_res = d

            if min_lig < np.inf:
                ha_dists_lig.append(min_lig * 10.0)   # nm → Å
            if min_res < np.inf:
                ha_dists_res.append(min_res * 10.0)
            # Combined minimum for classification — consistent with paper Methods
            combined = min(min_lig, min_res)
            if combined < np.inf:
                ha_dists_min.append(combined * 10.0)

        results[label] = {
            'res_index':    ri,
            'resname':      row['resname'],
            'resSeq':       row['resSeq'],
            'R':            row['R'],
            'ha_dists':     np.array(ha_dists_min),   # kept for compatibility
            'ha_dists_lig': np.array(ha_dists_lig),
            'ha_dists_res': np.array(ha_dists_res),
        }

    # Summarise — classification uses combined minimum
    summary = []
    for label, d in results.items():
        arr = d['ha_dists']
        if len(arr) == 0:
            print(f"  {label}: no bridging frames — check ligand/water resnames")
            continue
        mean_d = arr.mean()
        std_d  = arr.std()
        strong = STRONG_MIN <= mean_d <= STRONG_MAX
        stable = std_d < STABLE_STD
        sig    = strong and stable

        # Also report per-component means for context
        mean_lig = d['ha_dists_lig'].mean() if len(d['ha_dists_lig']) else float('nan')
        mean_res = d['ha_dists_res'].mean() if len(d['ha_dists_res']) else float('nan')
        print(f"  {label:12s}  n={len(arr):5d}  "
              f"mean(combined)={mean_d:.2f} Å  std={std_d:.2f} Å  "
              f"mean(lig)={mean_lig:.2f} Å  mean(res)={mean_res:.2f} Å  "
              f"{'STRONG' if strong else 'weak  '}  "
              f"{'STABLE' if stable else 'unstable'}  "
              f"→ {'★ SIGNIFICANT' if sig else ''}")
        summary.append(dict(
            residue=label, **{k: d[k] for k in ('resname', 'resSeq', 'R')},
            n_bridging_frames=len(arr),
            mean_HA_dist_A        =round(mean_d,   3),
            std_HA_dist_A         =round(std_d,    3),
            mean_HA_dist_ligand_A =round(mean_lig, 3),
            mean_HA_dist_residue_A=round(mean_res, 3),
            strong=strong, stable=stable, significant=sig,
        ))
        d['mean']   = mean_d
        d['std']    = std_d
        d['strong'] = strong
        d['stable'] = stable

    return results, pd.DataFrame(summary)


# ─────────────────────────────────────────────────────────────────────────────
# 6. PLOTTING
# ─────────────────────────────────────────────────────────────────────────────
def plot_R_scores(df, seq_id, out_dir, R_DOM):
    """Bar chart of R-scores for all residues with non-trivial contact (I > 0.02)."""
    df_plot = df[df['I'] > 0.02].dropna(subset=['R']).sort_values('resSeq')
    if df_plot.empty:
        print("No residues with contact occupancy > 0.02 — skipping R-score plot.")
        return

    fig, ax = plt.subplots(figsize=(max(10, len(df_plot) * 0.45), 4),
                           constrained_layout=True)
    colors = ['#d62728' if r < R_DOM else '#4472C4' for r in df_plot['R']]
    labels = [f"{row['resname']}{row['resSeq']}" for _, row in df_plot.iterrows()]

    ax.bar(range(len(df_plot)), df_plot['R'], color=colors,
           edgecolor='k', linewidth=0.4)
    ax.axhline(R_DOM, color='k', linestyle='--', linewidth=1.0,
               label=f'R = {R_DOM} (dominant water-mediated)')
    ax.axhline(0, color='grey', linewidth=0.6, alpha=0.5)

    ax.set_xticks(range(len(df_plot)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel('R-score', fontsize=10)
    ax.set_xlabel('Residue', fontsize=10)
    ax.set_title(f'{seq_id} — Water-mediated contact R-scores', fontsize=11)
    ax.set_ylim(-1.15, 1.15)
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=8)

    out = os.path.join(out_dir, f"{seq_id}_R_scores_{TAG}.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Saved → {out}")


def plot_HA_distributions(results, seq_id, out_dir,
                          STRONG_MIN, STRONG_MAX, STABLE_STD):
    """
    Replicates Figure 5C/D from the paper.

    Each dominant residue gets one panel showing TWO curves:
      - ligand    (pink)  : H-acceptor distance from bridging water H
                            to acceptors on the LIGAND
      - residue   (teal)  : H-acceptor distance from bridging water H
                            to acceptors on the PROTEIN RESIDUE

    This matches the paper's Figure 5C/D layout and makes clear which end
    of the water bridge is the stronger/more stable interaction.
    The green shaded region marks the strong H-bond range [STRONG_MIN, STRONG_MAX].
    The orange dotted line marks STRONG_MAX + STABLE_STD as a visual stability guide.
    """
    dom = {k: v for k, v in results.items() if len(v['ha_dists']) > 0}
    if not dom:
        print("No dominant residues with data — skipping H-bond distribution plot.")
        return

    COLOR_LIG = '#d62728'   # red   — ligand
    COLOR_RES = '#1f77b4'   # blue  — protein residue

    n = len(dom)
    fig, axes = plt.subplots(1, n, figsize=(3.8 * n, 3.8), constrained_layout=True)
    if n == 1:
        axes = [axes]

    # Fixed bandwidth of 0.1 Å — Scott's rule over-smooths the sharp peaks
    # typical of genuine H-bond distance distributions (~1.8–2.5 Å range)
    KDE_BW = 0.10   # Å

    for ax, (label, d) in zip(axes, dom.items()):

        def _kde_curve(arr, color, curve_label):
            """KDE with fixed bandwidth, plotted over the data range."""
            if len(arr) < 5:
                return
            try:
                # gaussian_kde takes bandwidth as fraction of std;
                # convert fixed Å value accordingly
                bw = KDE_BW / arr.std() if arr.std() > 0 else 'scott'
                kde = gaussian_kde(arr, bw_method=bw)
            except Exception:
                return
            x = np.linspace(0, 8, 1000)
            y = kde(x)
            ax.plot(x, y, color=color, linewidth=1.5, label=curve_label)
            ax.fill_between(x, y, alpha=0.15, color=color)
            mean_d = arr.mean()
            ax.axvline(mean_d, color=color, linestyle='--',
                       linewidth=0.9, alpha=0.8,
                       label=f'Mean={mean_d:.2f} Å')

        _kde_curve(d['ha_dists_lig'], COLOR_LIG, f'water → LCA')
        _kde_curve(d['ha_dists_res'], COLOR_RES, f'water → {d["resname"]}{d["resSeq"]}')

        # Set x-axis to cover the actual data range with a small margin
        all_vals = np.concatenate([
            d['ha_dists_lig'] if len(d['ha_dists_lig']) else np.array([]),
            d['ha_dists_res'] if len(d['ha_dists_res']) else np.array([]),
        ])
        if len(all_vals):
            x_lo = max(0, all_vals.min() - 0.5)
            x_hi = min(8, all_vals.max() + 0.5)
            ax.set_xlim(x_lo, x_hi)

        ax.set_title(label, fontsize=9,
                     fontweight='bold' if d.get('significant') else 'normal')
        ax.set_xlabel('H–Acceptor Distance (Å)', fontsize=8)
        ax.set_ylabel('Density', fontsize=8)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.4)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc='upper right')

    fig.suptitle(
        f'{seq_id} — H-acceptor distance distributions\n'
        f'(dominant water-mediated residues, R < {R_DOM})\n'
        f'Red = water→{LIG_RESNAME}   |   Blue = water→protein residue',
        fontsize=9)

    out = os.path.join(out_dir, f"{seq_id}_HA_distributions_{TAG}.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Saved → {out}")


def plot_pocket_overview(df, seq_id, out_dir):
    """
    2D scatter of D vs W occupancy, coloured by R-score — gives a quick
    overview of the direct/water contact balance across the binding pocket.
    """
    df_c = df[df['I'] > 0.02].dropna(subset=['R'])
    if df_c.empty:
        return

    fig, ax = plt.subplots(figsize=(5, 4.5), constrained_layout=True)
    sc = ax.scatter(df_c['D'], df_c['W'], c=df_c['R'],
                    cmap='RdBu', vmin=-1, vmax=1,
                    s=60, edgecolors='k', linewidths=0.4)
    plt.colorbar(sc, ax=ax, label='R-score')

    # Annotate dominant residues
    for _, row in df_c[df_c['R'] < -0.7].iterrows():
        ax.annotate(f"{row['resname']}{row['resSeq']}",
                    xy=(row['D'], row['W']),
                    xytext=(4, 2), textcoords='offset points', fontsize=7)

    ax.set_xlabel('Direct contact occupancy (D)', fontsize=10)
    ax.set_ylabel('Water-mediated contact occupancy (W)', fontsize=10)
    ax.set_title(f'{seq_id} — Direct vs water-mediated contacts', fontsize=10)
    ax.grid(True, alpha=0.4)

    out = os.path.join(out_dir, f"{seq_id}_DW_scatter_{TAG}.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Saved → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # 1. Read R-scores from compute_r_scores.py output (Step 1)
    if not os.path.exists(rscores_path):
        raise FileNotFoundError(
            f"R-scores CSV not found:\n  {rscores_path}\n\n"
            "Run compute_r_scores.py first (Step 1) to generate this file.")

    df_R   = pd.read_csv(rscores_path)
    # Run H-acceptor analysis on all key pocket residues
    # regardless of R-score, to get the full water network picture
    POCKET_RESIDUES = [58, 59, 62, 83, 87, 88, 89, 92, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164]

    df_dom = df_R[df_R['resSeq'].isin(POCKET_RESIDUES)].dropna(subset=['R'])
    print(f"Analysing {len(df_dom)} pocket residues for H-acceptor distributions")
    print(df_dom[['resname', 'resSeq', 'chain', 'D', 'W', 'I', 'R']].to_string(index=False))

    if df_dom.empty:
        print("No dominant residues found - nothing to analyse. Exiting.")
        raise SystemExit(0)

    # 2. Load trajectory and parse topology
    traj = load_trajectory(traj_path, top_path, STRIDE)

    # ── Restrict to analysis window ────────────────────────────────
    mask = (traj.time >= START_PS) & (traj.time <= END_PS)
    traj = traj[mask]
    print(f"  Time window: {START_NS:.0f}–{END_NS:.0f} ns  "
          f"({mask.sum()} of {mask.shape[0]} frames retained after striding)")
    # ──────────────────────────────────────────────────────────────

    top  = traj.topology

    (lig_heavy_idx, lig_acc_idx,
     wat_O_idx, wat_H_idx,
     prot_res, prot_heavy_by_res,
     prot_acc_by_res) = parse_topology(top, LIG_RESNAME,
                                       WATER_RESNAMES, ION_RESNAMES)

    # 3. Compute wmed_occ only for dominant residues (not all residues)
    target_indices = list(df_dom['res_index'].astype(int))
    wmed_occ = compute_wmed_occ_for_residues(
        traj, lig_heavy_idx, wat_O_idx,
        target_indices, prot_heavy_by_res, HEAVY_CUT)

    # 4. Load thresholds from hbond_threshold_calibration.py (Step 2)
    print("\nLoading H-bond thresholds...")
    STRONG_MIN, STRONG_MAX, STABLE_STD = load_thresholds(thresh_path)

    # 5. H-bond stability classification
    print("\nComputing H-acceptor distance distributions...")
    results, df_summary = compute_hbond_stability(
        traj, df_dom, wmed_occ,
        lig_acc_idx, wat_O_idx, wat_H_idx,
        prot_heavy_by_res, prot_acc_by_res, lig_heavy_idx,
        HEAVY_CUT, STRONG_MIN, STRONG_MAX, STABLE_STD)

    if not df_summary.empty:
        df_summary.to_csv(
            os.path.join(out_dir, f"{seq_id}_hbond_summary_{TAG}.csv"), index=False)
        print(f"\nH-bond summary saved -> {out_dir}/{seq_id}_hbond_summary_{TAG}.csv")

    # 6. Plot
    print("\nGenerating plots...")
    plot_HA_distributions(results, seq_id, out_dir,
                          STRONG_MIN, STRONG_MAX, STABLE_STD)

    print("\n=== Done ===")
    print("Significant water-mediated H-bonds (strong AND stable):")
    if not df_summary.empty:
        sig = df_summary[df_summary['significant'] == True]
        print(sig[['residue', 'mean_HA_dist_A', 'std_HA_dist_A', 'R']
                   ].to_string(index=False) if not sig.empty else "  None found.")
