"""
Water-mediated H-bond stability analysis
=========================================
Implements the R-score and H-bond stability methodology from:
  Leonard et al., ACS Chem. Biol. 2024, 19, 1757-1772.

Key quantities
--------------
R-score  = (D - W) / (D + W) * I
  D  = per-residue occupancy of *direct*        heavy-atom contacts  (<4 Å)
  W  = per-residue occupancy of *water-mediated* heavy-atom contacts  (<4 Å water bridge)
  I  = occupancy of either contact type (union)
  → R = +1  : purely direct;   R = -1 : purely water-mediated
  Threshold: R < -0.7 → "dominant" water-mediated interaction (Figure 5 A/B)

H-bond stability (for dominant residues, Figure 5 C/D)
  "Strong" : mean H-acceptor distance  2.0 – 2.5 Å
  "Stable" : std  H-acceptor distance  < 0.45 Å
  Both criteria must hold for a water-mediated H-bond to be
  classified as functionally significant (Methods, paper).

Usage
-----
1. Set the path / sequence variables at the top of the CONFIG block.
2. Run interactively in a Jupyter cell or as a script on Alpine:
       conda activate IS_env && python water_hbond_stability.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
base    = "/scratch/alpine/ivta1597/LCA_boltz_models"
seq_id  = "pair_XXXX"                              # ← set your sequence ID
prod    = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

traj_path = os.path.join(base, seq_id, prod, "prod_md_500ns.xtc")
top_path  = os.path.join(base, seq_id, prod, "prod_md_500ns.gro")
out_dir   = os.path.join(base, seq_id, "water_contacts")
os.makedirs(out_dir, exist_ok=True)

LIG_RESNAME    = "LIG"          # residue name of your ligand in the .gro/.pdb
WATER_RESNAMES = {"HOH", "WAT", "SOL"}
ION_RESNAMES   = {"NA", "CL", "NA+", "CL-"}

HEAVY_CUT    = 0.40    # nm  (= 4.0 Å)  heavy-atom contact threshold
R_DOM        = -0.70   # R-score threshold for dominant water-mediated

# Thresholds below are NOT set manually — they are determined at runtime by
# calibrate_thresholds(), which mirrors the paper's approach:
#   "numerical ranges were set to include 95% of analyzed water-residue H-bonds"
# These variables are populated by calibrate_thresholds() and used downstream.
STRONG_MIN   = None    # Å  2.5th  percentile of per-residue mean H-acceptor dist
STRONG_MAX   = None    # Å  97.5th percentile of per-residue mean H-acceptor dist
STABLE_STD   = None    # Å  95th   percentile of per-residue std  H-acceptor dist

CALIB_STRIDE   = 50    # stride for calibration scan (coarser is fine — large sample)
STRIDE         = 10    # stride for main analysis; use 1 for publication quality


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD TRAJECTORY
# ─────────────────────────────────────────────────────────────────────────────
def load_trajectory(traj_path, top_path, stride):
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
# 4. COMPUTE R-SCORES
# ─────────────────────────────────────────────────────────────────────────────
def compute_R_scores(prot_res, direct_occ, wmed_occ):
    """
    R = (D - W) / (D + W) * I
    where D, W, I are scalar occupancies (mean over frames).
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

        denom = D + W
        R = (D - W) / denom * I if denom > 1e-9 else np.nan

        records.append(dict(
            res_index=ri, chain=r.chain.index,
            resSeq=r.resSeq, resname=r.name,
            D=round(D, 4), W=round(W, 4), I=round(I, 4),
            R=round(R, 4) if not np.isnan(R) else np.nan,
        ))

    df = pd.DataFrame(records)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 5a. CALIBRATE STRONG / STABLE THRESHOLDS FROM THE SIMULATION
# ─────────────────────────────────────────────────────────────────────────────
def calibrate_thresholds(traj, top, wat_O_idx, wat_H_idx,
                         prot_res, prot_heavy_by_res, prot_acc_by_res,
                         HEAVY_CUT, out_dir, seq_id, calib_stride=50):
    """
    Derive STRONG_MIN, STRONG_MAX, and STABLE_STD from the simulation itself,
    following the paper's approach:

        "The numerical ranges for classifying strong and stable H-bonds were
         set to include 95% of analyzed water-residue H-bonds. We included
         H-bonds with solvent-exposed residues as well as residues within the
         binding pocket."    — Leonard et al. Methods

    Algorithm
    ---------
    For every protein residue (pocket AND solvent-exposed):
      1. Find frames where any water O is within HEAVY_CUT of a heavy atom.
      2. Among those waters, compute the minimum H-acceptor distance to any
         acceptor atom on that residue for each such frame.
      3. Record the per-residue mean and std dev of those distances.

    Thresholds are then set from the distribution of per-residue statistics:
      STRONG_MIN = 2.5th  percentile of per-residue means   (central 95%)
      STRONG_MAX = 97.5th percentile of per-residue means
      STABLE_STD = 95th  percentile of per-residue std devs

    Parameters
    ----------
    calib_stride : int
        Sub-sample the trajectory for calibration. A coarser stride (default 50)
        is fine — we need distributional statistics, not per-frame precision.

    Returns
    -------
    STRONG_MIN, STRONG_MAX, STABLE_STD  (all in Å)
    """
    print(f"\nCalibrating H-bond thresholds (calib_stride={calib_stride})…")
    sub = traj[::calib_stride]
    nf  = sub.n_frames

    watO_xyz = sub.xyz[:, wat_O_idx, :]   # (F, W, 3)
    H_by_O   = {wat_O_idx[i]: wat_H_idx[i] for i in range(len(wat_O_idx))}

    per_res_means = []
    per_res_stds  = []

    res_indices = [r.index for r in prot_res if r.index in prot_heavy_by_res]
    n_res = len(res_indices)

    for k, ri in enumerate(res_indices):
        if k % 200 == 0:
            print(f"  residue {k:>4d}/{n_res}", end='\r', flush=True)

        rh_idx = prot_heavy_by_res[ri]
        ra_idx = prot_acc_by_res.get(ri, np.array([], dtype=int))
        if len(ra_idx) == 0:
            continue

        ha_dists_res = []

        for f in range(nf):
            rp = sub.xyz[f, rh_idx, :]   # (R, 3)
            wp = watO_xyz[f]             # (W, 3)

            # Waters within HEAVY_CUT of any residue heavy atom
            dWR = np.linalg.norm(
                wp[:, None, :] - rp[None, :, :], axis=-1)  # (W, R)
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
                    H_pos = sub.xyz[f, H_idx, :]
                    for acc_idx in ra_idx:
                        d = np.linalg.norm(H_pos - sub.xyz[f, acc_idx, :])
                        if d < min_d:
                            min_d = d

            if min_d < np.inf:
                ha_dists_res.append(min_d * 10.0)   # nm → Å

        if len(ha_dists_res) < 5:   # skip residues with barely any water contact
            continue

        arr = np.array(ha_dists_res)
        per_res_means.append(arr.mean())
        per_res_stds.append(arr.std())

    print(f"\n  Collected statistics from {len(per_res_means)} residues.")

    if len(per_res_means) < 10:
        raise RuntimeError(
            "Too few residues with water contacts for reliable calibration.\n"
            "Check LIG_RESNAME, WATER_RESNAMES, or lower calib_stride.")

    means_arr = np.array(per_res_means)
    stds_arr  = np.array(per_res_stds)

    strong_min = float(np.percentile(means_arr, 2.5))
    strong_max = float(np.percentile(means_arr, 97.5))
    stable_std = float(np.percentile(stds_arr,  95.0))

    print(f"\n  Calibrated thresholds (95% coverage over {len(per_res_means)} residues):")
    print(f"    STRONG_MIN = {strong_min:.3f} Å  (2.5th  pctile of per-residue means)")
    print(f"    STRONG_MAX = {strong_max:.3f} Å  (97.5th pctile of per-residue means)")
    print(f"    STABLE_STD = {stable_std:.3f} Å  (95th   pctile of per-residue std devs)")

    # ── Diagnostic plot (mirrors Figure S4 in the paper) ─────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), constrained_layout=True)

    axes[0].hist(means_arr, bins=40, color='#4472C4', edgecolor='k',
                 linewidth=0.3, alpha=0.8)
    axes[0].axvline(strong_min, color='#d62728', linestyle='--', linewidth=1.2,
                    label=f'2.5th pctile  = {strong_min:.2f} Å')
    axes[0].axvline(strong_max, color='#d62728', linestyle='-.', linewidth=1.2,
                    label=f'97.5th pctile = {strong_max:.2f} Å')
    axes[0].set_xlabel('Per-residue mean H-acceptor dist (Å)', fontsize=9)
    axes[0].set_ylabel('Count', fontsize=9)
    axes[0].set_title('Strong H-bond calibration (central 95%)', fontsize=9)
    axes[0].legend(fontsize=7)
    axes[0].grid(True, alpha=0.4)

    axes[1].hist(stds_arr, bins=40, color='#70AD47', edgecolor='k',
                 linewidth=0.3, alpha=0.8)
    axes[1].axvline(stable_std, color='#d62728', linestyle='--', linewidth=1.2,
                    label=f'95th pctile = {stable_std:.2f} Å')
    axes[1].set_xlabel('Per-residue std H-acceptor dist (Å)', fontsize=9)
    axes[1].set_ylabel('Count', fontsize=9)
    axes[1].set_title('Stable H-bond calibration (95th pctile)', fontsize=9)
    axes[1].legend(fontsize=7)
    axes[1].grid(True, alpha=0.4)

    fig.suptitle(f'{seq_id} — H-bond threshold calibration\n'
                 f'(all water–residue contacts, {len(per_res_means)} residues)',
                 fontsize=10)
    calib_out = os.path.join(out_dir, f"{seq_id}_hbond_calibration.png")
    fig.savefig(calib_out, dpi=150)
    plt.show()
    print(f"  Calibration plot → {calib_out}")

    return strong_min, strong_max, stable_std


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
    For each residue with R < R_DOM:
      - Restrict to frames where the water-mediated contact is present
        (uses wmed_occ from compute_contact_occupancies — no closure needed)
      - Find the bridging water(s) in each such frame
      - Compute the minimum H-acceptor distance across all
        (water-H, acceptor-on-ligand-or-residue) pairs
      - Classify using data-derived thresholds from calibrate_thresholds()

    Classification:
      Strong : STRONG_MIN <= mean H-acceptor dist <= STRONG_MAX
      Stable : std  H-acceptor dist < STABLE_STD
      Significant : strong AND stable
    """
    nf = traj.n_frames
    lig_xyz  = traj.xyz[:, lig_heavy_idx, :]   # (F, L, 3)
    watO_xyz = traj.xyz[:, wat_O_idx,     :]   # (F, W, 3)

    # Build O→H lookup  {O_global_idx: (H1_global, H2_global)}
    H_by_O = {wat_O_idx[i]: wat_H_idx[i] for i in range(len(wat_O_idx))}

    results = {}

    for _, row in df_dominant.iterrows():
        ri     = int(row['res_index'])
        label  = f"{row['resname']}{row['resSeq']}"
        rh_idx = prot_heavy_by_res.get(ri, np.array([], dtype=int))
        ra_idx = prot_acc_by_res.get(ri, np.array([], dtype=int))
        all_acc_idx = list(lig_acc_idx) + list(ra_idx)

        ha_dists = []

        for f in range(nf):
            # Only consider frames where water-mediated contact was recorded
            if not wmed_occ[ri][f]:
                continue

            lp = lig_xyz[f]   # (L, 3)
            wp = watO_xyz[f]  # (W, 3)
            rp = traj.xyz[f, rh_idx, :] if len(rh_idx) else np.empty((0, 3))

            # ── Find bridging waters ────────────────────────────────────────
            dWL = np.linalg.norm(wp[:, None, :] - lp[None, :, :], axis=-1)
            near_lig = dWL.min(axis=1) < HEAVY_CUT  # (W,)

            if len(rp) == 0 or near_lig.sum() == 0:
                continue

            dWR = np.linalg.norm(
                wp[near_lig][:, None, :] - rp[None, :, :], axis=-1)  # (W', R)
            near_res = dWR.min(axis=1) < HEAVY_CUT                    # (W',)

            bridging_global_O = wat_O_idx[near_lig][near_res]
            if len(bridging_global_O) == 0:
                continue

            # ── Min H-acceptor distance across all bridging waters ──────────
            min_d = np.inf
            for O_idx in bridging_global_O:
                H_pair = H_by_O.get(O_idx)
                if H_pair is None:
                    continue
                for H_idx in H_pair:
                    H_pos = traj.xyz[f, H_idx, :]     # (3,)
                    for acc_idx in all_acc_idx:
                        acc_pos = traj.xyz[f, acc_idx, :]
                        d = np.linalg.norm(H_pos - acc_pos)
                        if d < min_d:
                            min_d = d

            if min_d < np.inf:
                ha_dists.append(min_d * 10.0)   # nm → Å

        results[label] = {
            'res_index': ri,
            'resname':   row['resname'],
            'resSeq':    row['resSeq'],
            'R':         row['R'],
            'ha_dists':  np.array(ha_dists),
        }

    # Summarise
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
        print(f"  {label:12s}  n={len(arr):5d}  mean={mean_d:.2f} Å  "
              f"std={std_d:.2f} Å  "
              f"{'STRONG' if strong else 'weak  '}  "
              f"{'STABLE' if stable else 'unstable'}  "
              f"→ {'★ SIGNIFICANT' if sig else ''}")
        summary.append(dict(
            residue=label, **{k: d[k] for k in ('resname', 'resSeq', 'R')},
            n_bridging_frames=len(arr),
            mean_HA_dist_A=round(mean_d, 3),
            std_HA_dist_A =round(std_d, 3),
            strong=strong, stable=stable, significant=sig,
        ))
        d['mean'] = mean_d
        d['std']  = std_d
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

    out = os.path.join(out_dir, f"{seq_id}_R_scores.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Saved → {out}")


def plot_HA_distributions(results, seq_id, out_dir,
                          STRONG_MIN, STRONG_MAX, STABLE_STD):
    """
    Replicates Figure 5C/D from the paper:
    Probability distribution of H-acceptor distances for each dominant residue.
    """
    dom = {k: v for k, v in results.items() if len(v['ha_dists']) > 0}
    if not dom:
        print("No dominant residues with data — skipping H-bond distribution plot.")
        return

    n = len(dom)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.5), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, (label, d) in zip(axes, dom.items()):
        arr   = d['ha_dists']
        color = '#1f77b4' if d.get('significant') else '#aec7e8'

        counts, edges = np.histogram(arr, bins=40, range=(0, 6))
        probs = counts / counts.max() if counts.max() > 0 else counts
        centers = (edges[:-1] + edges[1:]) / 2

        ax.plot(centers, probs, color=color, linewidth=1.5, label=seq_id)
        ax.axvspan(STRONG_MIN, STRONG_MAX, alpha=0.12, color='green',
                   label='Strong range')
        ax.axvline(STRONG_MAX + STABLE_STD, color='orange', linestyle=':',
                   linewidth=0.8, alpha=0.7)

        mean_d = d['ha_dists'].mean()
        ax.axvline(mean_d, color='k', linestyle='--', linewidth=0.9,
                   label=f'Mean={mean_d:.2f} Å')

        ax.set_title(label, fontsize=9)
        ax.set_xlabel('H–Acceptor Distance (Å)', fontsize=8)
        ax.set_ylabel('Probability', fontsize=8)
        ax.set_xlim(0, 6)
        ax.grid(True, alpha=0.4)
        ax.tick_params(labelsize=7)
        if ax == axes[0]:
            ax.legend(fontsize=6)

        # Annotate strong/stable
        tag = []
        if d.get('strong'):  tag.append('Strong')
        if d.get('stable'):  tag.append('Stable')
        if tag:
            ax.text(0.97, 0.95, '/'.join(tag), transform=ax.transAxes,
                    ha='right', va='top', fontsize=7, color='darkgreen',
                    fontweight='bold')

    fig.suptitle(f'{seq_id} — H-acceptor distance distributions\n'
                 f'(dominant water-mediated residues, R < {-0.7})',
                 fontsize=10)

    out = os.path.join(out_dir, f"{seq_id}_HA_distributions.png")
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

    out = os.path.join(out_dir, f"{seq_id}_DW_scatter.png")
    fig.savefig(out, dpi=150)
    plt.show()
    print(f"Saved → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # 1. Load
    traj = load_trajectory(traj_path, top_path, STRIDE)
    top  = traj.topology

    # 2. Parse topology
    (lig_heavy_idx, lig_acc_idx,
     wat_O_idx, wat_H_idx,
     prot_res, prot_heavy_by_res,
     prot_acc_by_res) = parse_topology(top, LIG_RESNAME,
                                       WATER_RESNAMES, ION_RESNAMES)

    # 3. Contact occupancies
    direct_occ, wmed_occ = compute_contact_occupancies(
        traj, lig_heavy_idx, wat_O_idx,
        prot_res, prot_heavy_by_res, HEAVY_CUT)

    # 4. R-scores
    df_R = compute_R_scores(prot_res, direct_occ, wmed_occ)
    df_R.to_csv(os.path.join(out_dir, f"{seq_id}_R_scores.csv"), index=False)
    print(f"\nR-score table saved → {out_dir}/{seq_id}_R_scores.csv")

    df_dom = df_R[df_R['R'] < R_DOM].dropna(subset=['R'])
    print(f"\nDominant water-mediated residues (R < {R_DOM}):")
    print(df_dom[['resname', 'resSeq', 'chain', 'D', 'W', 'I', 'R']].to_string(index=False))

    # 5a. Calibrate thresholds from this simulation's water–residue H-bonds
    STRONG_MIN, STRONG_MAX, STABLE_STD = calibrate_thresholds(
        traj, top, wat_O_idx, wat_H_idx,
        prot_res, prot_heavy_by_res, prot_acc_by_res,
        HEAVY_CUT, out_dir, seq_id,
        calib_stride=CALIB_STRIDE)

    # 5b. H-bond stability for dominant residues
    print("\nComputing H-acceptor distance distributions…")
    results, df_summary = compute_hbond_stability(
        traj, df_dom, wmed_occ,
        lig_acc_idx, wat_O_idx, wat_H_idx,
        prot_heavy_by_res, prot_acc_by_res, lig_heavy_idx,
        HEAVY_CUT, STRONG_MIN, STRONG_MAX, STABLE_STD)

    if not df_summary.empty:
        df_summary.to_csv(
            os.path.join(out_dir, f"{seq_id}_hbond_summary.csv"), index=False)
        print(f"\nH-bond summary saved → {out_dir}/{seq_id}_hbond_summary.csv")

    # 6. Plots
    print("\nGenerating plots…")
    plot_R_scores(df_R, seq_id, out_dir, R_DOM)
    plot_HA_distributions(results, seq_id, out_dir,
                          STRONG_MIN, STRONG_MAX, STABLE_STD)
    plot_pocket_overview(df_R, seq_id, out_dir)

    print("\n=== Done ===")
    print(f"Significant water-mediated H-bonds (strong AND stable):")
    if not df_summary.empty:
        sig = df_summary[df_summary['significant'] == True]
        print(sig[['residue', 'mean_HA_dist_A', 'std_HA_dist_A', 'R']
                   ].to_string(index=False) if not sig.empty else "  None found.")