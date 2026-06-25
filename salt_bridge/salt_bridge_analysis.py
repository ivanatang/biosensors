#!/usr/bin/env python
"""
salt_bridge_analysis.py
Usage: python salt_bridge_analysis.py <config.yaml> <seq_id> <seq_type> <custom_path>

Scans all basic (chargeable) residue sidechains for proximity to the LCA
carboxylate, with no assumption about which position participates.
Outputs per-residue occupancy across the full trajectory, plus distance
time series for the top-N most-contacted residues.
"""
import sys, os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import yaml

config_path, seq_id, seq_type = sys.argv[1:4]

with open(config_path) as f:
    cfg = yaml.safe_load(f)

# ── Resolve paths ──────────────────────────────────────────────────────────────
base   = os.path.expandvars(cfg["paths"]["base"])
runrel = cfg["paths"]["runrel"]
type_subdir = cfg["paths"]["type_subdir"]

rundir = os.path.join(base, type_subdir[seq_type], seq_id, runrel)

gro = os.path.join(rundir, cfg["trajectory"]["tpr"])
xtc = os.path.join(rundir, cfg["trajectory"]["xtc"])

u = mda.Universe(gro, xtc)

# ── Ligand carboxylate ─────────────────────────────────────────────────────────
lig_resname  = cfg["ligand"]["resname"]
carbox_idx   = cfg["ligand"]["carboxylate_indices"]
hydroxyl_idx = cfg["ligand"]["hydroxyl_index"]

lig      = u.select_atoms(f"resname {lig_resname}")
carbox   = lig[carbox_idx]
hydroxyl = lig[[hydroxyl_idx]]

assert len(carbox) == 2, f"Expected 2 carboxylate atoms, found {len(carbox)}"
assert set(carbox.names) == {"O"}, f"Expected O atoms, got {list(carbox.names)}"

# ── Build candidate atom selection from config (no position assumption) ───────
sb_cfg     = cfg["salt_bridge"]
CUTOFF     = sb_cfg["cutoff_angstrom"]
MIN_OCC    = sb_cfg["min_occupancy_pct"]
N_TOP      = sb_cfg["n_top_for_timeseries"]

selections = []
for entry in sb_cfg["basic_residues"]:
    atom_str = " ".join(entry["atoms"])
    sel = f"resname {entry['resname']} and name {atom_str}"
    selections.append(sel)

candidate_sel = u.select_atoms(" or ".join(f"({s})" for s in selections))
print(f"{seq_id}: scanning {len(candidate_sel)} candidate basic atoms "
      f"across {candidate_sel.n_residues} residues")

if len(candidate_sel) == 0:
    print(f"{seq_id}: no basic residue atoms found — skipping")
    sys.exit(0)

# ── Group candidate atoms by residue for per-residue min-distance ─────────────
residues = list(candidate_sel.residues)
n_res    = len(residues)
res_atom_idx = [candidate_sel.select_atoms(f"resid {r.resid} and resname {r.resname}")
                for r in residues]

# ── Single pass over trajectory: per-residue min distance to carboxylate ──────
n_frames = len(u.trajectory)
dist_matrix = np.full((n_frames, n_res), np.nan)

for ts_i, ts in enumerate(u.trajectory):
    carbox_pos = carbox.positions  # (2, 3)
    for r_i, atoms in enumerate(res_atom_idx):
        if len(atoms) == 0:
            continue
        d = np.linalg.norm(
            carbox_pos[:, None, :] - atoms.positions[None, :, :], axis=-1
        )
        dist_matrix[ts_i, r_i] = d.min()

# ── Per-residue occupancy summary ──────────────────────────────────────────────
rows = []
for r_i, res in enumerate(residues):
    d = dist_matrix[:, r_i]
    if np.all(np.isnan(d)):
        continue
    occ = np.nanmean(d < CUTOFF) * 100
    rows.append({
        "seq_id":        seq_id,
        "resid":         res.resid,
        "resname":       res.resname,
        "mean_dist_A":   round(np.nanmean(d), 4),
        "min_dist_A":    round(np.nanmin(d), 4),
        "occupancy_pct": round(occ, 2),
    })

occ_df = pd.DataFrame(rows).sort_values("occupancy_pct", ascending=False)
occ_df_filtered = occ_df[occ_df["occupancy_pct"] >= MIN_OCC]

print(f"{seq_id}: {len(occ_df_filtered)}/{len(occ_df)} residues "
      f"≥ {MIN_OCC}% occupancy")
if len(occ_df_filtered) > 0:
    print(occ_df_filtered.head(N_TOP).to_string(index=False))

# ── Save occupancy table ───────────────────────────────────────────────────────
out_dir = os.path.join(rundir, sb_cfg["output_subdir"])
os.makedirs(out_dir, exist_ok=True)
occ_df_filtered.to_csv(
    os.path.join(out_dir, sb_cfg["output_files"]["occupancy"]), index=False)

# ── Save distance time series for top-N residues ──────────────────────────────
top_idx = occ_df.head(N_TOP).index.tolist()
ts_df = pd.DataFrame({"frame": np.arange(n_frames)})
for idx in top_idx:
    res = residues[occ_df.index.get_loc(idx)] if idx in occ_df.index else None
for r_i in range(n_res):
    res = residues[r_i]
    label = f"resid{res.resid}_{res.resname}"
    if occ_df.loc[occ_df["resid"] == res.resid, "occupancy_pct"].values[0] >= \
       occ_df.head(N_TOP)["occupancy_pct"].min():
        ts_df[label] = dist_matrix[:, r_i]

ts_df.to_csv(
    os.path.join(out_dir, sb_cfg["output_files"]["timeseries"]), index=False)

print(f"{seq_id}: done. Outputs written to {out_dir}")
