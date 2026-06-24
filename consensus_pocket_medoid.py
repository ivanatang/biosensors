#!/usr/bin/env python
"""
consensus_pocket_medoid.py — identify the most structurally central consensus binder.

For each consensus-binder sequence, reads the pre-existing medoid_PL.pdb from its
run directory. Computes pairwise pocket-Calpha RMSD across all sequences and selects
the structural medoid: the sequence with minimum mean RMSD to all others.

Outputs:
  pairwise_pocket_ca_rmsd.csv  — full N×N RMSD matrix
  medoid_ranking.csv           — sequences ranked by mean RMSD to others
  pairwise_rmsd_heatmap.png    — visual heatmap of the matrix
  reference_structure.pdb      — the chosen medoid's medoid_PL.pdb, copied here

Usage:
  python consensus_pocket_medoid.py config.yaml
"""
import os, sys, shutil
import numpy as np
import pandas as pd
import yaml
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd

# ── Config helpers ─────────────────────────────────────────────────────────────
def load_cfg(path):
    with open(path) as f:
        return yaml.safe_load(f)

def resolve_paths(cfg):
    base   = os.path.expandvars(cfg["paths"]["base"])
    runrel = cfg["paths"]["runrel"]
    tsub   = cfg["paths"]["type_subdir"]
    return base, runrel, tsub

def rundir_for(cfg, seq_id, seq_type):
    base, runrel, tsub = resolve_paths(cfg)
    overrides = cfg["paths"].get("overrides", {})
    if seq_id in overrides:
        return os.path.expandvars(overrides[seq_id])
    return os.path.join(base, tsub[seq_type], seq_id, runrel)

def consensus_binders(cfg):
    """Return list of consensus-binder folder names (binder + full motif)."""
    base = os.path.expandvars(cfg["paths"]["base"])
    feat = cfg["paths"].get("feat_table", os.path.join(base, "feat_table.xlsx"))
    motif = {int(k): v for k, v in cfg["medoid"]["binder_motif"].items()}
    n_motif = len(motif)

    ft = pd.read_excel(feat, sheet_name="all_feats_500ns")
    ft.columns = ft.columns.str.strip()
    ft = ft.rename(columns={c: c.lower() for c in ft.columns
                            if c.lower() in ["name", "group", "sequence"]})
    ft["group"]    = ft["group"].astype(str).str.strip().str.lower()
    ft["sequence"] = ft["sequence"].astype(str).str.strip()
    ft = ft[ft["group"] == "binder"].reset_index(drop=True)

    def score(seq):
        return sum(seq[p - 1].upper() == aa for p, aa in motif.items())
    ft["motif_score"] = ft["sequence"].apply(score)
    cb = ft[ft["motif_score"] == n_motif]
    cb = cb[~cb["name"].str.contains("open", case=False)]
    return cb["name"].tolist()

# ── Cross-sequence pocket medoid ───────────────────────────────────────────────
def find_group_pocket_medoid(cfg):
    base   = os.path.expandvars(cfg["paths"]["base"])
    mcfg   = cfg["medoid"]
    ca_sel = f"name CA and resid {mcfg['pocket_resids']}"
    out_dir = os.path.join(base, mcfg["aggregate_outdir"])
    os.makedirs(out_dir, exist_ok=True)

    # gather pocket-CA coordinates from each sequence's medoid_PL.pdb
    pos, missing = {}, []
    for seq_id in consensus_binders(cfg):
        rundir = rundir_for(cfg, seq_id, "binder")
        pdb = os.path.join(rundir, "medoid_PL.pdb")
        if not os.path.exists(pdb):
            missing.append(seq_id); continue
        ag = mda.Universe(pdb).select_atoms(ca_sel)
        pos[seq_id] = ag.positions.copy()

    if missing:
        print(f"Missing medoid_PL.pdb ({len(missing)}): {missing[:5]}"
              + (f" ... +{len(missing)-5}" if len(missing) > 5 else ""))

    sids = sorted(pos)
    n = len(sids)
    if n < 2:
        sys.exit("Need >=2 medoid_PL.pdb files to compute a group medoid.")

    # verify identical pocket-CA atom counts across all sequences
    counts = {len(v) for v in pos.values()}
    if len(counts) != 1:
        raise ValueError(f"Pocket-CA atom counts differ across sequences: {counts}. "
                         "Check pocket_resids vs each structure.")

    # pairwise pocket-CA RMSD with optimal superposition
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = rmsd(pos[sids[i]], pos[sids[j]],
                                     center=True, superposition=True)

    mean_to_others = D.sum(axis=1) / (n - 1)
    order = np.argsort(mean_to_others)
    medoid_seq = sids[order[0]]

    # ── Save outputs ──────────────────────────────────────────────────────────
    Dmat = pd.DataFrame(D, index=sids, columns=sids)
    Dmat.to_csv(os.path.join(out_dir, "pairwise_pocket_ca_rmsd.csv"))

    rank = pd.DataFrame({
        "seq_id": [sids[k] for k in order],
        "mean_rmsd_to_others_A": np.round(mean_to_others[order], 4),
    })
    rank.to_csv(os.path.join(out_dir, "medoid_ranking.csv"), index=False)

    # copy chosen medoid to a stable reference path
    src = os.path.join(rundir_for(cfg, medoid_seq, "binder"), "medoid_PL.pdb")
    ref_out = os.path.join(out_dir, "reference_structure.pdb")
    shutil.copyfile(src, ref_out)

    # heatmap
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(max(6, n * 0.3), max(5, n * 0.3)))
        im = ax.imshow(D, cmap="viridis")
        ax.set_xticks(range(n)); ax.set_xticklabels(sids, rotation=90, fontsize=5)
        ax.set_yticks(range(n)); ax.set_yticklabels(sids, fontsize=5)
        plt.colorbar(im, label="Pocket Calpha RMSD (Å)")
        ax.set_title("Pairwise pocket-Cα RMSD — consensus binders")
        plt.tight_layout()
        fig.savefig(os.path.join(out_dir, "pairwise_rmsd_heatmap.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"(heatmap skipped: {e})")

    # ── Report + optional sanity check vs existing metrics ────────────────────
    print("\n=== Pocket medoid ranking (lowest mean RMSD = most central) ===")
    print(rank.head(10).to_string(index=False))
    print(f"\nGroup pocket medoid: {medoid_seq}")
    print(f"Reference structure written to: {ref_out}")

    metrics_csv = mcfg.get("metrics_csv")
    if metrics_csv and os.path.exists(metrics_csv):
        m = pd.read_csv(metrics_csv)
        key = "folder_name" if "folder_name" in m.columns else "seq_id"
        m = m[m[key].isin(sids)]
        row = m[m[key] == medoid_seq]
        if not row.empty:
            print("\n=== Sanity check: medoid percentile within consensus binders ===")
            for col in m.select_dtypes("number").columns:
                pct = (m[col] < row[col].values[0]).mean() * 100
                print(f"  {col}: value={row[col].values[0]:.3f}  (~{pct:.0f}th percentile)")

    return medoid_seq

# ── CLI ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Find the group pocket medoid across all consensus-binder sequences."
    )
    ap.add_argument("config", help="Path to config YAML")
    args = ap.parse_args()
    find_group_pocket_medoid(load_cfg(args.config))
