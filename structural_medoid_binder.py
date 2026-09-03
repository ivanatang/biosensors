#!/usr/bin/env python
"""
structural_medoid_binder.py — determine the representative (structural medoid) consensus binder.

Stage 1: per-sequence representative frame (frame closest to the windowed average,
         aligned on pocket Calpha), written as a medoid PDB.
Stage 2: pairwise pocket-Calpha RMSD across all consensus-binder medoids; the
         sequence with minimum mean RMSD to all others is the structural medoid.

Usage:
  # one sequence (SLURM-parallelizable)
  python structural_medoid_binder.py config.yaml stage1 --seq_id pair_3062_binder
  # all consensus binders, serially
  python structural_medoid_binder.py config.yaml stage1_all
  # aggregate to pick the medoid
  python structural_medoid_binder.py config.yaml stage2
  # do everything
  python structural_medoid_binder.py config.yaml all
"""
import os, sys, argparse, shutil
import numpy as np
import pandas as pd
import yaml
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

# ── Config / classification helpers ───────────────────────────────────────────
def load_cfg(path):
    """Loads the YAML config file.

    Args:
        path (str): Path to the config YAML.

    Returns:
        dict: Parsed config.
    """
    with open(path) as f:
        return yaml.safe_load(f)

def resolve_paths(cfg):
    """Expands the config's base path variables.

    Args:
        cfg (dict): Loaded config (see load_cfg).

    Returns:
        tuple[str, str, dict]: (base, runrel, type_subdir).
    """
    base   = os.path.expandvars(cfg["paths"]["base"])
    runrel = cfg["paths"]["runrel"]
    tsub   = cfg["paths"]["type_subdir"]
    return base, runrel, tsub

def rundir_for(cfg, seq_id, seq_type):
    """Resolves one sequence's production run directory.

    Args:
        cfg (dict): Loaded config (see load_cfg).
        seq_id (str): Sequence ID.
        seq_type (str): Key into cfg["paths"]["type_subdir"].

    Returns:
        str: Run directory path, or the config's per-seq_id override if set.
    """
    base, runrel, tsub = resolve_paths(cfg)
    overrides = cfg["paths"].get("overrides", {})
    if seq_id in overrides:
        return os.path.expandvars(overrides[seq_id])
    return os.path.join(base, tsub[seq_type], seq_id, runrel)

def consensus_binders(cfg):
    """Finds consensus-binder sequences (Binder group + full pocket motif).

    Args:
        cfg (dict): Loaded config (see load_cfg); reads cfg["medoid"]["binder_motif"]
            and cfg["paths"]["feat_table"].

    Returns:
        list[str]: Folder names of consensus binders, excluding "open"
        conformer variants.
    """
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
    # exclude any "open" conformer variants by convention
    cb = cb[~cb["name"].str.contains("open", case=False)]
    return cb["name"].tolist()

# ── Stage 1: per-sequence representative frame ────────────────────────────────
def stage1_one(cfg, seq_id, seq_type="binder"):
    """Writes one sequence's medoid frame (closest to the windowed average).

    Aligns the equilibration-trimmed, strided trajectory on pocket
    C-alpha, then picks the frame with minimum RMSD to the average
    aligned structure.

    Args:
        cfg (dict): Loaded config (see load_cfg).
        seq_id (str): Sequence ID.
        seq_type (str): Key into cfg["paths"]["type_subdir"] (default:
            "binder").

    Returns:
        str | None: Path to the written medoid PDB, or None if the medoid
        already exists or the trajectory files are missing.
    """
    base, runrel, tsub = resolve_paths(cfg)
    mcfg   = cfg["medoid"]
    rundir = rundir_for(cfg, seq_id, seq_type)
    tpr    = os.path.join(rundir, cfg["trajectory"]["tpr"])
    xtc    = os.path.join(rundir, cfg["trajectory"]["xtc"])

    out_dir = os.path.join(rundir, mcfg["output_subdir"])
    os.makedirs(out_dir, exist_ok=True)
    out_pdb = os.path.join(out_dir, mcfg["output_pdb"])

    if os.path.exists(out_pdb):
        print(f"{seq_id}: medoid PDB already exists, skipping")
        return out_pdb
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        print(f"{seq_id}: MISSING trajectory ({tpr} / {xtc}) — skipping")
        return None

    spacing  = cfg["trajectory"]["frame_spacing_ps"]
    start    = int(mcfg["equil_ns"] * 1000 / spacing)
    stride   = mcfg["stride"]
    ca_sel   = f"name CA and resid {mcfg['pocket_resids']}"

    u = mda.Universe(tpr, xtc)
    # load only the windowed, strided frames into memory (stride controls memory)
    u.transfer_to_memory(start=start, step=stride)

    # align in-memory frames on pocket CA to remove rigid-body motion
    align.AlignTraj(u, u, select=ca_sel, ref_frame=0, in_memory=True).run()

    ca = u.select_atoms(ca_sel)
    coords = np.array([ca.positions.copy() for _ in u.trajectory])   # (nF, nCA, 3)
    avg    = coords.mean(axis=0)
    d2avg  = np.sqrt(((coords - avg) ** 2).sum(axis=2).mean(axis=1))  # per-frame RMSD to avg
    k      = int(d2avg.argmin())

    u.trajectory[k]
    u.select_atoms("protein or resname LIG").write(out_pdb)
    print(f"{seq_id}: medoid frame {k}/{len(coords)} "
          f"(RMSD-to-avg {d2avg[k]:.3f} A) -> {out_pdb}")
    return out_pdb

def stage1_all(cfg):
    """Runs stage1_one for every consensus binder, logging any failures.

    Args:
        cfg (dict): Loaded config (see load_cfg).
    """
    for seq_id in consensus_binders(cfg):
        try:
            stage1_one(cfg, seq_id, "binder")
        except Exception as e:
            print(f"{seq_id}: ERROR in stage1 — {e}")

# ── Stage 2: cross-sequence structural medoid ─────────────────────────────────
def stage2(cfg):
    """Finds the structural medoid across all consensus-binder medoid PDBs.

    Computes pairwise pocket-C-alpha RMSD (optimal superposition) between
    every consensus binder's Stage 1 medoid, picks the sequence with
    minimum mean RMSD to all others, and saves the RMSD matrix, ranking,
    a heatmap, and the chosen reference structure.

    Args:
        cfg (dict): Loaded config (see load_cfg).

    Returns:
        str | None: The structural medoid's seq_id, or None if fewer than
        2 medoid PDBs are available.

    Raises:
        ValueError: Pocket-C-alpha atom counts differ across sequences
            (pocket_resids mismatch).
    """
    base   = os.path.expandvars(cfg["paths"]["base"])
    mcfg   = cfg["medoid"]
    ca_sel = f"name CA and resid {mcfg['pocket_resids']}"
    out_dir = os.path.join(base, mcfg["aggregate_outdir"])
    os.makedirs(out_dir, exist_ok=True)

    # gather medoid PDBs
    pos, missing = {}, []
    for seq_id in consensus_binders(cfg):
        rundir = rundir_for(cfg, seq_id, "binder")
        pdb = os.path.join(rundir, mcfg["output_subdir"], mcfg["output_pdb"])
        if not os.path.exists(pdb):
            missing.append(seq_id); continue
        ag = mda.Universe(pdb).select_atoms(ca_sel)
        pos[seq_id] = ag.positions.copy()

    if missing:
        print(f"Missing medoid PDBs ({len(missing)}): {missing[:5]}"
              + (f" ... +{len(missing)-5}" if len(missing) > 5 else ""))

    sids = sorted(pos)
    n = len(sids)
    if n < 2:
        print("Need >=2 medoid PDBs for Stage 2."); return

    # verify identical atom counts
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
    src = os.path.join(rundir_for(cfg, medoid_seq, "binder"),
                       mcfg["output_subdir"], mcfg["output_pdb"])
    ref_out = os.path.join(out_dir, "reference_structure.pdb")
    shutil.copyfile(src, ref_out)

    # heatmap
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(max(6, n * 0.3), max(5, n * 0.3)))
        im = ax.imshow(D, cmap="viridis")
        ax.set_xticks(range(n)); ax.set_xticklabels(sids, rotation=90, fontsize=5)
        ax.set_yticks(range(n)); ax.set_yticklabels(sids, fontsize=5)
        plt.colorbar(im, label="Pocket Calpha RMSD (A)")
        ax.set_title("Pairwise pocket-Calpha RMSD (consensus binders)")
        plt.tight_layout()
        fig.savefig(os.path.join(out_dir, "pairwise_rmsd_heatmap.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"(heatmap skipped: {e})")

    # ── Report + optional sanity check vs existing metrics ────────────────────
    print("\n=== Structural medoid ranking (lowest mean RMSD = most central) ===")
    print(rank.head(10).to_string(index=False))
    print(f"\nStructural medoid: {medoid_seq}")
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
def main():
    """Dispatches to stage1/stage1_all/stage2 per the `mode` CLI argument."""
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("mode", choices=["stage1", "stage1_all", "stage2", "all"])
    ap.add_argument("--seq_id")
    ap.add_argument("--seq_type", default="binder")
    args = ap.parse_args()
    cfg = load_cfg(args.config)

    if args.mode == "stage1":
        if not args.seq_id:
            sys.exit("stage1 requires --seq_id")
        stage1_one(cfg, args.seq_id, args.seq_type)
    elif args.mode == "stage1_all":
        stage1_all(cfg)
    elif args.mode == "stage2":
        stage2(cfg)
    elif args.mode == "all":
        stage1_all(cfg)
        stage2(cfg)

if __name__ == "__main__":
    main()
