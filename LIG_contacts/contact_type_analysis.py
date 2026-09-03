"""Computes protein-ligand contact-type features for one PYR1 sequence.

Computes hydrophobic/polar/charged contact-type features for a single
sequence trajectory on Alpine HPC. seq_id is always passed as a positional
CLI argument, either directly for local testing or by
contact_type_worker.sh for SLURM array runs.

Usage:
    python contact_type_analysis.py <seq_id>

Output:
    <output_dir>/<seq_id>_contact_perframe.csv: per-frame contact counts.
    <output_dir>/<seq_id>_residue_occupancy.csv: per-residue occupancy.
    <output_dir>/<seq_id>_contact_summary.csv: scalar features for feat_table.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

# ─────────────────────────────────────────────
# PATHS  — edit these for your environment
# ─────────────────────────────────────────────
# Trajectory inputs are read from the PetaLibrary archive, not scratch —
# scratch auto-deletes after 90 days and older runs' xtc/gro are already gone.
# Layout mirrors scratch's LCA_boltz_models: binders/nonbinders/neg_low_pkt/neg_fail_gate.
base    = "/pl/active/shirts_archive/IvanaTang/biosensors"
runrel  = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

# Results are written to the persistent repo location (not scratch), so
# they survive scratch's 90-day auto-deletion.
output_base = "/projects/ivta1597/biosensors/LIG_contacts"
output_dir  = os.path.join(output_base, "contact_type_results")   # default; overridden by TAG in main()

# ─────────────────────────────────────────────
# TYPE SUBDIRECTORY MAP
# ─────────────────────────────────────────────
TYPE_SUBDIR = {
    "binder":    "binders",
    "nb":        "nonbinders",
    "low_pkt":   "neg_low_pkt",
    "fail_gate": "neg_fail_gate",
}

def get_type_subdir(seq_id):
    """Maps a seq_id's naming suffix to its data subdirectory.

    Checks two-token suffixes (low_pkt, fail_gate) before single-token ones
    (binder, nb) to avoid partial matches.

    Examples:
        pair_3059_binder -> binders
        pair_0777_nb -> nonbinders
        pair_0008_low_pkt -> neg_low_pkt
        pair_0715_fail_gate -> neg_fail_gate

    Args:
        seq_id (str): Sequence identifier.

    Returns:
        str: Subdirectory name (a value from TYPE_SUBDIR).

    Raises:
        ValueError: `seq_id` doesn't end with a recognized suffix.
    """
    tokens = seq_id.split("_")

    # Check two-token suffix first
    if len(tokens) >= 2:
        two_token = "_".join(tokens[-2:])
        if two_token in TYPE_SUBDIR:
            return TYPE_SUBDIR[two_token]

    # Fall back to single-token suffix
    if tokens[-1] in TYPE_SUBDIR:
        return TYPE_SUBDIR[tokens[-1]]

    raise ValueError(
        f"Cannot determine type subdirectory from seq_id '{seq_id}'. "
        f"Expected suffix from: {list(TYPE_SUBDIR.keys())}"
    )


def run_dir(seq_id, end_ns):
    """Finds the directory containing a sequence's medoid_PL.pdb / xtc.

    Newer pipeline runs write these under runrel/{end_ns}ns/; older ones
    write directly into runrel/. Prefers whichever location actually has
    medoid_PL.pdb rather than guessing from seq_id (see the equivalent
    rmsf_run_dir() in extract_rmsf_feats.py for why guessing is unreliable).

    Args:
        seq_id (str): Sequence identifier.
        end_ns (float): End of the analysis window in ns, used to build the
            nested subdirectory name.

    Returns:
        str: Path to the run directory (nested if medoid_PL.pdb is there,
        otherwise flat).
    """
    flat_dir   = os.path.join(base, get_type_subdir(seq_id), seq_id, runrel)
    nested_dir = os.path.join(flat_dir, f"{int(end_ns)}ns")
    if os.path.exists(os.path.join(nested_dir, "medoid_PL.pdb")):
        return nested_dir
    return flat_dir

# ─────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────
LIGAND_RESNAME  = "LIG"
CUTOFF_NM       = 0.45     # 4.5 Å heavy-atom contact cutoff
STRIDE          = 1        # set >1 only for quick debugging; use 1 for publication

# Ligand heavy-atom split between the C20-C24 pentanoic-acid tail (side
# chain off the steroid D-ring, ending in the carboxylate) and everything
# else (fused 4-ring core, angular methyls, ring hydroxyls/sulfate).
#
# Expressed as 0-indexed ORDINAL POSITIONS in topology order, not atom
# names: the production .itp/.gro use generic per-element names ("O", "O",
# "C", ...), not the unique names ("O41", "C44", ...) in the standalone
# ligand_<seq_id>.pdb used to derive this split, so name-based matching
# would silently select the wrong atoms. Atom order is preserved between a
# sequence's standalone PDB and its .gro, but NOT across Boltz-2 prediction
# batches (the original 95-seq LCA cohort and the 4-sequence new-ligand
# test batch order the same LCA molecule differently), and CDCA/GLCA/LCA3S
# are chemically distinct ligands (extra ring hydroxyl, glycine conjugate,
# sulfate ester in place of the free 3-OH, respectively).
#
# Each pattern was derived once from that ligand's own standalone PDB (with
# CONECT records): prune degree-1 leaves to isolate the fused steroid-ring
# core, then take the largest connected component off it as the tail (the
# solvent-facing side chain; for GLCA this extends through the amide and
# glycine conjugate, the analogous terminus once LCA's free acid becomes an
# amide). get_tail_ordinals() checks this at runtime so a sequence matching
# no known pattern fails loudly instead of silently mis-splitting.
LIGAND_PATTERNS = [
    # (label, expected_elements, tail_ordinals)
    ('LCA (original 95-seq cohort)',
     ['O', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
     {0, 1, 3, 4, 7, 9, 19}),
    ('LCA (new-ligand test batch)',
     ['C', 'C', 'C', 'C', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
     {0, 1, 2, 3, 4, 5, 6}),
    ('CDCA (new-ligand test batch)',
     ['C', 'C', 'C', 'C', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
     {0, 1, 2, 3, 4, 5, 6}),
    ('GLCA (new-ligand test batch, tail includes glycine conjugate)',
     ['C', 'C', 'C', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
    ('LCA3S (new-ligand test batch)',
     ['C', 'C', 'C', 'C', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'S', 'O', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
     {0, 1, 2, 3, 4, 5, 6}),
]


def get_tail_ordinals(seq_id, got_elements):
    """Looks up the tail-atom ordinal positions for a ligand's element order.

    Args:
        seq_id (str): Sequence ID, used only in the error message.
        got_elements (list[str]): Heavy-atom element symbols for this
            ligand, in topology order.

    Returns:
        set[int]: Ordinal positions (0-indexed) of the tail's heavy atoms.

    Raises:
        ValueError: `got_elements` doesn't match any pattern in
            LIGAND_PATTERNS (a new ligand or a new Boltz-2 atom ordering).
    """
    for label, expected, tail_ordinals in LIGAND_PATTERNS:
        if got_elements == expected:
            return tail_ordinals
    known = "\n".join(f"  {label}: {pat}" for label, pat, _ in LIGAND_PATTERNS)
    raise ValueError(
        f"[{seq_id}] Ligand heavy-atom element order doesn't match any known "
        f"pattern (likely a new ligand, or a new Boltz-2 batch with yet "
        f"another atom ordering -- derive its own pattern and add it to "
        f"LIGAND_PATTERNS).\nGot: {got_elements}\nKnown patterns:\n{known}"
    )

# ─────────────────────────────────────────────
# RESIDUE CHEMISTRY MAP
# ─────────────────────────────────────────────
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR       = {"SER", "THR", "CYS", "TYR", "ASN", "GLN"}
POS_CHARGED = {"ARG", "LYS", "HIS"}
NEG_CHARGED = {"ASP", "GLU"}

def classify_residue(resname):
    """Classifies a residue's side chain as hydrophobic, polar, or charged.

    Args:
        resname (str): Three-letter residue name.

    Returns:
        str: "hydrophobic", "polar", "pos_charged", "neg_charged", or
        "other" if unrecognized.
    """
    if resname in HYDROPHOBIC:  return "hydrophobic"
    if resname in POLAR:        return "polar"
    if resname in POS_CHARGED:  return "pos_charged"
    if resname in NEG_CHARGED:  return "neg_charged"
    return "other"


# ─────────────────────────────────────────────
# ARGUMENTS
# ─────────────────────────────────────────────
def parse_args():
    """Parses CLI arguments for a single-sequence contact-type run.

    Returns:
        argparse.Namespace: Parsed arguments (seq_id, start_ns, end_ns,
        ligand_region, suffix).
    """
    parser = argparse.ArgumentParser(
        description="Compute protein–ligand contact type features for one sequence."
    )
    parser.add_argument("seq_id",
                        help="Sequence identifier (e.g. pair_3059_binder)")
    parser.add_argument("--start-ns", type=float, default=40.0,
                        help="Start of analysis window in ns (default: 40)")
    parser.add_argument("--end-ns",   type=float, default=500.0,
                        help="End of analysis window in ns (default: 500)")
    parser.add_argument("--ligand-region", choices=["whole", "core", "tail"], default="whole",
                        help="Restrict ligand atoms to the steroid ring core, the "
                             "C20-C24 carboxylate tail, or the whole ligand "
                             "(default: whole)")
    parser.add_argument("--suffix", default="",
                        help="Run-directory/output suffix, e.g. '_qfix' for the "
                             "bond-order/charge-fix systems (default: '', the "
                             "standard production directory)")
    return parser.parse_args()


# ─────────────────────────────────────────────
# LOAD TRAJECTORY
# ─────────────────────────────────────────────
def load_trajectory(seq_id, start_ps, end_ps, end_ns):
    """Loads and time-windows a sequence's medoid trajectory.

    Args:
        seq_id (str): Sequence identifier.
        start_ps (int): Start of the analysis window in ps.
        end_ps (int): End of the analysis window in ps.
        end_ns (float): End of the analysis window in ns (used to locate
            the run directory; see run_dir).

    Returns:
        mdtraj.Trajectory: Trajectory restricted to [start_ps, end_ps].

    Raises:
        FileNotFoundError: The expected topology or trajectory file is
            missing.
    """
    seq_dir  = run_dir(seq_id, end_ns)
    top_path = os.path.join(seq_dir, "medoid_PL.pdb")
    xtc_path = os.path.join(seq_dir, f"PL_only_40_{int(end_ns)}ns.xtc")

    for path in (top_path, xtc_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Expected file not found: {path}")

    print(f"[{seq_id}] Loading topology: {top_path}")
    print(f"[{seq_id}] Loading trajectory: {xtc_path}  (stride={STRIDE})")

    traj = md.load(xtc_path, top=top_path, stride=STRIDE)
    print(f"[{seq_id}] Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")

    # ── Restrict to analysis window ────────────────────────────────────
    mask = (traj.time >= start_ps) & (traj.time <= end_ps)
    traj = traj[mask]
    print(f"[{seq_id}] Time window: {start_ps/1000:.0f}–{end_ps/1000:.0f} ns  "
          f"({mask.sum()} frames retained)")
    # ──────────────────────────────────────────────────────────────────

    return traj


# ─────────────────────────────────────────────
# CONTACT TYPE COMPUTATION
# ─────────────────────────────────────────────
def compute_contact_type_features(traj, seq_id, ligand_region="whole"):
    """Computes per-frame and per-residue ligand contact-type features.

    For each protein residue, tests per-frame heavy-atom contact
    (CUTOFF_NM) with the ligand and tallies it by residue chemistry class.

    Args:
        traj (mdtraj.Trajectory): Trajectory to analyze.
        seq_id (str): Sequence identifier, used in error messages.
        ligand_region (str): "whole", "core", or "tail" (default: "whole");
            see get_tail_ordinals.

    Returns:
        tuple[pandas.DataFrame, pandas.DataFrame]: (per_frame_df, with one
        row per frame and a column per contact-type count; residue_df, with
        one row per residue and its contact occupancy).

    Raises:
        ValueError: No atoms found with residue name LIGAND_RESNAME.
    """
    top = traj.topology

    # ── Ligand heavy atoms, in topology order, optionally restricted to
    # core/tail by ordinal position (see LIGAND_PATTERNS above) ────────────
    lig_heavy_atoms_all = [
        a for a in top.atoms
        if a.residue.name == LIGAND_RESNAME and a.element.symbol != "H"
    ]
    if not lig_heavy_atoms_all:
        raise ValueError(
            f"[{seq_id}] No atoms found with residue name '{LIGAND_RESNAME}'. "
            f"Check LIGAND_RESNAME. Available residue names: "
            f"{sorted({r.name for r in top.residues})}"
        )

    if ligand_region != "whole":
        got_elements = [a.element.symbol for a in lig_heavy_atoms_all]
        tail_ordinals = get_tail_ordinals(seq_id, got_elements)
        if ligand_region == "core":
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i not in tail_ordinals]
        elif ligand_region == "tail":
            lig_heavy_atoms_all = [a for i, a in enumerate(lig_heavy_atoms_all)
                                    if i in tail_ordinals]

    lig_atoms = [a.index for a in lig_heavy_atoms_all]
    print(f"[{seq_id}] Ligand '{LIGAND_RESNAME}': {len(lig_atoms)} heavy atoms "
          f"(region={ligand_region})")

    # ── Protein residues ─────────────────────────────────────────────────
    prot_residues = [
        (r.index, r.name, classify_residue(r.name))
        for r in top.residues if r.is_protein
    ]
    print(f"[{seq_id}] Protein residues: {len(prot_residues)}")

    # ── Per-residue contact calculation ──────────────────────────────────
    records = []
    for res_idx, res_name, res_class in prot_residues:
        prot_atoms = [
            a.index for a in top.residue(res_idx).atoms
            if a.element.symbol != "H"
        ]
        if not prot_atoms:
            continue

        # All heavy-atom pairs between this residue and LCA
        pairs = np.array(
            [[pa, la] for pa in prot_atoms for la in lig_atoms],
            dtype=int
        )

        # Minimum distance to LCA per frame  →  shape (n_frames,)
        dists    = md.compute_distances(traj, pairs)   # (n_frames, n_pairs)
        min_dist = dists.min(axis=1)
        in_contact = (min_dist < CUTOFF_NM).astype(np.int8)

        records.append({
            "res_idx":    res_idx,
            "res_name":   res_name,
            "res_class":  res_class,
            "occupancy":  in_contact.mean(),
            "in_contact": in_contact,        # (n_frames,) — kept for summing
        })

    residue_df = pd.DataFrame([
        {k: v for k, v in r.items() if k != "in_contact"}
        for r in records
    ])

    # ── Per-frame type counts ─────────────────────────────────────────────
    contact_matrix = np.vstack([r["in_contact"] for r in records])  # (n_res, n_frames)
    type_labels    = np.array([r["res_class"] for r in records])

    contact_types = ["hydrophobic", "polar", "pos_charged", "neg_charged", "other"]
    per_frame = {}
    for t in contact_types:
        mask = (type_labels == t)
        per_frame[f"n_{t}"] = contact_matrix[mask].sum(axis=0).astype(int)

    total = sum(per_frame[f"n_{t}"] for t in contact_types)
    per_frame["n_total"] = total
    per_frame["frac_hydrophobic"] = np.where(
        total > 0, per_frame["n_hydrophobic"] / total, 0.0
    )

    per_frame_df = pd.DataFrame(per_frame)
    per_frame_df.index.name = "frame"

    return per_frame_df, residue_df


# ─────────────────────────────────────────────
# SUMMARISE → SCALAR FEATURES FOR feat_table
# ─────────────────────────────────────────────
def summarise_features(per_frame_df, seq_id):
    """Collapses per-frame contact-type arrays into scalar summary features.

    Args:
        per_frame_df (pandas.DataFrame): Per-frame contact counts (see
            compute_contact_type_features).
        seq_id (str): Sequence identifier, added as a column.

    Returns:
        pandas.DataFrame: One row, with mean/std of each column and, for
        count columns, the fraction of frames with >=1 contact.
    """
    rows = {"seq_id": seq_id}

    numeric_cols = [c for c in per_frame_df.columns
                    if c not in ("frame",)]

    for col in numeric_cols:
        arr = per_frame_df[col].values
        rows[f"mean_{col}"] = arr.mean()
        rows[f"std_{col}"]  = arr.std()
        if col.startswith("n_"):
            rows[f"occ_{col}_gt0"] = (arr > 0).mean()

    return pd.DataFrame([rows])


# ─────────────────────────────────────────────
# DIAGNOSTIC PRINT
# ─────────────────────────────────────────────
def print_diagnostics(per_frame_df, residue_df, seq_id):
    """Prints a summary of contact-type means and top-occupancy residues.

    Args:
        per_frame_df (pandas.DataFrame): Per-frame contact counts.
        residue_df (pandas.DataFrame): Per-residue occupancy.
        seq_id (str): Sequence identifier, used in the printed header.
    """
    print(f"\n── Contact summary for {seq_id} ─────────────────────────")
    print(f"  Frames analysed : {len(per_frame_df)}")
    print(f"  Mean total contacts     : {per_frame_df['n_total'].mean():.2f}")
    print(f"  Mean hydrophobic        : {per_frame_df['n_hydrophobic'].mean():.2f}")
    print(f"  Mean polar              : {per_frame_df['n_polar'].mean():.2f}")
    print(f"  Mean pos_charged        : {per_frame_df['n_pos_charged'].mean():.2f}")
    print(f"  Mean neg_charged        : {per_frame_df['n_neg_charged'].mean():.2f}")
    print(f"  Mean frac_hydrophobic   : {per_frame_df['frac_hydrophobic'].mean():.3f}")

    top_contacts = (
        residue_df[residue_df["occupancy"] > 0.10]
        .sort_values("occupancy", ascending=False)
        .head(10)
    )
    print(f"\n  Top residues by occupancy (>10%):")
    print(top_contacts[["res_name", "res_class", "occupancy"]].to_string(index=False))
    print("─" * 55)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    """Runs the contact-type analysis for one sequence end to end."""
    global runrel

    args     = parse_args()
    seq_id   = args.seq_id
    start_ns = args.start_ns
    end_ns   = args.end_ns
    start_ps = int(start_ns * 1000)
    end_ps   = int(end_ns   * 1000)
    TAG      = f"{int(start_ns)}_{int(end_ns)}ns"
    ligand_region = args.ligand_region

    # "" reproduces the original runrel exactly; e.g. "_qfix" points run_dir()
    # at the bond-order/charge-fix systems' differently-named run directory.
    runrel = f"{runrel}{args.suffix}"

    # Region suffix keeps core/tail runs from overwriting the whole-ligand
    # results; "whole" reproduces the original, unsuffixed path exactly.
    REGION_TAG = "" if ligand_region == "whole" else f"_{ligand_region}"

    # Tagged output directory — one per time window, consistent with
    # water_contacts_{TAG}/ naming used by the water analysis pipeline
    tagged_output_dir = os.path.join(output_base, f"contact_type_results_{TAG}{REGION_TAG}{args.suffix}")
    os.makedirs(tagged_output_dir, exist_ok=True)

    print(f"\n=== contact_type_analysis.py  |  seq_id={seq_id}  "
          f"window={start_ns:.0f}–{end_ns:.0f} ns  region={ligand_region} ===")

    # Check if already done (safe to re-queue without reprocessing)
    summary_out = os.path.join(tagged_output_dir,
                               f"{seq_id}_contact_summary_{TAG}{REGION_TAG}.csv")
    if os.path.exists(summary_out):
        print(f"[{seq_id}] Output already exists, skipping: {summary_out}")
        sys.exit(0)

    traj = load_trajectory(seq_id, start_ps, end_ps, end_ns)

    per_frame_df, residue_df = compute_contact_type_features(traj, seq_id, ligand_region)

    print_diagnostics(per_frame_df, residue_df, seq_id)

    # ── Write outputs ─────────────────────────────────────────────────────
    per_frame_out = os.path.join(tagged_output_dir,
                                 f"{seq_id}_contact_perframe_{TAG}{REGION_TAG}.csv")
    residue_out   = os.path.join(tagged_output_dir,
                                 f"{seq_id}_residue_occupancy_{TAG}{REGION_TAG}.csv")

    per_frame_df.to_csv(per_frame_out)
    residue_df.to_csv(residue_out, index=False)
    summarise_features(per_frame_df, seq_id).to_csv(summary_out, index=False)

    print(f"[{seq_id}] Wrote: {per_frame_out}")
    print(f"[{seq_id}] Wrote: {residue_out}")
    print(f"[{seq_id}] Wrote: {summary_out}")
    print(f"[{seq_id}] Done.")


if __name__ == "__main__":
    main()
